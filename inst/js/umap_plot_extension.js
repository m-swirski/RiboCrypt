(elem, _, data) => {
  if (elem._fullLayout?.dragmode !== "select") {
    Plotly.relayout(elem, { dragmode: "select" });
  }
  const valuesInputId = data;
  const sendSelection = (selection) => {
    Shiny.setInputValue(valuesInputId, selection, { priority: "event" });
  };
  let suppressSelectionEvents = 0;
  const registry = window.__ribocryptUmapSelectionRegistry || {
    widgets: {},
    handlersInstalled: false
  };
  window.__ribocryptUmapSelectionRegistry = registry;

  const releaseProgrammaticSelectionSync = () => {
    const release = () => {
      suppressSelectionEvents = Math.max(0, suppressSelectionEvents - 1);
    };

    if (typeof window.requestAnimationFrame === "function") {
      window.requestAnimationFrame(() => window.requestAnimationFrame(release));
    } else {
      window.setTimeout(release, 50);
    }
  };
  const runProgrammaticSelectionSync = (sync) => {
    suppressSelectionEvents += 1;
    let result = null;
    try {
      result = sync();
    } catch (error) {
      releaseProgrammaticSelectionSync();
      throw error;
    }

    if (result && typeof result.then === "function") {
      result.then(releaseProgrammaticSelectionSync, releaseProgrammaticSelectionSync);
    } else {
      releaseProgrammaticSelectionSync();
    }
  };
  const suppressEvent = (evt) => {
    if (!evt) return;
    if (typeof evt.preventDefault === "function") evt.preventDefault();
    if (typeof evt.stopPropagation === "function") evt.stopPropagation();
    if (typeof evt.stopImmediatePropagation === "function") evt.stopImmediatePropagation();
  };
  const normalizeRunId = (run) => {
    if (run == null) return null;
    const value = String(run).trim();
    return value.length === 0 ? null : value;
  };
  const runIdFromText = (text) => {
    const value = normalizeRunId(text);
    if (!value) return null;
    const parts = value.split(/<br\s*\/?>/i);
    for (let i = 0; i < parts.length; i += 1) {
      const part = parts[i].replace(/<[^>]*>/g, "").trim();
      const match = part.match(/^Run(?: ID)?:\s*(.+)$/i);
      if (match) return normalizeRunId(match[1]);
    }
    // New compact traces carry the Run ID directly; retain compatibility with
    // cached and older traces whose text contains the full hover markup.
    if (!/<[^>]*>|<br\s*\/?>/i.test(value) && !/^Run(?: ID)?:/i.test(value)) return value;
    return null;
  };
  const messagePayload = (message) => {
    if (
      message &&
      typeof message === "object" &&
      !Array.isArray(message) &&
      Object.prototype.hasOwnProperty.call(message, "runs")
    ) {
      return message.runs;
    }
    return message;
  };
  const selectedRunSetFromMessage = (message) => {
    const payload = messagePayload(message);
    const selected = Array.isArray(payload) ? payload : [payload];
    return new Set(selected.map(normalizeRunId).filter((run) => run !== null));
  };

  // Single-trace layout (R/page_organism_UMAP.R): all points are trace 0, the legend entries are empty traces 1..n.
  // layout.meta.rcPointGroup holds each point's 0-based group, so legend clicks hide/show groups through per-point
  // opacity on trace 0 (one WebGL trace keeps box/lasso dragging smooth), and hidden points are never selected.
  const pointGroup = (elem.layout && elem.layout.meta && Array.isArray(elem.layout.meta.rcPointGroup))
    ? elem.layout.meta.rcPointGroup : null;
  const hiddenGroups = new Set();
  const baseOpacity = (elem.data[0] && elem.data[0].marker && typeof elem.data[0].marker.opacity === "number")
    ? elem.data[0].marker.opacity : 0.75;
  const applyGroupVisibility = () => {
    const nGroups = elem.data.length - 1;
    const opacity = pointGroup.map((g) => (hiddenGroups.has(g) ? 0 : baseOpacity));
    const legendVisible = Array.from({ length: nGroups }, (_, i) => (hiddenGroups.has(i) ? "legendonly" : true));
    const legendIdx = Array.from({ length: nGroups }, (_, i) => i + 1);
    Plotly.restyle(elem, { "marker.opacity": [opacity] }, [0]);
    Plotly.restyle(elem, { visible: legendVisible }, legendIdx);
  };
  if (pointGroup) {
    elem.on("plotly_legendclick", (ev) => {
      const g = ev.curveNumber - 1;
      if (g < 0) return false;
      if (hiddenGroups.has(g)) hiddenGroups.delete(g); else hiddenGroups.add(g);
      applyGroupVisibility();
      return false;
    });
    elem.on("plotly_legenddoubleclick", (ev) => {
      const g = ev.curveNumber - 1;
      if (g < 0) return false;
      const nGroups = elem.data.length - 1;
      const isolated = hiddenGroups.size === nGroups - 1 && !hiddenGroups.has(g);
      hiddenGroups.clear();
      if (!isolated) for (let i = 0; i < nGroups; i += 1) if (i !== g) hiddenGroups.add(i);
      applyGroupVisibility();
      return false;
    });
  }
  const pointIsHidden = (point) => pointGroup && point.curveNumber === 0 && hiddenGroups.has(pointGroup[point.pointIndex]);

  const onSelected = (e) => {
    if (suppressSelectionEvents > 0) return;
    const selectedPoints = (e && Array.isArray(e.points) ? e.points : []).filter((point) => !pointIsHidden(point));
    const pointValuesToSet = selectedPoints.map((point) => {
      return runIdFromText(point.text);
    }).filter((run) => run !== null);
    sendSelection(pointValuesToSet);
  };
  const onDeselected = () => {
    if (suppressSelectionEvents > 0) return;
    sendSelection([]);
  };
  const triggerDeselected = (evt) => {
    suppressEvent(evt);
    if (evt && evt.event) suppressEvent(evt.event);
    onDeselected();
    return false;
  };

  // Live selection while dragging: plotly_selecting fires on every mouse move; send at most every 150 ms and only
  // when the selected Run set changed, so the table follows the box. plotly_selected (mouse-up) sends the final set.
  let lastLiveKey = null; let liveTimer = null; let pendingLive = null;
  const flushLive = () => {
    liveTimer = null;
    if (!pendingLive) return;
    const runs = (Array.isArray(pendingLive.points) ? pendingLive.points : []).filter((p) => !pointIsHidden(p))
      .map((p) => runIdFromText(p.text)).filter((r) => r !== null);
    pendingLive = null;
    const key = runs.length + ":" + runs.slice(0, 50).join(",") + ":" + runs.slice(-50).join(",");
    if (key === lastLiveKey) return;
    lastLiveKey = key;
    sendSelection(runs);
  };
  elem.on("plotly_selecting", (e) => {
    if (suppressSelectionEvents > 0) return;
    pendingLive = e;
    if (liveTimer === null) liveTimer = window.setTimeout(flushLive, 150);
  });
  elem.on("plotly_selected", (e) => { lastLiveKey = null; pendingLive = null;
    if (liveTimer !== null) { window.clearTimeout(liveTimer); liveTimer = null; } onSelected(e); });
  elem.on("plotly_deselect", onDeselected);
  elem.on("plotly_doubleclick", triggerDeselected);

  elem.addEventListener("dblclick", triggerDeselected, true);

  const attachInnerDblclick = () => {
    const inner = elem.querySelectorAll(".gl-container, .nsewdrag, .plotly .user-select-none, canvas");
    Array.prototype.forEach.call(inner, (node) => {
      if (!node || node.__rcUmapDblclickBound) return;
      node.__rcUmapDblclickBound = true;
      node.addEventListener("dblclick", triggerDeselected, true);
    });
  };

  attachInnerDblclick();
  elem.on("plotly_afterplot", attachInnerDblclick);

  const onSelectionChanged = (message) => {
    const selectedRuns = selectedRunSetFromMessage(message);
    const tracesToUpdate = Array.from({ length: elem.data.length }, (_, i) => i + 1);
    let indicesToUpdate = Array.from({ length: elem.data.length }, (_, i) => []);

    elem.data.forEach((trace, index) => {
      if (Array.isArray(trace.text)) {
        trace.text.forEach((text, pointIndex) => {
          const run = runIdFromText(text);
          if (run && selectedRuns.has(run)) {
            indicesToUpdate[index].push(pointIndex);
          }
        });
      } else {
        const run = runIdFromText(trace.text);
        if (run && selectedRuns.has(run)) {
          indicesToUpdate[index].push(0);
        }
      }
    });

    let updatedData = [...elem.data]
    tracesToUpdate.forEach((t, index) => {
      updatedData[index]["selectedpoints"] = indicesToUpdate[index]
    });

    let updatedLayout = { ...elem.layout };
    updatedLayout["selections"] = [];

    runProgrammaticSelectionSync(() => Plotly.react(elem, updatedData, updatedLayout));
  };

  const onSelectionReset = (_) => {
    const hasSelection = elem.data.some(trace => trace.selectedpoints != null) ||
      elem._fullData?.some(trace => trace.selectedpoints != null);
    const hasRegion = elem.layout.selections?.length || elem._fullLayout?.selections?.length;
    if (!hasSelection && !hasRegion) return;
    const tracesToUpdate = Array.from({ length: elem.data.length }, (_, i) => i + 1);

    let updatedData = [...elem.data]
    tracesToUpdate.forEach((t, index) => {
      updatedData[index]["selectedpoints"] = null
    });

    let updatedLayout = { ...elem.layout };
    updatedLayout["selections"] = [];

    runProgrammaticSelectionSync(() => Plotly.react(elem, updatedData, updatedLayout));
  };

  registry.widgets[valuesInputId] = {
    elem,
    onSelectionChanged,
    onSelectionReset
  };

  const dispatchSelectionMessage = (message, method) => {
    const target = message && typeof message === "object" && !Array.isArray(message)
      ? message.target
      : null;
    const widgets = target ? [registry.widgets[target]] : Object.values(registry.widgets);
    widgets.forEach((widget) => {
      if (!widget || !document.body.contains(widget.elem)) return;
      widget[method](message);
    });
  };

  if (!registry.handlersInstalled) {
    Shiny.addCustomMessageHandler("librariesActiveSelectionChanged", (message) => {
      dispatchSelectionMessage(message, "onSelectionChanged");
    });
    Shiny.addCustomMessageHandler("librariesActiveSelectionReset", (message) => {
      dispatchSelectionMessage(message, "onSelectionReset");
    });
    registry.handlersInstalled = true;
  }
  // Deferred legend (layout.meta.rcDeferredLegend, R/page_organism_UMAP.R): the points are painted first and the
  // native legend is added in a later task, so it no longer delays the first paint.
  const meta = elem.layout && elem.layout.meta;
  if (meta && meta.rcDeferredLegend && elem.layout.showlegend === false && typeof window.setTimeout === "function") {
    const showLegend = () => {
      if (document.body.contains(elem) && elem.layout && elem.layout.showlegend === false) {
        Plotly.relayout(elem, { showlegend: true });
      }
    };
    const later = () => window.setTimeout(showLegend, 0);
    if (typeof window.requestAnimationFrame === "function") window.requestAnimationFrame(later); else later();
  }
  Shiny.setInputValue(valuesInputId + "_ready", Math.random(), { priority: "event" });
};
