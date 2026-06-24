(elem, _, data) => {
  Plotly.relayout(elem, { dragmode: "select" });
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

  const onSelected = (e) => {
    if (suppressSelectionEvents > 0) return;
    const selectedPoints = e && Array.isArray(e.points) ? e.points : [];
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

  elem.on("plotly_selected", onSelected);
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
};
