(elem, _, data) => {
  Plotly.relayout(elem, { dragmode: "select" });
  const valuesInputId = data;
  const sendSelection = (selection) => {
    Shiny.setInputValue(valuesInputId, selection, { priority: "event" });
  };
  let suppressSelectionEvents = 0;

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

  const onSelected = (e) => {
    if (suppressSelectionEvents > 0) return;
    const selectedPoints = e && Array.isArray(e.points) ? e.points : [];
    const pointValuesToSet = selectedPoints.map((elem) => {
      return elem.text.split("<br />")[1].split(": ")[1]
    });
    sendSelection(pointValuesToSet);
    console.log("onSelected fired");
    console.log(pointValuesToSet);
  };
  const onDeselected = () => {
    if (suppressSelectionEvents > 0) return;
    sendSelection([]);
    console.log("onDeselected fired");
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
    let selected = null
    if (Array.isArray(message)) {
      selected = message;
    } else {
      selected = [message];
    };

    let runId_counter = 0;
    const tracesToUpdate = Array.from({ length: elem.data.length }, (_, i) => i + 1);
    let indicesToUpdate = Array.from({ length: elem.data.length }, (_, i) => []);

    while (runId_counter < selected.length) {
      elem.data.forEach((trace, index) => {
        if (Array.isArray(trace.text)) {
          trace.text.forEach((t, tIndex) => {
            if (t.includes(selected[runId_counter])) {
              indicesToUpdate[index].push(tIndex);
              runId_counter += 1;
            };
          });
        } else {
          if (trace.text.includes(selected[runId_counter])) {
            indicesToUpdate[index].push(1)
            runId_counter += 1;
          }
        }
      });
      runId_counter += 1;
    };

    let updatedData = [...elem.data]
    tracesToUpdate.forEach((t, index) => {
      updatedData[index]["selectedpoints"] = indicesToUpdate[index]
    });

    let updatedLayout = { ...elem.layout };
    updatedLayout["selections"] = [];

    runProgrammaticSelectionSync(() => Plotly.react(elem, updatedData, updatedLayout));
  };
  Shiny.addCustomMessageHandler("librariesActiveSelectionChanged", onSelectionChanged);

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
  Shiny.addCustomMessageHandler("librariesActiveSelectionReset", onSelectionReset);
};
