(elem, _, data) => {
  Plotly.relayout(elem, { dragmode: "select" });
  const valuesInputId = data;
  const sendSelection = (selection) => {
    Shiny.setInputValue(valuesInputId, selection, { priority: "event" });
  };

  const onSelected = (e) => {
    const selectedPoints = e && Array.isArray(e.points) ? e.points : [];
    const pointValuesToSet = selectedPoints.map((elem) => {
      return elem.text.split("<br />")[1].split(": ")[1]
    });
    sendSelection(pointValuesToSet);
    console.log("onSelected fired");
    console.log(pointValuesToSet);
  };
  const onDeselected = () => {
    sendSelection([]);
    console.log("onDeselected fired");
  };

  elem.on("plotly_selected", onSelected);
  elem.on("plotly_deselect", onDeselected);
  elem.on("plotly_doubleclick", onDeselected);

  onSelectionChanged = (message) => {
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

    Plotly.react(elem, updatedData, updatedLayout);
  };
  Shiny.addCustomMessageHandler("librariesActiveSelectionChanged", onSelectionChanged);

  onSelectionReset = (_) => {
    const tracesToUpdate = Array.from({ length: elem.data.length }, (_, i) => i + 1);

    let updatedData = [...elem.data]
    tracesToUpdate.forEach((t, index) => {
      updatedData[index]["selectedpoints"] = []
    });

    let updatedLayout = { ...elem.layout };
    updatedLayout["selections"] = [];

    Plotly.react(elem, updatedData, updatedLayout);
  };
  Shiny.addCustomMessageHandler("librariesActiveSelectionReset", onSelectionReset);
};
