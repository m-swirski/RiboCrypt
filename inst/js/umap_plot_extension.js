(elem, _, data) => {
  const onSelected = (e) => {
    const selectedPoints = e.points;
    const valuesInputId = data;
    const pointValuesToSet = selectedPoints.map((elem) => {
      return elem.text.split("<br />")[1].split(": ")[1]
    });
    Shiny.setInputValue(valuesInputId, pointValuesToSet);
    console.log("onSelected fired");
    console.log(pointValuesToSet);
  };
  elem.on("plotly_selected", onSelected);

  onChangedSelection = (message) => {
    let selected = null
    if(Array.isArray(message)) {
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
  Shiny.addCustomMessageHandler("samplesActiveSelectionChanged", onChangedSelection);
  Shiny.addCustomMessageHandler("samplesActiveFilteredSelectionChanged", onChangedSelection);
};
