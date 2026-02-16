(elem, _, data) => {
  const onSelected = (e) => {
    const selectedPoints = e.points;
    const valuesInputId = data;
    const pointValuesToSet = selectedPoints.map((elem) => {
      return elem.text.split("<br />")[1].split(": ")[1]
    });
    Shiny.setInputValue(valuesInputId, pointValuesToSet);
  };
  elem.on("plotly_selected", onSelected);

  onChangedSelection = (message) => {
    let runId_counter = 0;
    const tracesToUpdate = Array.from({ length: elem.data.length }, (_, i) => i + 1);
    let indicesToUpdate = Array.from({ length: elem.data.length }, (_, i) => []);

    while (runId_counter < message.length) {
      elem.data.forEach((trace, index) => {
        let indices_buffer = [];
        if (Array.isArray(trace.text)) {
          trace.text.forEach((t, tIndex) => {
            if (t.includes(message[runId_counter])) {
              indicesToUpdate[index].push(tIndex);
              runId_counter += 1;
            };
          });
        } else {
          if (trace.text.includes(message[runId_counter])) {
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
