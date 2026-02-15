(elem, _, data) => {
  const onSelected = (e) => {
    const selectedPoints = e.points;
    const valuesInputId = data;
    const pointValuesToSet = {
      "sample": selectedPoints.map((elem) => {
        return elem.text.split("<br />")[1].split(": ")[1]
      }),
      "curveIndex": selectedPoints.map((elem) => {
        return elem.curveNumber
      }),
      "pointIndex": selectedPoints.map((elem) => {
        return elem.pointIndex
      })
    };
    Shiny.setInputValue(valuesInputId, pointValuesToSet);
  };
  elem.on("plotly_selected", onSelected);

  onChangedSelection = (message) => {
    const currentData = elem.data;

    currentData.forEach((trace, index) => {
      const foundIndex = message.curveIndex.findIndex((i) => i === index)
      if (foundIndex > -1) {
        trace.selectedpoints = message.pointIndex[foundIndex];
      } else {
        trace.selectedpoints = [];
      }
    });

    Plotly.react(elem, currentData, elem.layout)
  };
  Shiny.addCustomMessageHandler("samplesActiveSelectionChanged", onChangedSelection);
  Shiny.addCustomMessageHandler("samplesActiveFilteredSelectionChanged", onChangedSelection);
};
