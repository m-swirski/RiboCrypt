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
    console.log(message);
  };
  Shiny.addCustomMessageHandler("samplesSelectionChanged", onChangedSelection);
};
