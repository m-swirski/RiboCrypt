(elem, _, data) => {
  const onSelected = (e) => {
    console.log(e);
    const inputId = data;
    const selectedPoints = e.points;
    const valueToSet = selectedPoints.map((elem) => {
      return elem.text.split("<br />")[1].split(": ")[1];    
    });
    console.log(valueToSet);
    Shiny.setInputValue(inputId, valueToSet);
  };
  elem.on("plotly_selected", onSelected);
};
