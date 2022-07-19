(elem, _, data) => {
  let targetAxis = "y".concat(2 + data["nplots"]);
  let onClick = (ed) => {
    let point = ed["points"][0];
    let xaxis = point["data"]["xaxis"];
    let yaxis = point["data"]["yaxis"];
    if (xaxis === "x" && yaxis === targetAxis && "text" in point) {
      Shiny.setInputValue("selectedRegion", point["text"]);
    } else {
      Shiny.setInputValue("selectedRegion", null);
    }
  };
  elem.on('plotly_click', onClick);
}
