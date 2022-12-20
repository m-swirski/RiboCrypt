(elem, _, data) => {
  console.log(data);
  let targetAxis = "y".concat(2 + data["nplots"]);
  let onClick = (ed) => {
    let point = ed["points"][0];
    let xaxis = point["data"]["xaxis"];
    let yaxis = point["data"]["yaxis"];
    if (xaxis === "x" && yaxis === targetAxis && "text" in point["data"]) {
      Shiny.setInputValue(data["input_id"], point["data"]["text"], {priority: "event"});
    } else {
      Shiny.setInputValue(data["input_id"], null);
    }
  };
  elem.on('plotly_click', onClick);
}
