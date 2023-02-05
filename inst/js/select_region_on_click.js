(elem, _, data) => {
  console.log(data);
  let targetAxis = "y".concat(2 + data["nplots"]);
  let onClick = (ed) => {
    let point = ed["points"][0];
    let xaxis = point["data"]["xaxis"];
    let yaxis = point["data"]["yaxis"];
    let region = null;
    if ("text" in point["data"] && point["data"]["text"].constructor === String) {
      region = point["data"]["text"];
    } else if ("text" in point) {
      region = point["text"];
    } else {
      region = null;
    }
    if (xaxis === "x" && yaxis === targetAxis && region !== null) {
      Shiny.setInputValue(data["input_id"], region, {priority: "event"});
    } else {
      Shiny.setInputValue(data["input_id"], null);
    }
  };
  elem.on('plotly_click', onClick);
}
