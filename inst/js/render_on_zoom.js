(elem, _, data) => {
  let tracesVisible = 0;
  let switchDistance = data[0]["distance"];
  let tracesToAdd = data.map(traceDef => {
      return {
          "x": traceDef["x"],
          "y": traceDef["y"],
          "text": traceDef["text"],
          "textfont": { "color": traceDef["color"] },
          "xaxis": traceDef["xaxis"],
          "yaxis": traceDef["yaxis"],
          "mode": "text",
      }
  });
  let indexesToDelete = Array(data.length).fill(0).map((_, index) => {
      return index - data.length
  });
  let onRelayout = (ed) => {
      if ("xaxis.range[0]" in ed && "xaxis.range[1]" in ed) {
          start = Math.floor(ed["xaxis.range[0]"]);
          end = Math.ceil(ed["xaxis.range[1]"]);
      } else {
          start = 0;
          end = Math.ceil(ed["width"]);
      }
      let distance = end - start;
      if (distance <= switchDistance && tracesVisible == 0) {
          Plotly.addTraces(elem, tracesToAdd);
          tracesVisible = 1;
      }
      if (distance > switchDistance && tracesVisible == 1) {
          Plotly.deleteTraces(elem, indexesToDelete);
          tracesVisible = 0;
      }
  };
  elem.on('plotly_relayout', onRelayout);
}
