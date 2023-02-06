(elem, _, data) => {
  let tracesVisible = 0;
  let switchDistance = data[0]["distance"];
  let tracesToAdd = (data, start, end) => { 
    return data.map((traceDef, index) => {
      let frame = index % 3;
      
      let startFrame = (start - 1) % 3;
      let startShift = 0;
      if (startFrame > frame) {
        startShift = 1
      };
      let frameAdjustedStart = Math.floor(start / 3) + startShift;
      
      let endFrame = (end - 1) % 3;
      let endShift = 0;
      if (endFrame < frame) {
        startShift = -1
      };
      let frameAdjustedEnd = Math.floor(end / 3) + endShift;
      
      return {
          "x": traceDef["x"].slice(frameAdjustedStart, frameAdjustedEnd),
          "y": traceDef["y"].slice(frameAdjustedStart, frameAdjustedEnd),
          "text": traceDef["text"].slice(frameAdjustedStart, frameAdjustedEnd),
          "textfont": { "color": traceDef["color"] },
          "xaxis": traceDef["xaxis"],
          "yaxis": traceDef["yaxis"],
          "mode": "text",
      }
    });
  };
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
          Plotly.addTraces(elem, tracesToAdd(data, start - 100, end + 100));
          tracesVisible = 1;
      }
      if (distance > switchDistance && tracesVisible == 1) {
          Plotly.deleteTraces(elem, indexesToDelete);
          tracesVisible = 0;
      }
  };
  elem.on('plotly_relayout', onRelayout);
}
