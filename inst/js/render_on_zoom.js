(elem, _, data) => {
  let tracesVisible = 0;
  let switchDistance = data.traces[0].distance;
  let clickHandlerAttached = false;

  const tracesToAdd = (traces, start, end) => {
    return traces.map((traceDef, index) => {
      let frame = index % 3;

      let startFrame = (start - 1) % 3;
      let startShift = (startFrame > frame) ? 1 : 0;
      let frameAdjustedStart = Math.max(0, Math.floor(start / 3) + startShift);

      let endFrame = (end - 1) % 3;
      let endShift = (endFrame < frame) ? -1 : 0;
      let frameAdjustedEnd = Math.floor(end / 3) + endShift;

      return {
        x: traceDef.x.slice(frameAdjustedStart, frameAdjustedEnd),
        y: traceDef.y.slice(frameAdjustedStart, frameAdjustedEnd),
        text: traceDef.text.slice(frameAdjustedStart, frameAdjustedEnd),
        textfont: { color: traceDef.color },
        xaxis: traceDef.xaxis,
        yaxis: traceDef.yaxis,
        mode: "text",
        type: "scatter"
      };
    });
  };

  const indexesToDelete = Array(data.traces.length).fill(0).map((_, index) => index - data.traces.length);

  const copySequenceFrom = (eventData) => {
    if (!data || !data.sequence || data.sequence.length === 0) {
      console.warn("No sequence data to copy");
      return;
    }

    let yaxis = null;
    if (eventData && eventData.points && eventData.points.length > 0) {
      const trace = eventData.points[0].data;
      const yaxis = trace.yaxis || "";
      const expected_yaxis = data.traces[0].yaxis || "y";
      if (yaxis !== expected_yaxis) {
        console.log("Click ignored: not on DNA track");
        return;
      }
    }


    const xRange = elem._fullLayout.xaxis.range;
    const start = Math.max(0, Math.floor(xRange[0]) - 1);
    const end = Math.min(data.sequence.length, Math.ceil(xRange[1]));
    const visibleSeq = data.sequence.slice(start, end);

    const textarea = document.createElement("textarea");
    textarea.value = visibleSeq;
    textarea.setAttribute("readonly", "");
    textarea.style.position = "absolute";
    textarea.style.left = "-9999px";
    document.body.appendChild(textarea);
    textarea.select();

    try {
      document.execCommand("copy");
      console.log(`Copied ${visibleSeq.length} nt [${start}:${end}]`);
      if (data.input_id) {
        Shiny.setInputValue(data.input_id, visibleSeq, { priority: "event" });
      }
    } catch (err) {
      console.error("Clipboard copy failed:", err);
    }

    document.body.removeChild(textarea);
  };

  const attachClickHandler = () => {
    if (clickHandlerAttached) return;
    clickHandlerAttached = true;

    elem.on("plotly_click", (eventData) => {
      console.log("plotly_click triggered");
      copySequenceFrom(eventData);
    });
  };

  const onRelayout = (ed) => {
    let start, end;
    if ("xaxis.range[0]" in ed && "xaxis.range[1]" in ed) {
      start = Math.floor(ed["xaxis.range[0]"]);
      end = Math.ceil(ed["xaxis.range[1]"]);
    } else if ("xaxis.autorange" in ed && ed["xaxis.autorange"]) {
      const fullRange = elem.layout.xaxis.range;
      start = Math.floor(fullRange[0]);
      end = Math.ceil(fullRange[1]);
    } else {
      const fullRange = elem.layout.xaxis.range;
      start = Math.floor(fullRange[0]);
      end = Math.ceil(fullRange[1]);
    }

    const distance = end - start;
    if (distance <= switchDistance && tracesVisible === 0) {
      Plotly.addTraces(elem, tracesToAdd(data.traces, start - 300, end + 300));
      tracesVisible = 1;
    }

    if (distance > switchDistance && tracesVisible === 1) {
      Plotly.deleteTraces(elem, indexesToDelete);
      tracesVisible = 0;
    }
  };

  attachClickHandler();
  elem.on("plotly_relayout", onRelayout);
  setTimeout(() => onRelayout({}), 100);
};
