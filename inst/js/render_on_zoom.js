(elem, _, data) => {
  let tracesVisible = 0;
  let switchDistance = 300;
  let targetAxis = data.traces[0]["yaxis"];
  console.log("Target Y axis:", targetAxis);
  let currentVisibleSequence = "";
  let lastRange = [null, null];

  // Trace for DNA letters when zoom is < switchDistance
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
        type: "scatter",
        name: "sequence"
      };
    });
  };
  // Logic to fire DNA letters or placeholder text
  const onRelayout = (ed) => {
    const fallbackRange = elem._fullLayout?.xaxis?.range || [0, data.sequence.length];
    const start = Math.floor("xaxis.range[0]" in ed ? ed["xaxis.range[0]"] : fallbackRange[0]);
    const end = Math.ceil("xaxis.range[1]" in ed ? ed["xaxis.range[1]"] : fallbackRange[1]);
    const distance = end - start;

    // ✅ Prevent infinite relayout loop
    if (lastRange[0] === start && lastRange[1] === end) {
      return;
    }
    lastRange = [start, end];

    const clampedStart = Math.max(0, start);
    const clampedEnd = Math.min(data.sequence.length, end);
    currentVisibleSequence = data.sequence.slice(clampedStart, clampedEnd);
    console.log("Updated visible sequence of length:", currentVisibleSequence.length);

    if (distance <= switchDistance && tracesVisible === 0) {
      const placeholderIndex = elem.data.findIndex((trace) => trace.name === "sequence_placeholder");
      if (placeholderIndex !== -1) {
        Plotly.deleteTraces(elem, [placeholderIndex]);
      }

      Plotly.addTraces(elem, tracesToAdd(data.traces, start - 300, end + 300));
      tracesVisible = 1;
    } else if (distance > switchDistance && tracesVisible === 1) {
      const indexesToRemove = elem.data
        .map((trace, idx) => (trace.type === "scatter" && trace.mode === "text" && trace.name === "sequence") ? idx : null)
        .filter((i) => i !== null);

      if (indexesToRemove.length > 0) {
        Plotly.deleteTraces(elem, indexesToRemove);
      }

      tracesVisible = 0;

      const xMid = (start + end) / 2;
      Plotly.addTraces(elem, [{
        x: [xMid],
        y: [0.5],
        text: "\u00A0\u00A0\u00A0Click to copy sequence\u00A0\u00A0\u00A0",
        mode: "markers+text",
        type: "scatter",
        textposition: "middle center",
        textfont: { color: "gray", size: 18 },
        marker: {
          size: 100,
          color: "rgba(0,0,0,0)",
          opacity: 0.01, // Very low opacity to keep it interactable
          line: { width: 0 }
        },
        name: "sequence_placeholder",
        hoverinfo: "text", // Make the point interactive on hover
        showlegend: false,
        xaxis: "x",
        yaxis: targetAxis
      }]);

      Plotly.relayout(elem, {
        "xaxis.showticklabels": false,
        "xaxis.ticks": "",
        "xaxis.showline": false,
        "xaxis.showgrid": false,
        "xaxis.zeroline": false
      });
    }
  };

  // Function to display a message in the bottom left of the plot
  const showCopiedMessage = (message) => {
    const messageDiv = document.createElement("div");
    messageDiv.innerHTML = message;
    messageDiv.style.position = "absolute";
    messageDiv.style.bottom = "10px";
    messageDiv.style.left = "10px";
    messageDiv.style.backgroundColor = "rgba(0, 0, 0, 0.7)";
    messageDiv.style.color = "white";
    messageDiv.style.padding = "5px 10px";
    messageDiv.style.borderRadius = "5px";
    messageDiv.style.fontSize = "14px";
    messageDiv.style.zIndex = "10";  // Ensure it's on top of the plot
    document.body.appendChild(messageDiv);

    // Remove the message after 3 seconds
    setTimeout(() => {
      document.body.removeChild(messageDiv);
    }, 3000);
  };

  // Modified onClick to trigger the copied message
  const onClick = (ed) => {
    const point = ed.points?.[0];
    if (!point) return;

    const xaxis = point.data.xaxis;
    const yaxis = point.data.yaxis;
    const traceName = point.data.name;

    console.log("Clicked axis:", xaxis, yaxis);
    console.log("Trace name:", traceName);

    if (xaxis === "x" && yaxis === targetAxis && traceName && (traceName === "sequence" || traceName === "sequence_placeholder")) {
      // Get visible window from layout
      const layout = elem._fullLayout || {};
      const range = layout.xaxis?.range || [0, data.sequence.length];
      const start = Math.max(0, Math.floor(range[0]));
      const end = Math.min(data.sequence.length, Math.ceil(range[1]));
      const visibleSequence = data.sequence.slice(start, end);

      console.log("Copied sequence of length:", visibleSequence.length);
      Shiny.setInputValue(data.input_id, visibleSequence, { priority: "event" });

      // Trigger the copy action to clipboard
      navigator.clipboard.writeText(visibleSequence)
        .then(() => {
          console.log("Copied to clipboard!");
          showCopiedMessage(`Sequence copied: ${visibleSequence.length} nt`); // Show the copied message
        })
        .catch(err => console.error("Clipboard copy failed:", err));
    } else {
      console.log("Not copying, not DNA sequence or placeholder");
      Shiny.setInputValue(data.input_id, null);
    }
  };

  // Add text on initial load: Ensure initial rendering of "Click to copy" marker
  Plotly.addTraces(elem, [{
    x: [data.sequence.length / 2],
    y: [0.5],
    text: "\u00A0\u00A0\u00A0Click to copy sequence\u00A0\u00A0\u00A0",
    mode: "markers+text",
    type: "scatter",
    textposition: "middle center",
    textfont: { color: "gray", size: 18 },
    marker: {
      size: 100,
      color: "rgba(0,0,0,0)",
      opacity: 0.01, // Keep it invisible but clickable
      line: { width: 0 }
    },
    name: "sequence_placeholder",
    hoverinfo: "text",
    showlegend: false,
    xaxis: "x",
    yaxis: targetAxis
  }]);

  // Initialize the plot with an initial range (force render text)
  const initialStart = 0;
  const initialEnd = data.sequence.length;
  onRelayout({ "xaxis.range[0]": initialStart, "xaxis.range[1]": initialEnd });

  elem.on('plotly_click', onClick);
  elem.on("plotly_relayout", onRelayout);
};
