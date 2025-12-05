(elem, _, data) => {
  let tracesVisible = 0;
  const switchDistance = data.traces[0]["distance"];
  const targetAxis = data.traces[0]["yaxis"];
  let currentVisibleSequence = "";
  let lastRange = [null, null];

  // ⚙️ Add DNA base traces
  const tracesToAdd = (traces, start, end) => {
    return traces.map((traceDef, index) => {
      const frame = index % 3;

      const startFrame = (start - 1) % 3;
      const startShift = (startFrame > frame) ? 1 : 0;
      const frameAdjustedStart = Math.max(0, Math.floor(start / 3) + startShift);

      const endFrame = (end - 1) % 3;
      const endShift = (endFrame < frame) ? -1 : 0;
      const frameAdjustedEnd = Math.floor(end / 3) + endShift;

      const xSeg = traceDef.x.slice(frameAdjustedStart, frameAdjustedEnd);
      const ySeg = traceDef.y.slice(frameAdjustedStart, frameAdjustedEnd);
      if (!xSeg.length || !ySeg.length) return null;   // skip empty traces
      return {
        x: xSeg, y: ySeg,
        text: traceDef.text.slice(frameAdjustedStart, frameAdjustedEnd),
        textfont: { color: traceDef.color },
        xaxis: traceDef.xaxis,
        yaxis: traceDef.yaxis,
        mode: "text",
        type: "scattergl",
        name: "sequence",
        hovertemplate: "Pos %{x:.0f}"
      };
    });
  };

  // 👁️ Placeholder trace
  const createPlaceholderTrace = (centerX) => ({
    x: [centerX],
    y: [0.5],
    text: "\u00A0\u00A0\u00A0Click to copy sequence\u00A0\u00A0\u00A0",
    mode: "markers+text",
    type: "scatter",
    textposition: "middle center",
    textfont: { color: "gray", size: 18 },
    marker: {
      size: 100,
      color: "rgba(0,0,0,0)",
      opacity: 0.01,
      line: { width: 0 }
    },
    name: "sequence_placeholder",
    hoverinfo: "text",
    showlegend: false,
    xaxis: "x",
    yaxis: targetAxis
  });

  // 🔁 Triggered on zoom/pan
  const onRelayout = (ed) => {
    const fallbackRange = elem._fullLayout?.xaxis?.range || [0, data.sequence.length];
    const start = Math.floor("xaxis.range[0]" in ed ? ed["xaxis.range[0]"] : fallbackRange[0]);
    const end = Math.ceil("xaxis.range[1]" in ed ? ed["xaxis.range[1]"] : fallbackRange[1]);
    const distance = end - start;

    if (lastRange[0] === start && lastRange[1] === end) return;
    lastRange = [start, end];

    const clampedStart = Math.max(0, start);
    const clampedEnd = Math.min(data.sequence.length, end);
    currentVisibleSequence = data.sequence.slice(clampedStart, clampedEnd);
    console.log("Updated visible sequence of length:", currentVisibleSequence.length);

    const placeholderIndex = elem.data.findIndex(trace => trace.name === "sequence_placeholder");
    const sequenceTraceIndexes = elem.data
      .map((trace, i) => trace.name === "sequence" ? i : null)
      .filter(i => i !== null);

    if (distance <= switchDistance) {
      // Show DNA text traces
      if (placeholderIndex !== -1) Plotly.deleteTraces(elem, [placeholderIndex]);
      if (tracesVisible === 0) {
        Plotly.addTraces(elem, tracesToAdd(data.traces, start - 300, end + 300));
        tracesVisible = 1;
      }
      Plotly.relayout(elem, {
        "xaxis.showticklabels": true,
        "xaxis.ticks": "outside",
        "xaxis.showline": false,
        "xaxis.showgrid": false,
        "xaxis.zeroline": false
      });
    } else {
      // Show placeholder trace
      const centerX = (start + end) / 2;

      // Remove DNA traces if needed
      const sequenceTraceIndexes = elem.data
        .map((trace, i) => trace.name === "sequence" ? i : null)
        .filter(i => i !== null);
      if (sequenceTraceIndexes.length > 0) {
        Plotly.deleteTraces(elem, sequenceTraceIndexes);
        tracesVisible = 0;
      }

      // Add or update placeholder
      if (placeholderIndex !== -1) {
        Plotly.restyle(elem, { x: [[centerX]] }, [placeholderIndex]);
      } else {
        Plotly.addTraces(elem, [createPlaceholderTrace(centerX)]);
      }

      // Clean up axes visuals
      Plotly.relayout(elem, {
        "xaxis.showticklabels": true,
        "xaxis.ticks": "outside",
        "xaxis.showline": false,
        "xaxis.showgrid": false,
        "xaxis.zeroline": false
      });
    }
  };

  // Show copied popup
  const showCopiedMessage = (message, sequence) => {
    const messageDiv = document.createElement("div");
    const encodedSeq = encodeURIComponent(sequence);

    // Build BLAST URL with pre-filled sequence
    const blastUrl =
      "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_SPEC=GeoBlast&PAGE_TYPE=BlastSearch&QUERY="
      + encodedSeq;


    messageDiv.innerHTML = `
      ${message}
      &nbsp; <a href="${blastUrl}" target="_blank" style="color:#4DA3FF; text-decoration:underline;">
        Blast
      </a>
    `;
    Object.assign(messageDiv.style, {
      position: "absolute",
      bottom: "10px",
      left: "10px",
      backgroundColor: "rgba(0, 0, 0, 0.7)",
      color: "white",
      padding: "5px 10px",
      borderRadius: "5px",
      fontSize: "14px",
      zIndex: "10"
    });
    document.body.appendChild(messageDiv);
    setTimeout(() => document.body.removeChild(messageDiv), 3000);
  };

  //  Click handler
  const onClick = (ed) => {
    const point = ed.points?.[0];
    if (!point) return;

    const { xaxis, yaxis, name } = point.data;
    if (xaxis === "x" && yaxis === targetAxis &&
        (name === "sequence" || name === "sequence_placeholder")) {
      const layout = elem._fullLayout || {};
      const range = layout.xaxis?.range || [0, data.sequence.length];
      const start = Math.max(0, Math.floor(range[0]));
      const end = Math.min(data.sequence.length, Math.ceil(range[1]));
      const visibleSequence = data.sequence.slice(start, end);

      console.log("Copied sequence of length:", visibleSequence.length);
      Shiny.setInputValue(data.input_id, visibleSequence, { priority: "event" });

      navigator.clipboard.writeText(visibleSequence)
        .then(() => {
          console.log("Copied to clipboard!");
          showCopiedMessage(
            `Sequence copied: ${visibleSequence.length} nt`,
            visibleSequence
          );
        })
        .catch(err => console.error("Clipboard copy failed:", err));
    }
  };

  // Initial render with placeholder
  Plotly.addTraces(elem, [createPlaceholderTrace(data.sequence.length / 2)]);
  const initialStart = 0;
  const initialEnd = data.sequence.length;
  onRelayout({ "xaxis.range[0]": initialStart, "xaxis.range[1]": initialEnd });

  // Bind events
  elem.on('plotly_click', onClick);
  elem.on("plotly_relayout", onRelayout);
}
