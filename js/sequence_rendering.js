(elem, _, data) => {
    let isSequenceVisible = 0;
    elem.on('plotly_relayout',
        (ed) => {
            if ("xaxis.range[0]" in ed && "xaxis.range[1]" in ed) {
                start = Math.floor(ed["xaxis.range[0]"]);
                end = Math.ceil(ed["xaxis.range[1]"]);
            } else {
                start = 0;
                end = Math.ceil(ed["width"]);
            }
            let distance = end - start;
            if (distance <= 50 && isSequenceVisible == 0) {
                let sequenceIndex = ["frame0", "frame1", "frame2"];
                let sequenceTraces = sequenceIndex.map((val, idx) => {
                    return {
                        "x": Array.from(data[val], (_, i) => i + idx),
                        "y": Array.from(data[val], (_) => 0.5),
                        "text": data[val],
                        "textfont": { "color": data["colors"][idx] },
                        "xaxis": "x",
                        "yaxis": "y2",
                        "mode": "text",
                    }
                });
                Plotly.addTraces(elem, sequenceTraces);
                isSequenceVisible = 1;
            }
            if (distance > 50 && isSequenceVisible == 1) {
                Plotly.deleteTraces(elem, [-1]);
                isSequenceVisible = 0;
            }
        }
    )
}