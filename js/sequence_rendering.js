(elem, data) => {
    let isSequenceVisible = 0;
    let replacedTrace = null;
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
                let sequencePositions = elem.data[1].x;
                let sequenceText = elem.data[1].text;
                let sequenceTrace = {
                    "x": sequencePositions,
                    "y": Array(sequencePositions.length).fill(0.5),
                    "text": sequenceText.map(
                        (x) => { return x.split("<br />")[2].split(": ")[1] }
                    ),
                    "xaxis": "x",
                    "yaxis": "y2",
                    "mode": "text"
                };
                replacedTrace = elem.data[5];
                Plotly.deleteTraces(elem, 5);
                Plotly.addTraces(elem, sequenceTrace, 5);
                isSequenceVisible = 1;
            }
            if (distance > 50 && isSequenceVisible == 1) {
                Plotly.deleteTraces(elem, 5);
                Plotly.addTraces(elem, replacedTrace, 5);
                replacedTrace = null;
                isSequenceVisible = 0;
            }
        }
    )
}