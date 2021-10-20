(elem, _, data) => {
    let howManyTraces = elem.data.length;
    let isTraceVisible = Array(data.length).fill(0);
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
            data.forEach(traceDef, index => {
                if (distance <= traceDef["distance"] && isTraceVisible[index] == 0) {
                    let traceToAdd = {
                        "x": traceDef["x"],
                        "y": traceDef["y"],
                        "text": traceDef["text"],
                        "textfont": { "color": traceDef["color"] },
                        "xaxis": traceDef["xaxis"],
                        "yaxis": traceDef["yaxis"],
                        "mode": "text",
                    };
                    Plotly.addTraces(elem, traceToAdd, howManyTraces + index);
                    isTraceVisible[index] = 1;
                }
                if (distance > traceDef["distance"] && isTraceVisible[index] == 1) {
                    Plotly.deleteTraces(elem, howManyTraces + index);
                    isTraceVisible[index] = 0;
                }
            })
        }
    )
}