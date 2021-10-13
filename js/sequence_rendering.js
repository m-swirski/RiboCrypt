(elem, data) => {
    var active = 0;
    elem.on('plotly_relayout',
        (ed) => {
            if ("xaxis.range[0]" in ed && "xaxis.range[1]" in ed) {
                distance = ed["xaxis.range[1]"] - ed["xaxis.range[0]"]
            } else {
                distance = ed["width"]
            }
            if (distance <= 50 && active == 0) {
                let seqLength = elem.data[1].x.length;
                let seqTrace = {
                    "text": elem.data[1].text.map(
                        (x) => { return x.split("<br />")[2].split(": ")[1] }
                    ),
                    "x": elem.data[1].x,
                    "xaxis": "x",
                    "y": Array(seqLength).fill(0),
                    "yaxis": "y2",
                    "mode": "text"
                };
                Plotly.addTraces(elem, seqTrace);
                active = 1;
            }
            if (distance > 50 && active == 1) {
                Plotly.deleteTraces(elem, -1);
            }
            console.log(elem.data);
        }
    )
}