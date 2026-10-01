(elem, _, config) => {
  if (elem._localYCleanup) elem._localYCleanup();
  if (config.manual || !config.axes.length) return;
  const checkbox = document.getElementById(config.input_id);
  const originals = new Map(config.axes.map(name => {
    const axis = elem._fullLayout[name];
    return [name, {range: axis.range.slice(), tickmode: axis.tickmode,
      tickvals: axis.tickvals || null, ticktext: axis.ticktext || null}];
  }));
  let scheduled = null;
  let disposed = false;
  const finite = value => typeof value === 'number' && Number.isFinite(value);

  // Include boundary intersections of lines and filled polygons, not just vertices.
  const traceMaximum = (trace, low, high) => {
    let maximum = 0;
    const x = trace.x || [], y = trace.y || [];
    const lines = trace.type !== 'bar' && (trace.fill && trace.fill !== 'none' ||
      (trace.mode || '').includes('lines'));
    for (let i = 0; i < x.length; i++) {
      if (!finite(x[i]) || !finite(y[i])) continue;
      if (x[i] >= low && x[i] <= high) maximum = Math.max(maximum, y[i]);
      if (!lines || !i || !finite(x[i - 1]) || !finite(y[i - 1]) || x[i] === x[i - 1]) continue;
      for (const edge of [low, high]) {
        if (edge < Math.min(x[i - 1], x[i]) || edge > Math.max(x[i - 1], x[i])) continue;
        maximum = Math.max(maximum, y[i - 1] +
          (y[i] - y[i - 1]) * (edge - x[i - 1]) / (x[i] - x[i - 1]));
      }
    }
    return maximum;
  };

  const update = () => {
    scheduled = null;
    if (disposed || !elem._fullLayout) return;
    const range = elem._fullLayout.xaxis.range;
    const low = Math.min(...range), high = Math.max(...range);
    const full = config.full_range;
    const local = (!checkbox || checkbox.checked) && !(low <= full[0] && high >= full[1]);
    const changes = {};
    for (const [name, original] of originals) {
      const axis = elem._fullLayout[name];
      const ref = name.replace('axis', '');
      let maximum = 0;
      if (local) for (const [index, trace] of elem.data.entries()) {
        if ((trace.yaxis || 'y') !== ref || trace.visible === false || trace.visible === 'legendonly') continue;
        const hover = Array.isArray(trace.hoverinfo) ? trace.hoverinfo : [trace.hoverinfo];
        if (hover.length && hover.every(value => value === 'skip')) continue;
        // Plotly's calculated stack coordinates include interpolation and cumulative height.
        const calculated = trace.stackgroup && elem.calcdata?.[index];
        const values = calculated ? {...trace, x: calculated.map(point => point.x),
          y: calculated.map(point => point.y)} : trace;
        maximum = Math.max(maximum, traceMaximum(values, low, high));
      }
      const desired = local ? {range: [0, maximum > 0 ? maximum * 1.05 : 1],
        tickmode: 'auto', tickvals: null, ticktext: null} : original;
      for (const [key, value] of Object.entries(desired)) {
        if (JSON.stringify(axis[key] ?? null) !== JSON.stringify(value)) changes[name + '.' + key] = value;
      }
      if (axis.autorange !== false) changes[name + '.autorange'] = false;
    }
    if (Object.keys(changes).length) Plotly.relayout(elem, changes);
  };
  const schedule = () => {
    if (scheduled === null && !disposed) scheduled = requestAnimationFrame(update);
  };
  const onRelayout = event => {
    if (Object.keys(event).some(key => /^xaxis\d*\.(range|autorange)/.test(key))) schedule();
  };
  const namespace = '.localY_' + elem.id.replace(/[^a-zA-Z0-9]/g, '_');
  if (checkbox) $(checkbox).on('change' + namespace, schedule);
  elem.on('plotly_relayout', onRelayout);
  elem.on('plotly_restyle', schedule);
  elem.on('plotly_animated', schedule);
  elem._localYCleanup = () => {
    disposed = true;
    if (scheduled !== null) cancelAnimationFrame(scheduled);
    if (checkbox) $(checkbox).off('change' + namespace, schedule);
    elem.removeListener('plotly_relayout', onRelayout);
    elem.removeListener('plotly_restyle', schedule);
    elem.removeListener('plotly_animated', schedule);
  };
  schedule();
}
