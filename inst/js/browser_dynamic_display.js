(elem, _, config) => {
  if (elem._rcDynamicCleanup) elem._rcDynamicCleanup();
  if (!config.active) return;
  const category = span => span > 10000 ? 'lines' : span > 500 ? 'area' : 'columns';
  let current = config.type;
  let timer;
  const onRelayout = event => {
    if (!Object.keys(event).some(k => /^xaxis\d*\.(range|autorange)/.test(k))) return;
    clearTimeout(timer);
    timer = setTimeout(() => {
      const range = elem._fullLayout?.xaxis?.range;
      if (!range || range.length !== 2) return;
      const a = Number(range[0]), b = Number(range[1]);
      if (!Number.isFinite(a) || !Number.isFinite(b) || a === b) return;
      const lo = Math.min(a, b), hi = Math.max(a, b);
      const next = category(hi - lo + 1);
      if (next === current) return;
      current = next;
      Shiny.setInputValue(config.input_id, {range: [lo, hi], context: config.context}, {priority: 'event'});
    }, 150);
  };
  elem.on('plotly_relayout', onRelayout);
  elem._rcDynamicCleanup = () => { clearTimeout(timer); elem.removeListener('plotly_relayout', onRelayout); };
};
