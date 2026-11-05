(elem, _, config) => {
  if (elem._rcDynamicCleanup) elem._rcDynamicCleanup();
  elem.__rcDynamicActive = true;
  const groups = {};
  for (let i = 0; i < (elem.data || []).length; i++) {
    const meta = elem.data[i].meta || {};
    if (!meta.rc_columns_switch) continue;
    const key = meta.rc_columns_group + ':' + meta.rc_columns_frame;
    const g = groups[key] ||= {track: Number(String(meta.rc_columns_group).replace('track_', '')) - 1,
      frame: String(meta.rc_columns_frame)};
    g[meta.rc_columns_switch === 'gl_subset' ? 'gl' : 'line'] = i;
  }
  const entries = Object.values(groups).filter(g => g.gl !== undefined && g.line !== undefined)
    .map(g => ({...g, x: elem.data[g.line].x, raw: elem.data[g.line].y,
      smooth: config.tracks[g.track]?.[g.frame]}));
  if (!entries.length || entries.some(g => !g.smooth)) {
    console.error('Dynamic browser profiles could not be matched to tracks');
    return;
  }
  let running = false, pending = false, disposed = false, lastKey = '', lastMode = '';
  const hidden = new Set();
  const lowerBound = (x, value) => { let lo=0, hi=x.length; while(lo<hi) { const mid=(lo+hi)>>1; if(x[mid]<value) lo=mid+1; else hi=mid; } return lo; };
  const update = async () => {
    if (disposed) return;
    if (running) { pending = true; return; }
    const range = elem._fullLayout?.xaxis?.range;
    if (!range || range.length !== 2) return;
    const lo = Math.min(...range), hi = Math.max(...range), span = hi-lo+1;
    if (!Number.isFinite(span) || span <= 0) return;
    const mode = span > 10000 ? 'lines' : span > 500 ? 'area' : 'columns';
    const key = `${mode}:${lo}:${hi}`;
    if (key === lastKey) return;
    lastKey = key; running = true;
    try {
      const indices = entries.map(g => g.line).concat(entries.map(g => g.gl));
      const width = Math.max(1, Math.min(12, Math.round((elem._fullLayout.xaxis._length || 800) / span * 8) / 10));
      const lineData = entries.map(g => {
        if (mode === 'columns') return {x:[null], y:[null]};
        // Keep neighbours for interpolation at viewport boundaries.
        const start = Math.max(0, lowerBound(g.x, lo) - 1);
        const end = Math.min(g.x.length, lowerBound(g.x, hi) + 1);
        return {x:g.x.slice(start, end), y:g.smooth.y.slice(start, end)};
      });
      const glData = entries.map(g => {
        if (mode !== 'columns') return {x:[null], y:[null], text:[null]};
        const x=[], y=[], text=[];
        for(let i=lowerBound(g.x, lo); i<g.x.length && g.x[i]<=hi; i++) {
          const px=g.x[i], value=g.raw[i];
          if (!Number.isFinite(px) || !Number.isFinite(value)) continue;
          const label=`position: ${px}<br>count: ${value}<br>frame: ${g.frame}`;
          x.push(px,px,null); y.push(0,value,null); text.push(label,label,null);
        }
        return x.length ? {x,y,text} : {x:[null],y:[null],text:[null]};
      });
      // One Plotly update avoids drawing the full old and new plot separately.
      await Plotly.restyle(elem, {
        x:lineData.map(v=>v.x).concat(glData.map(v=>v.x)),
        y:lineData.map(v=>v.y).concat(glData.map(v=>v.y)),
        text:entries.map(()=>null).concat(glData.map(v=>v.text)),
        fill:entries.map(()=>mode === 'area' ? 'tozeroy' : 'none').concat(entries.map(()=>'none')),
        fillcolor:entries.map(g=>g.smooth.fill).concat(entries.map(()=>'rgba(0,0,0,0)')),
        'line.color':entries.map(g=>mode === 'area' ? 'black' : g.smooth.color).concat(entries.map(g=>g.smooth.color)),
        'line.width':entries.map(()=>mode === 'area' ? 0.4 : 0.8).concat(entries.map(()=>width)),
        'line.simplify':lineData.map(v=>config.simplify !== false && mode !== 'columns' && v.x.length > 5000).concat(entries.map(()=>false)),
        visible:entries.map(g=>mode !== 'columns' && !hidden.has(g.frame)).concat(
          entries.map(g=>hidden.has(g.frame) ? 'legendonly' : true))
      }, indices);
      if (disposed) return;
      if (mode !== lastMode) {
        lastMode = mode;
        const status = document.getElementById(config.status_id);
        if (status) status.textContent = `Current: ${mode}, k=${mode === 'columns' ? 1 : 9}`;
        if (window.Shiny) Shiny.setInputValue(config.input_id, mode, {priority:'event'});
      }
    } finally {
      running = false;
      if (pending) { pending = false; requestAnimationFrame(update); }
    }
  };
  const onRelayout = event => {
    if (Object.keys(event).some(k => /^xaxis\d*\.(range|autorange)/.test(k))) requestAnimationFrame(update);
  };
  const onLegend = event => {
    const meta = elem.data?.[event.curveNumber]?.meta;
    if (!meta?.rc_columns_switch) return;
    const frame = String(meta.rc_columns_frame);
    if (hidden.has(frame)) hidden.delete(frame); else hidden.add(frame);
    lastKey = ''; requestAnimationFrame(update);
    return false;
  };
  elem.on('plotly_relayout', onRelayout);
  elem.on('plotly_legendclick', onLegend);
  elem._rcDynamicCleanup = () => { disposed = true; elem.removeListener('plotly_relayout', onRelayout); elem.removeListener('plotly_legendclick', onLegend); elem.__rcDynamicActive = false; };
  update();
};
