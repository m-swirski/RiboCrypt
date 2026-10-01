const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');
const source = fs.readFileSync(process.env.RIBOCRYPT_LOCAL_Y_JS ||
  path.resolve(__dirname, '../../../inst/js/browser_local_y_zoom.js'), 'utf8');

function fixture(options = {}) {
  const callbacks = {}, frames = new Map(), changes = [];
  let next = 0;
  const checkbox = {checked: options.checked !== false};
  const axis = () => ({range: [0, 1100], tickmode: 'array', tickvals: [0, 500, 1000], autorange: false});
  const elem = {id: 'browser-c', data: options.traces || [
    {x: [1, 2, 3, 4], y: [1000, 10, 5, 0], mode: 'lines', yaxis: 'y'},
    {x: [1, 2, 3, 4], y: [500, 20, 0, 0], type: 'bar', yaxis: 'y2'},
    {x: [1, 2, 3, 4], y: [1e6, 1e6, 1e6, 1e6], yaxis: 'y3'}],
    _fullLayout: {xaxis: {range: [1, 4]}, yaxis: axis(), yaxis2: axis(), yaxis3: axis()},
    on: (name, fn) => (callbacks[name] ||= new Set()).add(fn),
    removeListener: (name, fn) => callbacks[name].delete(fn)};
  const emit = (name, event = {}) => [...(callbacks[name] || [])].forEach(fn => fn(event));
  const run = vm.runInNewContext(source, {
    document: {getElementById: () => checkbox},
    $: () => ({on: (_, fn) => checkbox.change = fn, off: () => delete checkbox.change}),
    requestAnimationFrame: fn => {frames.set(++next, fn); return next;},
    cancelAnimationFrame: id => frames.delete(id),
    Plotly: {relayout: (_, patch) => {
      changes.push(patch);
      for (const [key, value] of Object.entries(patch)) {
        const [axis, field] = key.split('.'); elem._fullLayout[axis][field] = value;
      }
      emit('plotly_relayout', patch);
    }}
  });
  const config = {axes: ['yaxis', 'yaxis2'], input_id: 'browser-local_y_max',
    full_range: [1, 4], manual: !!options.manual};
  const flush = () => {const pending = [...frames.values()]; frames.clear(); pending.forEach(fn => fn());};
  run(elem, null, config); flush();
  const zoom = range => {elem._fullLayout.xaxis.range = range; emit('plotly_relayout', {'xaxis.range': range}); flush();};
  return {elem, checkbox, changes, zoom, flush, emit, frames, callbacks, run, config};
}

test('zoom and pan scale each coverage track, reset restores ticks and other axes', () => {
  const f = fixture();
  assert.equal(f.changes.length, 0);
  f.zoom([2, 3]);
  assert.equal(f.elem._fullLayout.yaxis.range[1], 10.5);
  assert.equal(f.elem._fullLayout.yaxis2.range[1], 21);
  assert.equal(f.elem._fullLayout.yaxis3.range[1], 1100);
  assert.equal(f.frames.size, 0, 'Y relayout must not recurse');
  f.zoom([3, 4]);
  assert.equal(f.elem._fullLayout.yaxis.range[1], 5.25);
  assert.equal(f.elem._fullLayout.yaxis2.range[1], 1);
  f.zoom([1, 4]);
  assert.equal(f.elem._fullLayout.yaxis.range[1], 1100);
  assert.equal(f.elem._fullLayout.yaxis.tickmode, 'array');
});

test('live checkbox, manual limits, and rerender cleanup', () => {
  const f = fixture({checked: false}); f.zoom([2, 3]);
  assert.equal(f.changes.length, 0);
  f.checkbox.checked = true; f.checkbox.change(); f.flush();
  assert.equal(f.elem._fullLayout.yaxis.range[1], 10.5);
  f.checkbox.checked = false; f.checkbox.change(); f.flush();
  assert.equal(f.elem._fullLayout.yaxis.range[1], 1100);
  f.run(f.elem, null, f.config); f.flush();
  assert.equal(f.callbacks.plotly_relayout.size, 1);
  const manual = fixture({manual: true}); manual.zoom([2, 3]);
  assert.equal(manual.changes.length, 0);
});

test('boundary interpolation, filled polygons, hidden and non-finite values', () => {
  const f = fixture({traces: [
    {x: [1, 4, 4, 1], y: [0, 30, 0, 0], fill: 'toself'},
    {x: [2], y: [1e9], visible: 'legendonly'},
    {x: [2, 3], y: [null, NaN]},
    {x: [2, 3], y: [Infinity, -Infinity]}
  ]});
  f.zoom([2, 3]); assert.equal(f.elem._fullLayout.yaxis.range[1], 21);
  f.zoom([3, 2]); assert.equal(f.elem._fullLayout.yaxis.range[1], 21);
  f.zoom([5, 6]); assert.equal(f.elem._fullLayout.yaxis.range[1], 1);
});

test('guide lines cannot inflate local limits and stacked totals are not clipped', () => {
  const f = fixture({traces: [
    {x: [1, 2, 3, 4], y: [1, 2, 3, 4], mode: 'lines', stackgroup: 'coverage'},
    {x: [2, 2], y: [0, 1e9], mode: 'lines', hoverinfo: ['skip', 'skip']}
  ]});
  f.elem.calcdata = [[{x: 1, y: 10}, {x: 2, y: 20}, {x: 3, y: 30}, {x: 4, y: 40}]];
  f.zoom([2, 3]); assert.equal(f.elem._fullLayout.yaxis.range[1], 31.5);
});
