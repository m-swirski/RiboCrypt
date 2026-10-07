const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const {test} = require('node:test');

test('Shiny updates tolerate document targets and resolve the output by name', () => {
  const handlers = {};
  let resized = 0;
  const workspace = {clientWidth:390, querySelectorAll:()=>[plot]};
  const plot = {offsetWidth:300, offsetHeight:300, _fullLayout:{margin:{l:30,r:10}},
    closest: selector => selector === '.mega-workspace' ? workspace : null};
  const document = {addEventListener:()=>{}, getElementById:id=>id==='heatmap'?plot:null};
  vm.runInNewContext(fs.readFileSync(process.env.RIBOCRYPT_MEGA_LAYOUT_JS, 'utf8'), {
    window:{addEventListener:()=>{}, Plotly:{Plots:{resize:()=>resized++}}}, document,
    jQuery:()=>({on:(name, handler)=>handlers[name]=handler}), WeakMap,
    setTimeout:callback=>callback(), requestAnimationFrame:callback=>{callback();return 1;},
    cancelAnimationFrame:()=>{}
  });
  assert.doesNotThrow(()=>handlers['shiny:value']({target:document}));
  assert.equal(resized, 0);
  handlers['shiny:value']({target:document, name:'heatmap'});
  assert.equal(resized, 1);
  handlers['shiny:value']({target:plot});
  assert.equal(resized, 2);
});
