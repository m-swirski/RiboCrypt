const {test}=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
const source=fs.readFileSync(process.env.RIBOCRYPT_TRANSLON_JS ||
  path.resolve(__dirname,'../../../inst/js/translon_table_controls.js'),'utf8');
function harness(checked=false) {
  const checkbox={checked},node={},events=new Map(),visible=[true,true,true,true,true];
  const calls=[];
  const $=element=>({off:event=>{events.delete(element===checkbox?event:'table:'+event);return $(element);},
    on:(event,fn)=>{events.set(element===checkbox?event:'table:'+event,fn);return $(element);}});
  const table={column:i=>({visible:()=>visible[i]}),table:()=>({node:()=>node})};
  table.columns=indices=>({visible:(value,redraw)=>{
    calls.push({indices,value,redraw});indices.forEach(i=>visible[i]=value);
  }});
  table.columns.adjust=()=>calls.push('adjust');
  const run=vm.runInNewContext(source,{$,document:{getElementById:()=>checkbox}});
  const bind=()=>run(table,'simplified',[3,4]);bind();
  return {visible,calls,events,bind,toggle:value=>{checkbox.checked=value;events.get('change.ribocryptSimplified')();}};
}
test('default full view does no work; toggles hide only excluded columns without drawing',()=>{
  const h=harness();assert.equal(h.calls.length,0);
  h.toggle(true);assert.deepEqual(h.visible,[true,true,true,false,false]);
  assert.equal(h.calls[0].redraw,false);
  const calls=h.calls.length;h.toggle(true);assert.equal(h.calls.length,calls);
  h.toggle(false);assert.deepEqual(h.visible,[true,true,true,true,true]);
});
test('checked state applies to new tables and old handlers are replaced and removed',()=>{
  const h=harness(true);assert.deepEqual(h.visible,[true,true,true,false,false]);
  h.bind();assert.equal(h.events.size,2);
  h.events.get('table:destroy.dt')();assert.equal(h.events.has('change.ribocryptSimplified'),false);
});
