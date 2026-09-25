const {test}=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs'), path=require('node:path'), vm=require('node:vm');
const source=fs.readFileSync(process.env.RIBOCRYPT_UMAP_JS ||
  path.resolve(__dirname,'../../../inst/js/umap_plot_extension.js'),'utf8');

function harness({mode='select',selection,region=[],fullSelection,fullRegion=[]}={}) {
  const calls=[],inputs=[],handlers={},events={};
  const elem={data:[{text:['Run: A<br>Tissue: kidney','Run: B'],selectedpoints:selection}],
    layout:{selections:region},_fullLayout:{dragmode:mode,selections:fullRegion},
    _fullData:[{selectedpoints:fullSelection}],
    on:(name,fn)=>{events[name]=fn;},addEventListener:()=>{},querySelectorAll:()=>[]};
  const window={requestAnimationFrame:fn=>fn()};
  const Plotly={
    relayout:(_,patch)=>{calls.push('relayout');Object.assign(elem._fullLayout,patch);},
    react:(_,data,layout)=>{
      calls.push('react');elem.data=data;elem.layout=layout;
      elem._fullData=data.map(t=>({...t}));elem._fullLayout.selections=layout.selections;
      return Promise.resolve();
    }
  };
  const run=vm.runInNewContext(source,{window,Plotly,document:{body:{contains:()=>true}},
    Shiny:{setInputValue:(id,value)=>inputs.push({id,value}),
      addCustomMessageHandler:(name,handler)=>{handlers[name]=handler;}}});
  run(elem,null,'umap_selection');
  return {elem,calls,inputs,events,reset:()=>handlers.librariesActiveSelectionReset({target:'umap_selection'}),
    select:runs=>handlers.librariesActiveSelectionChanged({target:'umap_selection',runs})};
}

test('initial select mode and repeated empty resets require no Plotly updates',()=>{
  const h=harness();h.reset();h.reset();
  assert.deepEqual(h.calls,[]);
  assert.equal(h.inputs[0].id,'umap_selection_ready');
  h.events.plotly_deselect();
  assert.equal(h.inputs.at(-1).id,'umap_selection');
  assert.equal(h.inputs.at(-1).value.length,0);
});

test('a different drag mode is still corrected',()=>{
  const h=harness({mode:'pan'});
  assert.deepEqual(h.calls,['relayout']);
  assert.equal(h.elem._fullLayout.dragmode,'select');
});

test('nonempty and empty selectedpoints both reset, but only once',()=>{
  for(const selection of [[0],[]]) {
    const h=harness({selection});h.reset();h.reset();
    assert.deepEqual(h.calls,['react']);
    assert.equal(h.elem.data[0].selectedpoints,null);
  }
});

test('native fullData selection and drawn regions must not be skipped',()=>{
  for(const initial of [{fullSelection:[0]},{region:[{x0:1,x1:2}]},{fullRegion:[{x0:1,x1:2}]}]) {
    const h=harness(initial);h.reset();h.reset();
    assert.deepEqual(h.calls,['react']);
    assert.equal(h.elem.layout.selections.length,0);
  }
});

test('programmatic selection, reset and subsequent user selection still propagate',async()=>{
  const h=harness();h.select(['B']);
  assert.equal(h.elem.data[0].selectedpoints[0],1);
  await Promise.resolve();h.reset();await Promise.resolve();
  assert.deepEqual(h.calls,['react','react']);
  h.events.plotly_selected({points:[{text:'Run: A<br>Tissue: kidney'}]});
  assert.equal(h.inputs.at(-1).value[0],'A');
});
