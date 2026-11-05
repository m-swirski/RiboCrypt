const {test} = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const source = fs.readFileSync(process.env.RIBOCRYPT_SEQUENCE_JS ||
  path.resolve(__dirname,'../../../inst/js/render_on_zoom.js'),'utf8');

function harness(range=[1,1000], axisOverrides={}, existing=[], suppliedData=null) {
  const calls=[],handlers={},copied=[],inputs=[];
  const axis={range,showticklabels:true,ticks:'outside',showline:false,showgrid:false,zeroline:false,...axisOverrides};
  const elem={data:[{name:'coverage',x:[1,1000]},...existing],_fullLayout:{xaxis:axis},
    on:(event,handler)=>{handlers[event]=handler;}};
  const data=suppliedData || {sequence:'ACGT'.repeat(250),input_id:'plot_copy',
    traces:Array.from({length:6},(_,index)=>{
      const x=Array.from({length:333},(_,i)=>i*3+index%3+1);
      return {x,y:x.map(()=>0),text:x.map(()=>'A'),color:'red',
        xaxis:'x',yaxis:index<3?'y2':'y4',distance:250};
    })};
  const Plotly={
    addTraces:(_,traces)=>{calls.push({method:'add',traces});elem.data.push(...traces);return Promise.resolve();},
    deleteTraces:(_,indices)=>{calls.push({method:'delete',indices});elem.data=elem.data.filter((t,i)=>!indices.includes(i));return Promise.resolve();},
    restyle:(_,patch,indices)=>{calls.push({method:'restyle',patch,indices});indices.forEach((index,i)=>{elem.data[index].x=patch.x[i];});return Promise.resolve();},
    relayout:(_,patch)=>{calls.push({method:'relayout',patch});Object.entries(patch).forEach(([key,value])=>{axis[key.split('.')[1]]=value;});return Promise.resolve();}
  };
  const run=vm.runInNewContext(source,{Plotly,console:{log:()=>{},error:()=>{}},
    Shiny:{setInputValue:(id,value)=>inputs.push({id,value})},
    navigator:{clipboard:{writeText:value=>{copied.push(value);return Promise.resolve();}}},
    document:{createElement:()=>({style:{}}),body:{appendChild:()=>{},removeChild:()=>{}}},
    setTimeout:()=>{}});
  run(elem,null,data);
  return {calls,elem,data,handlers,copied,inputs,move:(next,event={})=>{
    axis.range=next;handlers.plotly_relayout(event);
  }};
}

test('full-range startup adds one centered placeholder and no redundant updates',()=>{
  const h=harness();
  assert.deepEqual(h.calls.map(c=>c.method),['add']);
  const placeholder=h.elem.data.find(t=>t.name==='sequence_placeholder');
  assert.equal(placeholder.x[0],500.5);
  h.move([1,1000],{width:600});
  assert.equal(h.calls.length,1);
});

test('compact and legacy payloads create identical visible sequence traces',()=>{
  const sequence='ACGT'.repeat(301);
  const translations=['MATX*'.repeat(80),'RX*MT'.repeat(80),'TMRX*'.repeat(80)];
  const aa_codes={M:'Met',A:'Ala',T:'Thr',X:'Xaa','*':'*',R:'Arg'};
  for (const three_letters of [false,true]) {
    const common={sequence,translations,three_letters,aa_codes,input_id:'plot_copy'};
    const traces=Array.from({length:6},(_,i)=>{
      const frame=i%3, kind=i<3?'nt':'aa';
      const length=kind==='nt'?Math.max(0,Math.floor((sequence.length-frame-1)/3)+1):translations[frame].length;
      return {kind,frame,length,compact:true,y:kind==='nt'?0.5:2.4-(frame+1),
        color:kind==='nt'?['red','green','blue'][frame]:'grey45',xaxis:'x',yaxis:kind==='nt'?'y2':'y4',distance:250};
    });
    const legacy=traces.map(t=>{
      const x=Array.from({length:t.length},(_,i)=>t.frame+1+3*i);
      const text=x.map((pos,i)=>t.kind==='nt'?sequence[pos-1]:
        (three_letters?aa_codes[translations[t.frame][i]]:translations[t.frame][i]));
      return {x,text,y:x.map(()=>t.y),color:t.color,xaxis:t.xaxis,yaxis:t.yaxis,distance:t.distance};
    });
    for(const range of [[1,50],[299,345],[400,550],[1100,1190]]){
      const a=harness(range,{},[],{...common,traces:legacy});
      const b=harness(range,{},[],{...common,traces});
      assert.deepEqual(JSON.parse(JSON.stringify(b.calls)),JSON.parse(JSON.stringify(a.calls)));
    }
  }
});

test('initial zoom displays letters directly instead of a full-gene placeholder',()=>{
  const h=harness([100,150]);
  assert.deepEqual(h.calls.map(c=>c.method),['add']);
  assert.equal(h.elem.data.filter(t=>t.name==='sequence').length,6);
  assert.ok(!h.elem.data.some(t=>t.name==='sequence_placeholder'));
});

test('axis styling is applied once and only for changed properties',()=>{
  const h=harness([1,1000],{ticks:'',showline:true});
  const updates=h.calls.filter(c=>c.method==='relayout');
  assert.equal(updates.length,1);
  assert.deepEqual(JSON.parse(JSON.stringify(updates[0].patch)),{'xaxis.ticks':'outside','xaxis.showline':false});
  h.move([100,900]);
  assert.equal(h.calls.filter(c=>c.method==='relayout').length,1);
});

test('panning moves the placeholder only when its center changes',()=>{
  const h=harness();
  h.move([100,900]);
  assert.equal(h.calls.filter(c=>c.method==='restyle').length,1);
  h.move([101,899]);
  assert.equal(h.calls.filter(c=>c.method==='restyle').length,1);
  assert.equal(h.elem.data.find(t=>t.name==='sequence_placeholder').x[0],500);
});

test('zoom and reset retain one representation and do not restyle the axes',()=>{
  const h=harness();
  h.move([100,150],{'xaxis.range[0]':100,'xaxis.range[1]':150});
  assert.equal(h.elem.data.filter(t=>t.name==='sequence').length,6);
  assert.ok(!h.elem.data.some(t=>t.name==='sequence_placeholder'));
  h.move([120,170],{'xaxis.range':[120,170]});
  assert.equal(h.elem.data.filter(t=>t.name==='sequence').length,6);
  h.move([1,1000],{'xaxis.autorange':true});
  assert.equal(h.elem.data.filter(t=>t.name==='sequence_placeholder').length,1);
  assert.ok(!h.elem.data.some(t=>t.name==='sequence'));
  assert.equal(h.calls.filter(c=>c.method==='relayout').length,0);
});

test('an already initialized placeholder is not duplicated',()=>{
  const h=harness([1,1000],{},[{name:'sequence_placeholder',x:[500.5]}]);
  assert.equal(h.calls.length,0);
  assert.equal(h.elem.data.filter(t=>t.name==='sequence_placeholder').length,1);
});

test('sequence copying still uses the current visible range',async()=>{
  const h=harness();
  h.move([100,900]);
  const placeholder=h.elem.data.find(t=>t.name==='sequence_placeholder');
  h.handlers.plotly_click({points:[{data:placeholder}]});
  await Promise.resolve();
  assert.equal(h.copied[0],h.data.sequence.slice(99,900));
  assert.deepEqual(h.inputs,[{id:'plot_copy',value:h.copied[0]}]);
});

test('copy includes the first and last visible bases and handles reversed ranges',async()=>{
  for (const range of [[1,1000],[100,150],[150,100],[100.2,150.8]]) {
    const h=harness(range);
    const target=h.elem.data.find(t=>t.name==='sequence_placeholder') ||
      h.elem.data.find(t=>t.name==='sequence' && t.yaxis==='y2');
    h.handlers.plotly_click({points:[{data:target}]});
    await Promise.resolve();
    const start=Math.max(0,Math.ceil(Math.min(...range))-1);
    const end=Math.min(h.data.sequence.length,Math.floor(Math.max(...range)));
    assert.equal(h.copied[0],h.data.sequence.slice(start,end));
  }
});
