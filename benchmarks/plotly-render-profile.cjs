// Profile public Plotly calls and widget setup without changing application code.
// NODE_PATH=/path/to/node_modules node benchmarks/plotly-render-profile.cjs [URL]
const { chromium } = require('playwright');
const fs = require('node:fs');
const assert = require('node:assert/strict');
const base = process.argv[2] || 'http://127.0.0.1:7820/';

function instrument({mode}) {
  const record = window.plotProfile = {mode,calls:[],tasks:[],marks:{}};
  const stack = [];
  function wrap(obj, key, prefix) {
    const original = obj[key];
    if(typeof original !== 'function') return;
    obj[key] = function(...args) {
      const el = typeof args[0] === 'string' ? document.getElementById(args[0]) : args[0];
      if(el?.id !== 'browser-c') return original.apply(this,args);
      const call = {method:prefix+key,start:performance.now(),parent:stack.at(-1)?.index,
        index:record.calls.length,traceCount:el.data?.length,
        patch:!['newPlot','react','renderValue'].includes(key) ? args[1] : undefined,
        caller:new Error().stack.split('\n').slice(2,6).join('\n')};
      record.calls.push(call); stack.push(call);
      if(key === 'newPlot') {
        const data = args[1]?.data || args[1];
        // Diagnostic only: change the three small bottom-panel GL traces in memory.
        if(mode === 'svg-bottom') data.forEach(t=>{
          if(t.type === 'scattergl' && t.yaxis === 'y4') t.type='scatter';
        });
        record.initialTraces = data.map(t=>({type:t.type,mode:t.mode,name:t.name,
          visible:t.visible,points:t.x?.length,yaxis:t.yaxis}));
      }
      try {
        const value = original.apply(this,args);
        call.syncEnd = performance.now();
        if(value?.then) value.then(()=>{call.end=performance.now();},e=>{call.error=String(e);});
        else call.end=call.syncEnd;
        return value;
      } finally {stack.pop();}
    };
  }
  function interceptGlobal(name, install) {
    let value;
    Object.defineProperty(window,name,{configurable:true,get:()=>value,set:next=>{
      value=next; install(next);
    }});
  }
  interceptGlobal('Plotly',value=>{
    ['newPlot','react','relayout','restyle','addTraces','deleteTraces','redraw'].forEach(key=>wrap(value,key,'Plotly.'));
  });
  interceptGlobal('HTMLWidgets',value=>{
    let widget;
    Object.defineProperty(value,'widget',{configurable:true,get:()=>widget,set:original=>{
      widget=function(definition){
        if(definition.name === 'plotly') {
          ['initialize','renderValue','resize'].forEach(key=>wrap(definition,key,'widget.'));
        }
        return original.apply(this,arguments);
      };
    }});
  });
  new PerformanceObserver(list=>{
    for(const entry of list.getEntries()) record.tasks.push({start:entry.startTime,duration:entry.duration});
  }).observe({type:'longtask',buffered:true});
  const timer=setInterval(()=>{
    if(!window.jQuery) return;
    clearInterval(timer);
    jQuery(document).on('shiny:connected',()=>{record.marks.connected=performance.now();});
    jQuery(document).on('shiny:value',event=>{
      if(event.name === 'browser-c') record.marks.received=performance.now();
    });
  },1);
}

(async()=>{
  const browser=await chromium.launch({executablePath:process.env.CHROME || '/opt/google/chrome/chrome',
    headless:true,args:['--no-sandbox']});
  const results=[];
  try {
    const visits=Number(process.env.VISITS || 5);
    for(let visit=0;visit<visits;visit++) {
      const context=await browser.newContext({viewport:{width:1440,height:1000}});
      const page=await context.newPage();
      const errors=[];page.on('pageerror',e=>errors.push(e.message));
      const mode=process.env.COMPARE_SVG && visit>0 && visit%2===0 ? 'svg-bottom' : 'baseline';
      await page.addInitScript(instrument,{mode});
      const cdp=await context.newCDPSession(page);
      await cdp.send('Performance.enable');
      if(visit===visits-1) {await cdp.send('Profiler.enable');await cdp.send('Profiler.start');}
      await page.goto(base);
      await page.waitForFunction(()=>{
        const el=document.getElementById('browser-c');
        if(!el?._fullData?.length || !el.querySelector('.main-svg')) return false;
        window.plotProfile.marks.visible ||= performance.now();
        if(document.documentElement.classList.contains('shiny-busy')) return false;
        window.plotProfile.marks.idle ||= performance.now();return true;
      },null,{timeout:120000,polling:25});
      await page.waitForTimeout(1500);
      const result=await page.evaluate(()=>({...window.plotProfile,
        resources:performance.getEntriesByType('resource').map(r=>({name:r.name,
          start:r.startTime,end:r.responseEnd,duration:r.duration,bytes:r.transferSize})),
        finalTraces:document.getElementById('browser-c').data.map(t=>({type:t.type,mode:t.mode,
          name:t.name,visible:t.visible,points:t.x?.length})),
        gene:document.getElementById('browser-gene').value,
        tx:document.getElementById('browser-tx').value,
        canvasCount:document.getElementById('browser-c').querySelectorAll('canvas').length}));
      assert.equal(result.gene,'ATF4-ENSG00000128272');
      assert.equal(result.finalTraces.length,21);assert.deepEqual(errors,[]);
      if(visit===visits-1) {
        result.cdpMetrics=(await cdp.send('Performance.getMetrics')).metrics;
        const {profile}=await cdp.send('Profiler.stop');
        fs.writeFileSync('/tmp/plotly-render.cpuprofile',JSON.stringify(profile));
        await page.screenshot({path:'/tmp/plotly-render-profile.png'});
      }
      results.push({visit,errors,...result});
      console.log(JSON.stringify({visit,mode,marks:result.marks,calls:result.calls.map(c=>({
        method:c.method,start:c.start,sync:c.syncEnd-c.start,total:c.end-c.start,parent:c.parent,patch:c.patch}))}));
      await context.close();
    }
    fs.writeFileSync('/tmp/plotly-render-profile.json',JSON.stringify(results,null,2));
  } finally {await browser.close();}
})().catch(error=>{console.error(error);process.exitCode=1;});
