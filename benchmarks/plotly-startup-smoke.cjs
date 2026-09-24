// Verify early dependency registration, single loading and preserved WebGL zoom.
const {chromium} = require('playwright');
const assert = require('node:assert/strict');
const url = process.argv[2] || 'http://127.0.0.1:7821/';

(async()=>{
  const browser=await chromium.launch({executablePath:process.env.CHROME || '/opt/google/chrome/chrome',
    headless:true,args:['--no-sandbox']});
  try {
    const context=await browser.newContext({viewport:{width:1440,height:1000}});
    // The second page also checks an already populated browser HTTP cache.
    for (let visit=0;visit<2;visit++) {
      const page=await context.newPage();
      const errors=[];
      page.on('pageerror',error=>errors.push(error.message));
      await page.addInitScript(()=>{
        window.dependencyMarks={};
        const timer=setInterval(()=>{
          if(!window.jQuery) return;
          clearInterval(timer);
          jQuery(document).on('shiny:connected',()=>{
            window.dependencyMarks.readyAtConnection=typeof window.Plotly?.newPlot==='function';
          });
          jQuery(document).on('shiny:value',event=>{
            if(event.name==='browser-c') window.dependencyMarks.received=performance.now();
          });
        },1);
      });
      await page.goto(url);
      await page.waitForFunction(()=>document.getElementById('browser-c')?.data?.length===21,
        null,{timeout:120000});
      await page.evaluate(()=>Plotly.relayout('browser-c',{'xaxis.range':[1000,1050]}));
      await page.waitForFunction(()=>document.getElementById('browser-c').data
        .filter(t=>t.name==='sequence' && t.type==='scattergl').length===6);
      const state=await page.evaluate(()=>({marks:window.dependencyMarks,
        scripts:[...document.scripts].filter(s=>s.src.includes('/plotly-main-')).length,
        resources:performance.getEntriesByType('resource')
          .filter(r=>r.name.includes('/plotly-main-'))
          .map(r=>({start:r.startTime,end:r.responseEnd,initiator:r.initiatorType})),
        gene:document.getElementById('browser-gene').value,
        errors:[...document.querySelectorAll('.shiny-output-error')].map(e=>e.textContent)}));
      assert.equal(state.gene,'ATF4-ENSG00000128272');
      assert.equal(state.marks.readyAtConnection,true);
      assert.equal(state.scripts,1);
      assert.equal(state.resources.length,1);
      assert.equal(state.resources[0].initiator,'script');
      assert.ok(state.resources[0].end<state.marks.received);
      assert.deepEqual(state.errors,[]);
      assert.deepEqual(errors,[]);
      console.log(JSON.stringify({visit,...state,webglZoom:true,errors}));
      await page.close();
    }
  } finally {await browser.close();}
})().catch(error=>{console.error(error);process.exitCode=1;});
