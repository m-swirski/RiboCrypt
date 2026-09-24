// Exercise the actual app renderer, without modifying any Plotly trace types.
const {chromium} = require('playwright');
const assert = require('node:assert/strict');
const url = process.argv[2] || 'http://127.0.0.1:7820/';

(async()=>{
  const browser=await chromium.launch({executablePath:process.env.CHROME || '/opt/google/chrome/chrome',
    headless:true,args:['--no-sandbox']});
  try {
    const context=await browser.newContext({viewport:{width:1440,height:1000},
      permissions:['clipboard-read','clipboard-write']});
    const page=await context.newPage();
    const errors=[];page.on('pageerror',error=>errors.push(error.message));
    await page.goto(url);
    await page.waitForFunction(()=>document.getElementById('browser-c')?.data?.length===21,
      null,{timeout:120000});
    await page.waitForTimeout(500);
    const initial=await page.evaluate(()=>{
      const el=document.getElementById('browser-c');
      return {range:el._fullLayout.xaxis.range,types:el.data.filter(t=>t.yaxis==='y4').map(t=>t.type),
        canvases:el.querySelectorAll('canvas').length,gene:document.getElementById('browser-gene').value};
    });
    assert.equal(initial.gene,'ATF4-ENSG00000128272');
    assert.deepEqual(initial.types,['scatter','scatter','scatter']);
    assert.equal(initial.canvases,0);
    const hover=await page.evaluate(()=>{
      const el=document.getElementById('browser-c');
      const trace=el.data.find(t=>t.yaxis==='y4' && t.mode==='lines' && t.line.color==='white');
      const rect=el.getBoundingClientRect();
      const x=el._fullLayout.xaxis, y=el._fullLayout.yaxis4;
      return {x:rect.x+x._offset+x.l2p(trace.x[0]),
        y:rect.y+y._offset+y.l2p((trace.y[0]+trace.y[1])/2)};
    });
    await page.mouse.move(hover.x,hover.y);
    await page.waitForFunction(()=>document.querySelector('#browser-c .hoverlayer')?.textContent.includes('Start codon'));
    await page.mouse.move(20,20);
    const copyBox=await page.getByText('Click to copy sequence',{exact:false}).first().boundingBox();
    await page.mouse.click(copyBox.x+copyBox.width/2,copyBox.y+copyBox.height/2);
    const copied=await page.evaluate(()=>navigator.clipboard.readText());
    assert.equal(copied.length,initial.range[1]-initial.range[0]);
    assert.match(copied,/^[ACGTN]+$/);
    await page.evaluate(()=>Plotly.relayout('browser-c',{'xaxis.range':[1000,1050]}));
    await page.waitForFunction(()=>document.getElementById('browser-c').data.some(t=>t.name==='sequence' && t.text?.length));
    await page.evaluate(()=>Plotly.relayout('browser-c',{'xaxis.autorange':true}));
    await page.waitForFunction(()=>{
      const el=document.getElementById('browser-c');
      return el.data.some(t=>t.name==='sequence_placeholder') && !el.data.some(t=>t.name==='sequence');
    });
    const reset=await page.evaluate(()=>document.getElementById('browser-c')._fullLayout.xaxis.range);
    assert.deepEqual(reset,initial.range);
    await page.screenshot({path:'/tmp/aa-svg-smoke-desktop.png'});
    await page.setViewportSize({width:390,height:844});
    await page.waitForTimeout(1000);
    assert.deepEqual(await page.evaluate(()=>document.getElementById('browser-c')._fullLayout.xaxis.range),initial.range);
    await page.screenshot({path:'/tmp/aa-svg-smoke-mobile.png',fullPage:true});
    assert.deepEqual(errors,[]);
    console.log(JSON.stringify({initial,hover:true,copyLength:copied.length,zoom:true,reset:true,resize:true,errors}));
  } finally {await browser.close();}
})().catch(error=>{console.error(error);process.exitCode=1;});
