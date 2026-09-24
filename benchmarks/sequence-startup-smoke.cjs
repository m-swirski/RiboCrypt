// Check initial URL zoom against a source-loaded app with the default human data.
const {chromium} = require('playwright');
const assert = require('node:assert/strict');
const url = new URL(process.argv[2] || 'http://127.0.0.1:7821/');
Object.entries({dff:'human_all_merged_l50',gene:'ATF4-ENSG00000128272',
  tx:'ENST00000674920',zoom_range:'1000:1050',go:'TRUE'})
  .forEach(([key,value])=>url.searchParams.set(key,value));

(async()=>{
  const browser = await chromium.launch({executablePath:process.env.CHROME || '/opt/google/chrome/chrome',
    headless:true,args:['--no-sandbox']});
  try {
    const page = await browser.newPage({viewport:{width:1440,height:1000}});
    const errors=[];
    page.on('pageerror',error=>errors.push(error.message));
    await page.goto(url.href);
    await page.waitForFunction(()=>{
      const el=document.getElementById('browser-c');
      return el?.data?.filter(t=>t.name==='sequence').length===6 &&
        !document.documentElement.classList.contains('shiny-busy');
    },null,{timeout:120000});
    const initial=await page.evaluate(()=>{
      const el=document.getElementById('browser-c');
      return {range:el._fullLayout.xaxis.range,
        placeholders:el.data.filter(t=>t.name==='sequence_placeholder').length,
        populated:el.data.filter(t=>t.name==='sequence').every(t=>t.text.length>0),
        gene:document.getElementById('browser-gene').value};
    });
    assert.deepEqual(initial.range,[1000,1050]);
    assert.equal(initial.gene,'ATF4-ENSG00000128272');
    assert.equal(initial.placeholders,0);
    assert.equal(initial.populated,true);
    await page.evaluate(()=>Plotly.relayout('browser-c',{'xaxis.range':[1,2294]}));
    await page.waitForFunction(()=>{
      const el=document.getElementById('browser-c');
      return el.data.filter(t=>t.name==='sequence_placeholder').length===1 &&
        !el.data.some(t=>t.name==='sequence');
    });
    assert.deepEqual(errors,[]);
    console.log(JSON.stringify({initial,fullRangeRestored:true,errors}));
  } finally {await browser.close();}
})().catch(error=>{console.error(error);process.exitCode=1;});
