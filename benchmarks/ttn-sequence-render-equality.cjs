// Compare the browser's actual pre/post compact sequence traces on zoom.
// NODE_PATH=/tmp/ribocrypt_browser_bench/node_modules PLAYWRIGHT_BROWSERS_PATH=/tmp/ribocrypt_browser_bench/browsers node benchmarks/ttn-sequence-render-equality.cjs
const {chromium}=require('playwright'), assert=require('node:assert/strict');
(async()=>{const browser=await chromium.launch({headless:true,args:['--no-sandbox','--disable-dev-shm-usage','--enable-unsafe-swiftshader']});
async function get(mode,old){const page=await browser.newPage({viewport:{width:1400,height:900}}),errors=[];page.on('pageerror',e=>errors.push(e.message));
 const prefix=old?'/tmp/ribocrypt_ttn_before_compact/':'/tmp/';await page.goto('file://'+prefix+'ribocrypt_ttn_'+mode+'.html');
 await page.waitForFunction(()=>document.querySelector('.js-plotly-plot')?._fullData?.length,{timeout:120000});
 await page.evaluate(()=>Plotly.relayout(document.querySelector('.js-plotly-plot'),{'xaxis.range':[40000,40100]}));
 await page.waitForFunction(()=>document.querySelector('.js-plotly-plot').data.filter(t=>t.name==='sequence').length===6,{timeout:120000});
 const data=await page.evaluate(()=>document.querySelector('.js-plotly-plot').data.filter(t=>t.name==='sequence').map(t=>({x:Array.from(t.x),y:Array.from(t.y),text:Array.from(t.text),xaxis:t.xaxis,yaxis:t.yaxis,type:t.type})));
 assert.deepEqual(errors,[]);await page.close();return data}
for(const mode of ['lines','columns','area']){const before=await get(mode,true),after=await get(mode,false);assert.deepEqual(after,before);console.log(mode,'zoom sequence traces identical:',before.length,'letters',before.reduce((n,t)=>n+t.x.length,0))}
await browser.close()})().catch(e=>{console.error(e);process.exit(1)});
