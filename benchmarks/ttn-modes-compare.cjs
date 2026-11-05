// Alternate old/new standalone widgets; backend/network and genuine GPU are excluded.
const {chromium}=require('playwright');
(async()=>{const browser=await chromium.launch({headless:true,args:['--no-sandbox','--disable-dev-shm-usage','--enable-unsafe-swiftshader']});
for (const mode of ['lines','columns','area']) for(let i=0;i<3;i++) for(const version of (i%2?['new','old']:['old','new'])) {
 const ctx=await browser.newContext({viewport:{width:1400,height:900}});const page=await ctx.newPage();let errors=[];page.on('pageerror',e=>errors.push(e.message));
 await page.addInitScript(()=>{window.measure=[];let value;Object.defineProperty(window,'Plotly',{configurable:true,get(){return value},set(v){value=v;let f=v.newPlot;if(typeof f==='function')v.newPlot=function(...a){let q={start:performance.now()};window.measure.push(q);let r=f.apply(this,a);if(r?.then)r.then(()=>q.end=performance.now());else q.end=performance.now();return r}}})});
 let prefix=version==='old'?'/tmp/ribocrypt_ttn_before_compact/':'/tmp/';await page.goto('file://'+prefix+'ribocrypt_ttn_'+mode+'.html',{waitUntil:'load'});await page.waitForFunction(()=>window.measure?.[0]?.end,{timeout:120000});
 const r=await page.evaluate(()=>({call:window.measure[0],nav:performance.getEntriesByType('navigation')[0].loadEventEnd,types:document.querySelector('.js-plotly-plot')._fullData.map(x=>x.type)}));
 if(errors.length)throw Error(errors.join('; '));console.log(JSON.stringify({mode,version,iteration:i,draw:r.call.end-r.call.start,load:r.nav,types:r.types}));await ctx.close();
}await browser.close()})().catch(e=>{console.error(e);process.exit(1)});
