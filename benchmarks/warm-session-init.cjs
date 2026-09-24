// Fresh browser contexts against a running app; the first visit primes its cache.
// NODE_PATH=/path/to/node_modules node benchmarks/warm-session-init.cjs URL[,URL] [visits]
const { chromium } = require('playwright');
const fs = require('node:fs');
const assert = require('node:assert/strict');

(async () => {
  const browser = await chromium.launch({executablePath:process.env.CHROME || '/opt/google/chrome/chrome',
    headless:true,args:['--no-sandbox']});
  const results = [];
  const urls = process.argv[2].split(',');
  const visits = Number(process.argv[3] || 6);
  try {
    for (let trial=0; trial<visits*urls.length; trial++) {
      const visit = Math.floor(trial / urls.length);
      const url = urls[trial % urls.length];
      const context = await browser.newContext({viewport:{width:1440,height:1000}});
      const page = await context.newPage();
      const errors = [];
      page.on('pageerror',error=>errors.push(error.message));
      await page.addInitScript(() => {
        window.loadMarks = {};
        const timer = setInterval(() => {
          if (!window.jQuery) return;
          clearInterval(timer);
          jQuery(document).on('shiny:connected',()=>{window.loadMarks.connected=performance.now();});
          jQuery(document).on('shiny:value',event=>{
            if(event.name === 'browser-c') window.loadMarks.plotReceived=performance.now();
          });
        }, 1);
      });
      await page.goto(url,{waitUntil:'domcontentloaded'});
      await page.waitForFunction(() => {
        const plot = document.getElementById('browser-c');
        if(!plot?._fullData?.length || !plot.querySelector('.main-svg')) return false;
        window.loadMarks.visible ||= performance.now();
        return true;
      }, null, {timeout:120000,polling:25});
      await page.waitForFunction(() => {
        if(document.documentElement.classList.contains('shiny-busy')) return false;
        window.loadMarks.idle ||= performance.now(); return true;
      }, null, {timeout:120000,polling:25});
      const result = await page.evaluate(() => ({...window.loadMarks,
        responseEnd:performance.getEntriesByType('navigation')[0].responseEnd,
        traces:document.getElementById('browser-c')._fullData.length,
        gene:document.getElementById('browser-gene').value,
        tx:document.getElementById('browser-tx').value,
        errors:[...document.querySelectorAll('.shiny-output-error')].map(e=>e.textContent)}));
      assert.equal(result.gene,'ATF4-ENSG00000128272');
      assert.equal(result.tx,'ENST00000674920');
      assert.equal(result.traces,21);
      assert.deepEqual(errors,[]); assert.deepEqual(result.errors,[]);
      results.push({url,visit,...result}); console.log(JSON.stringify(results.at(-1)));
      if (visit === 0) await page.screenshot({path:'/tmp/warm-init-'+new URL(url).port+'.png'});
      await page.waitForTimeout(1000);
      await context.close();
    }
    if(process.env.RESULTS) fs.writeFileSync(process.env.RESULTS,JSON.stringify(results,null,2));
  } finally {await browser.close();}
})().catch(error=>{console.error(error);process.exitCode=1;});
