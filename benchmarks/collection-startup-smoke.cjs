// Check the real app's independent module startup and state after tab switches.
const {chromium}=require('playwright');
const assert=require('node:assert/strict');
const base=process.argv[2] || 'http://127.0.0.1:7821/';
const prefix='observatory-selector-';
async function observatoryReady(page, rows=3857) {
  await page.waitForFunction(({prefix,rows})=>{
    const table=document.querySelector('#'+prefix+'libraries_data_table .dataTables_scrollBody table');
    return table && jQuery.fn.dataTable.isDataTable(table) &&
      jQuery(table).DataTable().page.info().recordsTotal===rows &&
      document.getElementById(prefix+'libraries_umap_plot')?.data?.length>0;
  },{prefix,rows},{timeout:120000});
}
async function megaReady(page) {
  await page.waitForFunction(()=>performance.getEntriesByType('resource')
    .some(r=>r.name.includes('/dataobj/browser_allsamp-gene')) &&
    document.getElementById('browser_allsamp-tx')?.value,
    null,{timeout:60000});
}
(async()=>{
  const browser=await chromium.launch({executablePath:process.env.CHROME || '/opt/google/chrome/chrome',
    headless:true,args:['--no-sandbox']});
  try {
    for(const first of ['Observatory','MegaBrowser']) {
      const context=await browser.newContext({viewport:{width:1440,height:1000}});
      const page=await context.newPage(), errors=[];
      page.on('pageerror',error=>errors.push(error.message));
      await page.goto(base+'#'+first);
      if(first==='Observatory') {
        await observatoryReady(page);
        await page.waitForTimeout(750);
        assert.equal(await page.evaluate(()=>performance.getEntriesByType('resource')
          .filter(r=>r.name.includes('/dataobj/browser_allsamp-')).length),0);
      } else {
        await megaReady(page);
        await page.waitForTimeout(750);
        assert.equal(await page.evaluate(()=>performance.getEntriesByType('resource')
          .filter(r=>r.name.includes('/dataobj/observatory-')).length),0);
        await page.locator('a[data-value="Observatory"]').click();
        await observatoryReady(page);
      }
      await page.getByRole('button',{name:'Subset to page',exact:true}).click();
      await observatoryReady(page,15);
      await page.locator('a[data-value="Browse"]').click();
      await page.waitForFunction(()=>document.getElementById('observatory-browser_obs-tx')?.value,
        null,{timeout:60000});
      await page.locator('#observatory-browser_obs-go').click();
      await page.waitForFunction(()=>document.getElementById('observatory-browser_obs-browser_plot')?._fullData?.length>0,
        null,{timeout:120000});
      await page.locator('a[data-value="Select libraries"]').click();
      await observatoryReady(page,15);
      await page.locator('a[data-value="MegaBrowser"]').first().click();
      await megaReady(page);
      const gene=await page.locator('#browser_allsamp-gene').inputValue();
      const tx=await page.locator('#browser_allsamp-tx').inputValue();
      assert.ok(gene && tx);
      await page.locator('a[data-value="Observatory"]').click();
      await observatoryReady(page,15);
      await page.getByRole('button',{name:'Remove subset',exact:true}).click();
      await observatoryReady(page);
      await page.locator('a[data-value="MegaBrowser"]').first().click();
      assert.equal(await page.locator('#browser_allsamp-gene').inputValue(),gene);
      assert.equal(await page.locator('#browser_allsamp-tx').inputValue(),tx);
      assert.deepEqual(await page.locator('.shiny-output-error').allTextContents(),[]);
      assert.deepEqual(errors,[]);
      console.log(JSON.stringify({first,hiddenModuleDeferred:true,
        selectedPagePlotted:true,subsetRetained:true,gene,tx,errors}));
      await context.close();
    }
  } finally {await browser.close();}
})().catch(error=>{console.error(error);process.exitCode=1;});
