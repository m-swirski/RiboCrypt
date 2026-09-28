// Exercise real Shiny downloads from both browser modules.
const {chromium}=require('playwright');
const assert=require('node:assert/strict');
const fs=require('node:fs');
const base=process.argv[2] || 'http://127.0.0.1:7821/';

(async()=>{
  const browser=await chromium.launch({executablePath:'/opt/google/chrome/chrome',
    headless:true,args:['--no-sandbox']});
  try {
    for(const observatory of [false,true]) {
      const context=await browser.newContext({viewport:{width:1440,height:1000}});
      const page=await context.newPage(), errors=[];
      page.on('pageerror',e=>errors.push(e.message));
      const prefix=observatory?'observatory-browser_obs-':'browser-';
      await page.goto(base+(observatory?'#Observatory':'#browser'));
      if(observatory) {
        await page.waitForFunction(()=>{
          const t=document.querySelector('#observatory-selector-libraries_data_table .dataTables_scrollBody table');
          return t && jQuery.fn.dataTable.isDataTable(t) && jQuery(t).DataTable().page.info().recordsTotal===3857;
        },null,{timeout:120000});
        await page.getByRole('button',{name:'Subset to page',exact:true}).click();
        await page.waitForFunction(()=>jQuery('#observatory-selector-libraries_data_table .dataTables_scrollBody table')
          .DataTable().page.info().recordsTotal===15);
        await page.locator('a[data-value="Browse"]').click();
        await page.waitForFunction(p=>document.getElementById(p+'tx')?.value,prefix);
        await page.locator('#'+prefix+'go').click();
      }
      await page.waitForFunction(id=>document.getElementById(id)?.data?.length>0,
        prefix+(observatory?'browser_plot':'c'),{timeout:120000});
      await page.locator('#'+prefix+'toggle_settings').click();
      await page.locator('#'+prefix+'floating_settings a[data-value="Settings"]').click();
      const button=page.locator('#'+prefix+'download_coverage');
      await button.scrollIntoViewIfNeeded();
      const [download]=await Promise.all([page.waitForEvent('download'),button.click()]);
      const file='/tmp/coverage-'+(observatory?'observatory':'browser')+'.csv';
      await download.saveAs(file);
      const traces=await page.evaluate(id=>document.getElementById(id).data.map(t=>({
        name:t.name,x:t.x,y:t.y,z:t.z,yaxis:t.yaxis,hovertemplate:t.hovertemplate,
        hoverinfo:t.hoverinfo,mode:t.mode,type:t.type
      })),prefix+(observatory?'browser_plot':'c'));
      fs.writeFileSync(file+'.plot.json',JSON.stringify(traces));
      assert.equal(await download.failure(),null);
      const csv=fs.readFileSync(file,'utf8');
      assert.ok(csv.startsWith('position,'));
      assert.ok(csv.split('\n').length>100);
      assert.match(download.suggestedFilename(),/^RiboCrypt_.+_coverage\.csv$/);
      assert.deepEqual(errors,[]);
      assert.deepEqual(await page.locator('.shiny-output-error').allTextContents(),[]);
      await page.screenshot({path:'/tmp/coverage-'+(observatory?'observatory':'browser')+'.png'});
      console.log(JSON.stringify({observatory,file,filename:download.suggestedFilename(),
        bytes:csv.length,header:csv.split('\n')[0],errors}));
      await context.close();
    }
  } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
