const {chromium}=require('playwright');
const assert=require('node:assert/strict');
const base=process.argv[2] || 'http://127.0.0.1:7821/';
(async()=>{
  const browser=await chromium.launch({executablePath:'/opt/google/chrome/chrome',
    headless:true,args:['--no-sandbox']});
  try {
    const page=await browser.newPage({viewport:{width:1440,height:1000}}), errors=[];
    page.on('pageerror',e=>errors.push(e.message));
    await page.goto(base+'#browser');
    await page.waitForFunction(()=>document.getElementById('browser-c')?._fullData?.length>0,
      null,{timeout:120000});
    const tab=page.locator('a[data-value="Predicted Translons"]');
    assert.equal(await tab.count(),1);
    assert.equal((await tab.innerText()).trim(),'Translons');
    assert.equal(await tab.evaluate(el=>el.closest('.dropdown-menu')),null);
    assert.equal(await tab.evaluate(el=>el.parentElement.previousElementSibling
      .querySelector('a').dataset.value),'Observatory');
    await tab.click();
    await page.locator('#predicted_translons-go').waitFor({state:'visible'});
    await page.waitForFunction(()=>document.getElementById('predicted_translons-dff')?.selectize);
    await page.evaluate(()=>document.getElementById('predicted_translons-dff').selectize.setValue('human_all_merged_l50'));
    await page.waitForFunction(()=>Shiny.shinyapp.$inputValues['predicted_translons-dff']==='human_all_merged_l50');
    await page.locator('#predicted_translons-go').click();
    await page.waitForFunction(()=>{
      const t=document.querySelector('#predicted_translons-translon_table .dataTables_scrollBody table');
      return t && jQuery.fn.dataTable.isDataTable(t) && jQuery(t).DataTable().page.info().recordsTotal>0;
    },null,{timeout:120000});
    const rows=await page.evaluate(()=>jQuery('#predicted_translons-translon_table .dataTables_scrollBody table').DataTable().page.info().recordsTotal);
    await page.locator('a[data-value="Observatory"]').click();
    await page.waitForFunction(()=>document.getElementById('observatory-selector-libraries_umap_plot')?._fullData?.length>0 &&
      !document.documentElement.classList.contains('shiny-busy'),null,{timeout:120000});
    await tab.click();
    await page.locator('#predicted_translons-translon_table .dataTables_scrollBody table').waitFor({state:'visible'});
    assert.equal(await page.evaluate(()=>jQuery('#predicted_translons-translon_table .dataTables_scrollBody table').DataTable().page.info().recordsTotal),rows);
    await page.screenshot({path:'/tmp/translons-navigation.png'});
    await page.goto(base+'#Predicted%20Translons');
    await page.locator('#predicted_translons-go').waitFor({state:'visible',timeout:60000});
    assert.deepEqual(errors,[]);
    assert.deepEqual(await page.locator('.shiny-output-error').allTextContents(),[]);
    console.log(JSON.stringify({topLevel:true,afterObservatory:true,rows,stateRetained:true,directUrl:true,errors}));
  } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
