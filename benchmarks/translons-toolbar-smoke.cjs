const {chromium}=require('playwright');
const assert=require('node:assert/strict'),fs=require('node:fs');
const base=process.argv[2] || 'http://127.0.0.1:7821/';
const selector='#predicted_translons-translon_table .dataTables_scrollBody table';
async function snapshot(page) {
  return page.evaluate(s=>{
    const t=jQuery(s).DataTable();
    return {columns:t.columns().header().toArray().map(x=>x.textContent.trim()),
      visible:t.columns().visible().toArray(),info:t.page.info()};
  },selector);
}
async function download(page,button,path) {
  const [file]=await Promise.all([page.waitForEvent('download',{timeout:120000}),button.click()]);
  await file.saveAs(path);assert.equal(await file.failure(),null);
  return fs.readFileSync(path,'utf8').split('\n')[0].toLowerCase();
}
(async()=>{
  const browser=await chromium.launch({executablePath:'/opt/google/chrome/chrome',headless:true,args:['--no-sandbox']});
  try {
    const page=await browser.newPage({viewport:{width:1440,height:1000}}),errors=[];
    page.on('pageerror',e=>errors.push(e.message));
    await page.goto(base+'#Predicted%20Translons');
    await page.waitForFunction(()=>Shiny.shinyapp.$inputValues['predicted_translons-dff']==='human_all_merged_l50',null,{timeout:120000});
    assert.equal(await page.locator('#predicted_translons-simplified').isChecked(),false);
    await page.locator('#predicted_translons-go').click();
    await page.waitForFunction(s=>{
      const t=document.querySelector(s);return t&&jQuery.fn.dataTable.isDataTable(t)&&jQuery(t).DataTable().page.info().recordsTotal>0;
    },selector,{timeout:120000});
    const full=await snapshot(page);assert.ok(full.visible.every(Boolean));
    const groups=await page.evaluate(()=>{
      const toolbar=document.querySelector('.rc-translon-toolbar');
      const study=toolbar.querySelector('.rc-translon-study-group');
      return {count:toolbar.children.length,searchWithStudy:!!study.querySelector('#predicted_translons-go'),
        color:getComputedStyle(document.getElementById('predicted_translons-trigger_download_excel')).backgroundColor,
        gap:parseFloat(getComputedStyle(toolbar).columnGap)};
    });
    assert.equal(groups.count,3);assert.equal(groups.searchWithStudy,true);
    assert.equal(groups.color,'rgb(33, 115, 70)');assert.ok(groups.gap>=32);
    await page.evaluate(s=>jQuery(s).DataTable().page(1).draw('page'),selector);
    await page.waitForFunction(s=>jQuery(s).DataTable().page.info().page===1,selector);
    await page.locator('#predicted_translons-simplified').check();
    const simplified=await snapshot(page);
    const expected=full.columns.filter((name,i)=>i<3 || ['external_gene_name','type','length','mapability_score','biotype'].includes(name.toLowerCase()) || name.startsWith('LOC_'));
    assert.deepEqual(simplified.columns.filter((_,i)=>simplified.visible[i]),expected);
    assert.equal(simplified.info.page,1);assert.equal(simplified.info.recordsTotal,full.info.recordsTotal);
    const header=await download(page,page.getByText('Download current page (CSV)',{exact:true}),'/tmp/translons-current.csv');
    assert.ok(!header.includes('alignment'));assert.ok(header.includes('mapability_score'));
    const allHeader=await download(page,page.locator('#predicted_translons-trigger_download_csv'),'/tmp/translons-full.csv');
    assert.ok(allHeader.includes('alignment'));assert.ok(allHeader.includes('sequence_aa'));
    await page.locator('#predicted_translons-simplified').uncheck();
    const restored=await snapshot(page);assert.ok(restored.visible.every(Boolean));assert.equal(restored.info.page,1);
    await page.screenshot({path:'/tmp/translons-toolbar-desktop.png'});
    await page.locator('#predicted_translons-simplified').check();
    // Rebuilding the table while checked must retain simplified visibility.
    await page.locator('#predicted_translons-go').click();
    await page.waitForTimeout(2000);
    const rebuilt=await snapshot(page);
    assert.deepEqual(rebuilt.columns.filter((_,i)=>rebuilt.visible[i]),expected);
    await page.setViewportSize({width:390,height:844});await page.waitForTimeout(500);
    const boxes=await page.evaluate(()=>{
      const r=e=>{const b=e.getBoundingClientRect();return {top:b.top,bottom:b.bottom,left:b.left,right:b.right};};
      return {toolbar:r(document.querySelector('.rc-translon-toolbar')),table:r(document.querySelector('#predicted_translons-translon_table')),
        controls:['dff','go','simplified','trigger_download_csv','trigger_download_excel'].map(id=>r(document.getElementById('predicted_translons-'+id))),width:innerWidth};
    });
    assert.ok(boxes.toolbar.bottom<=boxes.table.top);
    assert.ok(boxes.controls.every(b=>b.left>=0 && b.right<=boxes.width));
    await page.screenshot({path:'/tmp/translons-toolbar-mobile.png'});
    assert.deepEqual(errors,[]);assert.deepEqual(await page.locator('.shiny-output-error').allTextContents(),[]);
    console.log(JSON.stringify({rows:full.info.recordsTotal,fullColumns:full.columns.length,simplifiedColumns:expected,retainedPage:true,currentCsv:true,fullCsv:true,mobile:true,errors}));
  } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
