// Run against a locally loaded app with the human collection and real metadata.
// NODE_PATH=/path/to/node_modules node benchmarks/observatory-url-smoke.cjs [URL]
const { chromium } = require('playwright');
const zlib = require('node:zlib');
const assert = require('node:assert/strict');
const base = process.argv[2] || 'http://127.0.0.1:7820/';
const selector = 'observatory-selector-';
const browse = 'observatory-browser_obs-';
const settings = {gene:'AMD1-ENSG00000123505', kmer:9, frames_type:'lines',
  log_scale:true, colors:'Color_blind', frames_subset:['red','blue'],
  summary_track:false, go:true};
const state = {v:1, exp:'all_samples-Homo_sapiens', color_by:['tissue','cell_line'],
  view:'browser', browser:settings,
  selections:{active:'2', order:['1','2'], labels:{'1':'First group','2':'Saved group'},
    runs:['DRR277023'], p:{'1':[1],'2':[1]}, d:{'1':[1],'2':[1]}}};
const url = value => base + '#Observatory?obs_state=' +
  zlib.deflateSync(JSON.stringify(value)).toString('base64url');

async function ready(page, plot = true) {
  await page.waitForFunction(({selector,browse,plot}) => {
    const label = document.getElementById(selector+'library_selection-active_selection_label');
    return label?.value === 'Saved group' &&
      (!plot || document.getElementById(browse+'browser_plot')?.data?.length > 0);
  }, {selector,browse,plot}, {timeout:120000});
}

async function snapshot(page) {
  return page.evaluate(({selector,browse}) => {
    const table = document.querySelector('#'+selector+'libraries_data_table .dataTables_scrollBody table');
    const plot = document.getElementById(selector+'libraries_umap_plot');
    return {info:table ? jQuery(table).DataTable().page.info():null,
      label:document.getElementById(selector+'library_selection-active_selection_label')?.value,
      active:document.getElementById(selector+'library_selection-active_selection_id')?.value,
      highlighted:plot?.data?.reduce((n,t)=>n+(t.selectedpoints?.length || 0),0),
      traces:document.getElementById(browse+'browser_plot')?.data?.length || 0,
      settings:Object.fromEntries(Object.entries(Shiny.shinyapp.$inputValues)
        .filter(([key])=>key.startsWith(browse)).map(([key,value])=>[key.slice(browse.length).split(':')[0],value])),
      errors:[...document.querySelectorAll('.shiny-output-error')].map(e=>e.textContent)};
  }, {selector,browse});
}

async function checkSelection(page) {
  await page.locator('a[data-value="Select libraries"]').click();
  await page.waitForFunction(selector => {
    const table = document.querySelector('#'+selector+'libraries_data_table .dataTables_scrollBody table');
    return table && jQuery(table).DataTable().page.info()?.recordsDisplay === 1 &&
      document.getElementById(selector+'libraries_umap_plot')?.data?.some(t=>t.selectedpoints?.length);
  }, selector, {timeout:60000});
  const result = await snapshot(page);
  assert.equal(result.active,'2'); assert.equal(result.label,'Saved group');
  assert.equal(result.info.recordsTotal,1); assert.equal(result.highlighted,1);
  assert.deepEqual(result.errors,[]);
  return result;
}

(async()=>{
  const browser = await chromium.launch({executablePath:process.env.CHROME || '/opt/google/chrome/chrome',
    headless:true,args:['--no-sandbox']});
  try {
    const context = await browser.newContext({viewport:{width:1600,height:1000}});
    const page = await context.newPage();
    const errors = [];
    page.on('pageerror',e=>errors.push(e.message));
    let copied;
    page.on('websocket',socket=>socket.on('framereceived',event=>{
      try { const msg = JSON.parse(event.payload); copied = msg.custom?.['ribocrypt-copy-url']?.text || copied; }
      catch (_) { /* Ignore protocol frames that are not JSON. */ }
    }));
    await page.goto(url(state)); await ready(page);
    const first = await checkSelection(page);
    for (const [key,value] of Object.entries(settings)) {
      if(key !== 'go') assert.deepEqual(first.settings[key],value,key);
    }
    await page.locator('a[data-value="Browse"]').click();
    await page.locator('#'+browse+'toggle_settings').click();
    await page.locator('#'+browse+'floating_settings a[data-value="Settings"]').click();
    await page.locator('#'+browse+'clip_button').click();
    for(let i=0;!copied && i<100;i++) await page.waitForTimeout(100);
    assert.ok(copied,'Get URL must produce a URL');
    const shared = JSON.parse(zlib.inflateSync(Buffer.from(copied.split('obs_state=')[1],'base64url')));
    assert.deepEqual(shared.selections.labels,state.selections.labels);
    for(const [key,value] of Object.entries(settings)) assert.deepEqual(shared.browser[key],value,key);
    const second = await context.newPage();
    await second.goto(copied); await ready(second);
    const restored = await checkSelection(second);
    for(const [key,value] of Object.entries(shared.browser)) {
      if(key !== 'go') assert.deepEqual(restored.settings[key],value,key);
    }
    await second.getByRole('button',{name:'Remove subset',exact:true}).click();
    await second.waitForFunction(selector => jQuery('#'+selector+
      'libraries_data_table .dataTables_scrollBody table').DataTable().page.info()?.recordsTotal > 1000,
      selector,{timeout:30000});
    await second.getByRole('button',{name:'Subset to page',exact:true}).click();
    await second.waitForFunction(selector => jQuery('#'+selector+
      'libraries_data_table .dataTables_scrollBody table').DataTable().page.info()?.recordsTotal === 15,
      selector,{timeout:30000});
    await second.waitForTimeout(1500);
    assert.equal((await snapshot(second)).info.recordsDisplay,15);
    await second.getByRole('button',{name:'Remove subset',exact:true}).click();
    await second.waitForFunction(selector => jQuery('#'+selector+
      'libraries_data_table .dataTables_scrollBody table').DataTable().page.info()?.recordsTotal > 1000,
      selector,{timeout:30000});
    await second.locator('#'+selector+'libraries_data_table .dataTables_scrollBody tbody tr').first().click();
    await second.waitForTimeout(500);
    assert.ok((await snapshot(second)).info.recordsTotal > 1000,'Clicking a row must not subset');
    await second.getByRole('button',{name:'Subset to selected',exact:true}).click();
    await second.waitForFunction(selector => jQuery('#'+selector+
      'libraries_data_table .dataTables_scrollBody table').DataTable().page.info()?.recordsTotal === 1,
      selector,{timeout:30000});
    const stopped = structuredClone(state); stopped.browser.go = false;
    const third = await context.newPage(); await third.goto(url(stopped)); await ready(third,false);
    await third.waitForTimeout(3000); assert.equal((await snapshot(third)).traces,0);
    assert.deepEqual(errors,[]);
    await third.close();
    await page.close();
    await second.bringToFront();
    await second.screenshot({path:'/tmp/observatory-url-smoke.png',fullPage:true,timeout:60000});
    console.log(JSON.stringify({restoredLibraries:restored.info.recordsTotal,
      highlighted:restored.highlighted,browserTraces:restored.traces,
      settingsRoundtrip:true,resetAndPageSubset:true,selectedSubset:true,goFalse:true,errors}));
  } finally { await browser.close(); }
})().catch(error=>{console.error(error);process.exitCode=1;});
