const {chromium} = require('playwright');
const fs = require('node:fs'), assert = require('node:assert/strict');
const base = process.argv[2] || 'http://127.0.0.1:7823/';
const prefix = 'browser_allsamp-', heatmap = prefix + 'myPlotlyPlot';

function instrument() {
  window.mbProfile = {calls: [], tasks: [], messages: []};
  let plotly;
  Object.defineProperty(window, 'Plotly', {configurable: true, get: () => plotly, set: value => {
    plotly = value;
    for (const method of ['newPlot', 'react', 'relayout', 'restyle', 'redraw']) {
      const original = value[method];
      if (typeof original !== 'function') continue;
      value[method] = function(...args) {
        const element = typeof args[0] === 'string' ? document.getElementById(args[0]) : args[0];
        if (!element?.id?.startsWith('browser_allsamp-')) return original.apply(this, args);
        const call = {method, id: element.id, start: performance.now()};
        window.mbProfile.calls.push(call);
        const result = original.apply(this, args);
        if (result?.then) result.then(() => {call.end = performance.now();});
        else call.end = performance.now();
        return result;
      };
    }
  }});
  new PerformanceObserver(list => window.mbProfile.tasks.push(...list.getEntries().map(e => ({start: e.startTime, duration: e.duration})))).observe({type: 'longtask', buffered: true});
}

async function wait(page, rows) {
  await page.waitForFunction(({id, rows}) => {
    const el = document.getElementById(id);
    return el?._fullData?.[0]?.z?.length === rows && !document.documentElement.classList.contains('shiny-busy');
  }, {id: heatmap, rows}, {timeout: 180000});
  await page.waitForTimeout(150);
}

async function measure(page, name, action, rows) {
  await page.evaluate(() => {window.mbProfile.calls = []; window.mbProfile.tasks = []; window.mbProfile.messages = [];});
  const start = performance.now();
  await action(); await wait(page, rows);
  const elapsed = performance.now() - start;
  const profile = await page.evaluate(id => {
    const el = document.getElementById(id), side = document.getElementById('browser_allsamp-d');
    return {...window.mbProfile, rows: el._fullData[0].z.length, positions: el._fullData[0].z[0].length,
      yRange: el._fullLayout.yaxis.range, sidebarRows: side._fullData.find(t => t.y?.length)?.y.length};
  }, heatmap);
  return {name, elapsed_ms: elapsed, ...profile};
}

(async () => {
  const browser = await chromium.launch({executablePath: '/opt/google/chrome/chrome', headless: true, args: ['--no-sandbox']});
  const results = [], errors = [];
  try {
    const page = await browser.newPage({viewport: {width: 1440, height: 1000}});
    await page.addInitScript(instrument);
    page.on('pageerror', error => errors.push(error.message));
    page.on('websocket', socket => socket.on('framereceived', event => {
      const bytes = Buffer.byteLength(event.payload);
      let ids = [];
      try {ids = Object.keys(JSON.parse(event.payload).values || {});} catch {}
      page.evaluate(({bytes, ids}) => window.mbProfile?.messages.push({bytes, ids, time: performance.now()}), {bytes, ids}).catch(() => {});
    }));
    await page.goto(base + '#MegaBrowser');
    await page.waitForFunction(() => document.getElementById('browser_allsamp-tx')?.value, null, {timeout: 120000});
    await page.evaluate(() => document.getElementById('browser_allsamp-gene').selectize.setValue('AMD1-ENSG00000123505'));
    await page.waitForFunction(() => Shiny.shinyapp.$inputValues['browser_allsamp-gene'] === 'AMD1-ENSG00000123505' &&
      Shiny.shinyapp.$inputValues['browser_allsamp-tx'] === 'ENST00000451850' &&
      !document.documentElement.classList.contains('shiny-busy'));
    await page.evaluate(() => document.getElementById('browser_allsamp-metadata').selectize.setValue(['TISSUE', 'CELL_LINE']));
    const settings = await page.evaluate(() => Object.fromEntries(['gene', 'tx', 'region_type', 'kmer', 'normalization',
      'min_count', 'clusters', 'heatmap_color', 'color_mult', 'metadata', 'enrichment_term'].map(field =>
      [field, Shiny.shinyapp.$inputValues['browser_allsamp-' + field]])));
    await page.locator('#' + prefix + 'go').click();
    await page.waitForFunction(id => document.getElementById(id)?._fullData?.[0]?.z?.length > 5, heatmap, {timeout: 180000});
    const rows = await page.evaluate(id => document.getElementById(id)._fullData[0].z.length, heatmap);
    await wait(page, rows);
    async function collapsed(value) {
      await page.evaluate(({prefix, value}) => {
        const checkbox = document.getElementById(prefix + 'collapsed_clusters');
        checkbox.checked = value; jQuery(checkbox).trigger('change');
      }, {prefix, value});
    }
    for (let trial = 1; trial <= 3; trial++) {
      results.push(await measure(page, 'collapsed-' + trial, () => collapsed(true), 5));
      results.push(await measure(page, 'full-' + trial, () => collapsed(false), rows));
      console.log(JSON.stringify(results.slice(-2).map(({name, elapsed_ms, rows}) => ({name, elapsed_ms, rows}))));
    }
    await collapsed(true); await wait(page, 5);
    await page.screenshot({path: '/tmp/megabrowser-collapsed.png'});
    await page.locator('a[data-value="Statistics"]').click();
    await page.waitForFunction(() => jQuery.fn.dataTable?.isDataTable('#browser_allsamp-stats table'));
    const initialStatistics = await page.evaluate(() => JSON.stringify(jQuery('#browser_allsamp-stats table').DataTable().rows().data().toArray()));
    await page.evaluate(() => {window.mbProfile.calls = []; window.mbProfile.messages = [];});
    const statsStart = performance.now();
    await page.evaluate(() => document.getElementById('browser_allsamp-enrichment_metadata').selectize.setValue('CELL_LINE'));
    await page.waitForFunction(before => jQuery.fn.dataTable?.isDataTable('#browser_allsamp-stats table') &&
      jQuery('#browser_allsamp-stats table').DataTable().rows().data().length > 0 &&
      JSON.stringify(jQuery('#browser_allsamp-stats table').DataTable().rows().data().toArray()) !== before &&
      !document.documentElement.classList.contains('shiny-busy'), initialStatistics, {timeout: 120000});
    const stats_ms = performance.now() - statsStart;
    const stats = await page.evaluate(() => ({calls: window.mbProfile.calls, messages: window.mbProfile.messages}));
    assert.ok(!stats.calls.some(c => c.id === 'browser_allsamp-myPlotlyPlot'), 'Enrichment must not redraw coverage');
    await page.screenshot({path: '/tmp/megabrowser-statistics.png'});
    await page.locator('a[data-value="Result table"]').click();
    await page.waitForFunction(() => jQuery.fn.dataTable?.isDataTable('#browser_allsamp-result_table table') &&
      jQuery('#browser_allsamp-result_table table').DataTable().page.info().recordsTotal > 0);
    const samples = await page.evaluate(() => jQuery('#browser_allsamp-result_table table').DataTable().page.info().recordsTotal);
    assert.equal(samples, rows);
    assert.deepEqual(errors, []);
    assert.deepEqual(await page.locator('.shiny-output-error').allTextContents(), []);
    fs.writeFileSync('/tmp/megabrowser-browser.json', JSON.stringify({results, stats_ms, settings,
      stats, samples, errors}, null, 2));
    console.log(JSON.stringify({results: results.map(({name, elapsed_ms, rows, positions}) => ({name, elapsed_ms, rows, positions})), samples, errors}));
  } finally {await browser.close();}
})().catch(error => {console.error(error); process.exitCode = 1;});
