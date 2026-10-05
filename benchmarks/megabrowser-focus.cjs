const {chromium} = require('playwright');
const assert = require('node:assert/strict');
const fs = require('node:fs');

(async () => {
  const browser = await chromium.launch({executablePath: '/opt/google/chrome/chrome', args: ['--no-sandbox']});
  try {
    const page = await browser.newPage({viewport: {width: 1440, height: 1100}});
    const errors = [];
    page.on('pageerror', error => errors.push(error.message));
    await page.goto((process.argv[2] || 'http://127.0.0.1:7823/') + '#MegaBrowser');
    await page.waitForFunction(() => document.getElementById('browser_allsamp-tx')?.value);
    await page.evaluate(() => document.getElementById('browser_allsamp-gene').selectize.setValue('AMD1-ENSG00000123505'));
    await page.waitForFunction(() => Shiny.shinyapp.$inputValues['browser_allsamp-gene'] === 'AMD1-ENSG00000123505' &&
      Shiny.shinyapp.$inputValues['browser_allsamp-tx'] === 'ENST00000451850' &&
      !document.documentElement.classList.contains('shiny-busy'));
    await page.locator('#browser_allsamp-go').click();
    await page.waitForFunction(() => document.getElementById('browser_allsamp-myPlotlyPlot')?._fullData?.[0]?.z.length > 5);
    await page.locator('#browser_allsamp-visible_groups').locator('..').locator('.selectize-input').click();
    await page.waitForFunction(() => Object.keys(document.getElementById('browser_allsamp-visible_groups').selectize.options).length === 5);
    await page.keyboard.press('Escape');
    const full = await page.evaluate(() => ({
      rows: document.getElementById('browser_allsamp-myPlotlyPlot')._fullData[0].z.length,
      dtick: document.getElementById('browser_allsamp-myPlotlyPlot')._fullLayout.yaxis.dtick,
      groups: Object.keys(document.getElementById('browser_allsamp-visible_groups').selectize.options)
    }));
    if (full.rows > 100) assert.ok(full.dtick > 1, 'Large views must retain automatic tick spacing');
    const group = full.groups[1];
    const expected = await page.evaluate(group => {
      const label = document.getElementById('browser_allsamp-visible_groups').selectize.options[group].label;
      return Number(label.match(/\(([\d,]+) libraries\)/)[1].replaceAll(',', ''));
    }, group);
    const start = performance.now();
    await page.evaluate(group => document.getElementById('browser_allsamp-visible_groups').selectize.setValue([group]), group);
    await page.waitForFunction(rows => document.getElementById('browser_allsamp-myPlotlyPlot')._fullData[0].z.length === rows &&
      !document.documentElement.classList.contains('shiny-busy'), expected);
    const focused_ms = performance.now() - start;
    await page.locator('#browser_allsamp-collapsed_clusters').check();
    await page.waitForFunction(() => document.getElementById('browser_allsamp-myPlotlyPlot')._fullData[0].z.length === 1);
    const annotation = await page.evaluate(() => document.getElementById('browser_allsamp-d')._fullLayout.annotations[0].text);
    assert.equal(annotation, group);
    await page.screenshot({path: '/tmp/megabrowser-focused.png'});
    const second = full.groups[3];
    await page.evaluate(groups => document.getElementById('browser_allsamp-visible_groups').selectize.setValue(groups), [second, group]);
    await page.waitForFunction(() => document.getElementById('browser_allsamp-myPlotlyPlot')._fullData[0].z.length === 2);
    await page.waitForFunction(groups => {
      const annotations = document.getElementById('browser_allsamp-d')._fullLayout.annotations;
      return annotations?.length === 2 && annotations.every(a => groups.includes(a.text));
    }, [group, second]);
    await page.locator('#browser_allsamp-mb_tabs a[data-value="Result table"]').click();
    await page.waitForFunction(() => jQuery.fn.dataTable?.isDataTable('#browser_allsamp-result_table table') &&
      jQuery('#browser_allsamp-result_table table').DataTable().page.info().recordsTotal > 0);
    const samples = await page.evaluate(() => jQuery('#browser_allsamp-result_table table').DataTable().page.info().recordsTotal);
    assert.equal(samples, full.rows);
    await page.locator('#browser_allsamp-mb_tabs a[data-value="Heatmap"]').click();
    await page.evaluate(() => document.getElementById('browser_allsamp-visible_groups').selectize.clear());
    await page.waitForFunction(() => document.getElementById('browser_allsamp-myPlotlyPlot')._fullData[0].z.length === 5);
    await page.locator('#browser_allsamp-collapsed_clusters').uncheck();
    await page.waitForFunction(rows => document.getElementById('browser_allsamp-myPlotlyPlot')._fullData[0].z.length === rows, full.rows);
    assert.deepEqual(errors, []);
    assert.deepEqual(await page.locator('.shiny-output-error').allTextContents(), []);
    const result = {group, focused_libraries: expected, full_libraries: samples, focused_ms, errors};
    fs.writeFileSync('/tmp/megabrowser-focus.json', JSON.stringify(result, null, 2));
    console.log(result);
  } finally {await browser.close();}
})().catch(error => {console.error(error); process.exitCode = 1;});
