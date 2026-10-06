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
    await page.waitForFunction(() => Shiny.shinyapp.$inputValues['browser_allsamp-gene'] === 'AMD1-ENSG00000123505' &&
      Shiny.shinyapp.$inputValues['browser_allsamp-tx'] === 'ENST00000451850' &&
      !document.documentElement.classList.contains('shiny-busy'));
    await page.evaluate(() => {
      const select = document.getElementById('browser_allsamp-gene').selectize;
      select.addOption({value: 'ATF4-ENSG00000128272', label: 'ATF4-ENSG00000128272'});
      select.setValue('ATF4-ENSG00000128272');
    });
    await page.waitForFunction(() => Shiny.shinyapp.$inputValues['browser_allsamp-gene'] === 'ATF4-ENSG00000128272' &&
      document.getElementById('browser_allsamp-tx').value && !document.documentElement.classList.contains('shiny-busy'));
    await page.evaluate(() => {
      document.getElementById('browser_allsamp-tx').selectize.setValue('ENST00000674920');
      for (const id of ['add_translon', 'add_translons_transcode']) {
        const input = document.getElementById('browser_allsamp-' + id);
        input.checked = true; jQuery(input).trigger('change');
      }
    });
    await page.waitForFunction(() => Shiny.shinyapp.$inputValues['browser_allsamp-tx'] === 'ENST00000674920' &&
      Shiny.shinyapp.$inputValues['browser_allsamp-add_translons_transcode']);
    await page.locator('#browser_allsamp-go').click();
    await page.waitForFunction(() => document.getElementById('browser_allsamp-myPlotlyPlot')?._fullData?.[0]?.z.length > 5 &&
      !document.documentElement.classList.contains('shiny-busy'));
    const libraries = await page.evaluate(() => document.getElementById('browser_allsamp-myPlotlyPlot')._fullData[0].z.length);
    assert.deepEqual(await page.locator('#browser_allsamp-mb_tabs > li > a').allTextContents(),
      ['Heatmap', 'Factor enrichment', 'Translon enrichment']);
    await page.locator('a[data-value="Factor enrichment"]').click();
    assert.deepEqual(await page.locator('#browser_allsamp-factor_tabs > li > a').allTextContents(),
      ['Enrichment', 'Group summary', 'Statistics', 'Result table']);
    assert.equal(await page.evaluate(() => !!document.getElementById('browser_allsamp-translon_ratios_plot')?._fullLayout), false);
    const start = performance.now();
    await page.locator('a[data-value="Translon enrichment"]').click();
    await page.waitForFunction(() => document.getElementById('browser_allsamp-translon_ratios_plot')?._fullData?.[0]?.type === 'box');
    const elapsed = performance.now() - start;
    await page.locator('a[data-value="Regions"]').click();
    await page.waitForFunction(() => document.querySelectorAll('#browser_allsamp-translon_regions tbody tr').length === 6);
    const regionText = await page.locator('#browser_allsamp-translon_regions').innerText();
    assert.match(regionText, /T1, TC1/);
    assert.match(regionText, /clean_cds/);
    assert.match(regionText, /973/);
    await page.locator('a[data-value="Cluster statistics"]').click();
    await page.waitForFunction(() => document.querySelectorAll('#browser_allsamp-translon_statistics tbody tr').length === 5);
    const before = await page.locator('#browser_allsamp-translon_statistics tbody').innerText();
    await page.locator('#browser_allsamp-collapsed_clusters').check();
    await page.waitForTimeout(300);
    assert.equal(await page.locator('#browser_allsamp-translon_statistics tbody').innerText(), before);
    await page.locator('#browser_allsamp-visible_groups').locator('..').locator('.selectize-input').click();
    await page.waitForFunction(() => Object.keys(document.getElementById('browser_allsamp-visible_groups').selectize.options).length === 5);
    await page.keyboard.press('Escape');
    await page.evaluate(() => document.getElementById('browser_allsamp-visible_groups').selectize.setValue(['1']));
    await page.waitForTimeout(300);
    assert.equal(await page.locator('#browser_allsamp-translon_statistics tbody').innerText(), before);
    const download = page.waitForEvent('download');
    await page.locator('#browser_allsamp-translon_ratios_csv').click();
    const file = await download;
    await file.saveAs('/tmp/megabrowser-translon-ratios.csv');
    const rows = fs.readFileSync('/tmp/megabrowser-translon-ratios.csv', 'utf8').trim().split('\n');
    assert.equal(rows.length - 1, libraries * 15);
    assert.match(rows[0], /Numerator_label.*Denominator_label/);
    assert.match(rows[0], /Numerator_density.*Denominator_density.*Ratio/);
    await page.locator('a[data-value="Regions"]').click();
    await page.locator('#browser_allsamp-translon_add_region').click();
    await page.waitForSelector('#browser_allsamp-translon_region_label');
    assert.equal(await page.locator('#browser_allsamp-translon_region_label').inputValue(), 'UR1');
    await page.locator('#browser_allsamp-translon_region_coordinates').fill('40:39');
    await page.locator('#browser_allsamp-translon_region_save').click();
    await page.waitForFunction(() => document.getElementById('browser_allsamp-translon_region_error')?.textContent.includes('Start must not exceed'));
    await page.locator('#browser_allsamp-translon_region_label').fill('Custom ATF4');
    await page.locator('#browser_allsamp-translon_region_coordinates').fill('40:80');
    await page.locator('#browser_allsamp-translon_region_save').click();
    await page.waitForFunction(() => document.querySelectorAll('#browser_allsamp-translon_regions tbody tr').length === 7);
    assert.match(await page.locator('#browser_allsamp-translon_regions').innerText(), /Custom ATF4/);
    assert.match(await page.locator('#browser_allsamp-translon_regions').innerText(), /40-80/);
    await page.locator('#browser_allsamp-translon_add_region').click();
    await page.waitForSelector('#browser_allsamp-translon_region_label');
    assert.equal(await page.locator('#browser_allsamp-translon_region_label').inputValue(), 'UR2');
    await page.locator('#browser_allsamp-translon_region_coordinates').fill('40');
    await page.locator('#browser_allsamp-translon_region_save').click();
    await page.waitForFunction(() => document.querySelectorAll('#browser_allsamp-translon_regions tbody tr').length === 8);
    assert.match(await page.locator('#browser_allsamp-translon_regions').innerText(), /40-40/);
    const customDownload = page.waitForEvent('download');
    await page.locator('#browser_allsamp-translon_ratios_csv').click();
    await (await customDownload).saveAs('/tmp/megabrowser-custom-region-ratios.csv');
    const customRows = fs.readFileSync('/tmp/megabrowser-custom-region-ratios.csv', 'utf8').trim().split('\n');
    assert.equal(customRows.length - 1, libraries * 28);
    assert.match(customRows.join('\n'), /Custom ATF4/);
    await page.locator('a[data-value="Ratios"]').click();
    await page.screenshot({path: '/tmp/megabrowser-translons-desktop.png', fullPage: true});
    await page.setViewportSize({width: 390, height: 844});
    await page.waitForTimeout(750);
    await page.screenshot({path: '/tmp/megabrowser-translons-mobile.png', fullPage: true});
    assert.deepEqual(errors, []);
    assert.deepEqual(await page.locator('.shiny-output-error').allTextContents(), []);
    console.log({atf4_transcript: 'ENST00000674920', regions: 6, csv_ratios: rows.length - 1, first_tab_ms: elapsed, errors});
  } finally {await browser.close();}
})().catch(error => {console.error(error); process.exitCode = 1;});
