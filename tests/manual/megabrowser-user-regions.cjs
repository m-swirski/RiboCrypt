const { chromium } = require('playwright');
const assert = require('node:assert/strict');

(async () => {
  const browser = await chromium.launch({ executablePath: '/opt/google/chrome/chrome', args: ['--no-sandbox'] });
  try {
    const page = await browser.newPage({ viewport: { width: 1440, height: 1100 } });
    page.setDefaultTimeout(90000);
    const url = process.argv[2] || 'http://127.0.0.1:7826/#MegaBrowser';
    const deadline = Date.now() + 90000;
    while (true) {
      try { if ((await fetch(url)).ok) break; } catch (_) {}
      assert(Date.now() < deadline, 'App server did not start');
      await new Promise(resolve => setTimeout(resolve, 500));
    }
    await page.goto(url);
    await page.waitForFunction(() => document.getElementById('browser_allsamp-tx')?.value &&
      !document.documentElement.classList.contains('shiny-busy'));
    await page.locator('#browser_allsamp-go').click();
    await page.waitForFunction(() => document.getElementById('browser_allsamp-mb_bottom_gene')?._fullData &&
      !document.documentElement.classList.contains('shiny-busy'));
    const labels = () => page.evaluate(() => document.getElementById('browser_allsamp-mb_bottom_gene')._fullData.flatMap(t => t.text || []));
    const original = await labels();
    const heatmap = await page.evaluate(() => JSON.stringify(document.getElementById('browser_allsamp-myPlotlyPlot')._fullData[0].z));
    await page.locator('a[data-value="Translon enrichment"]').click();
    await page.locator('a[data-value="Regions"]').click();
    for (const [label, coordinates] of [['UR1', '40:80'], ['My region', '90']]) {
      await page.locator('#browser_allsamp-translon_add_region').click();
      await page.locator('#browser_allsamp-translon_region_label').fill(label);
      await page.locator('#browser_allsamp-translon_region_coordinates').fill(coordinates);
      await page.locator('#browser_allsamp-translon_region_save').click();
      await page.locator('#shiny-modal').waitFor({ state: 'hidden' });
    }
    await page.locator('#browser_allsamp-mb_tabs a[data-value="Heatmap"]').click();
    await page.waitForFunction(() => document.getElementById('browser_allsamp-mb_bottom_gene')._fullData
      .some(t => [].concat(t.text || []).some(label => String(label).includes('UR1'))));
    const updated = (await labels()).join(' ');
    assert.match(updated, /UR1/);
    assert.match(updated, /My region/);
    for (const label of original) assert(updated.includes(label));
    assert.equal(await page.evaluate(() => JSON.stringify(document.getElementById('browser_allsamp-myPlotlyPlot')._fullData[0].z)), heatmap);
    assert.deepEqual(await page.locator('.shiny-output-error:visible').allTextContents(), []);
    await page.screenshot({ path: '/tmp/ribocrypt-user-regions-heatmap.png', fullPage: true });
    console.log('PASS popup regions appear in Heatmap gene track; existing annotations and heatmap values preserved');
  } finally { await browser.close(); }
})().catch(error => { console.error(error); process.exitCode = 1; });
