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
    await page.locator('#browser_allsamp-mb_tabs a[data-value="Heatmap"]').click();
    await page.evaluate(() => document.getElementById('browser_allsamp-visible_groups').selectize.clear());
    await page.locator('#browser_allsamp-collapsed_clusters').uncheck();
    await page.locator('#browser_allsamp-collapsed_translons').check();
    await page.waitForFunction(n => {
      const p = document.getElementById('browser_allsamp-myPlotlyPlot');
      return p?._fullData?.[0]?.z.length === n && p._fullData[0].z[0].length === 7 &&
        p._fullLayout.xaxis.range[0] === 0.5 && !document.documentElement.classList.contains('shiny-busy');
    }, libraries);
    assert.equal(await page.evaluate(() => document.getElementById('browser_allsamp-myPlotlyPlot')._fullData[0].text[0].length), 7);
    assert.match(await page.evaluate(() => document.getElementById('browser_allsamp-myPlotlyPlot')._fullData[0].text[0][6]), /clean_cds/);
    await page.waitForFunction(() => document.getElementById('browser_allsamp-mb_bottom_gene')?._fullLayout?.xaxis.range[1] === 7.5);
    await page.waitForFunction(() => {
      const p = document.getElementById('browser_allsamp-myPlotlyPlot');
      return ['mb_top_summary', 'mb_bottom_gene'].every(id => {
        const track = document.getElementById('browser_allsamp-' + id);
        return Math.abs(track._fullLayout.xaxis._length - p._fullLayout.xaxis._length) < 2;
      });
    });
    await page.locator('#browser_allsamp-collapsed_clusters').check();
    await page.waitForFunction(() => {
      const p = document.getElementById('browser_allsamp-myPlotlyPlot');
      return p?._fullData?.[0]?.z.length === 5 && p._fullData[0].z[0].length === 7 &&
        !document.documentElement.classList.contains('shiny-busy');
    });
    const relative = await page.evaluate(() => {
      const trace = document.getElementById('browser_allsamp-myPlotlyPlot')._fullData[0];
      const baseline = document.getElementById('browser_allsamp-mb_top_summary')._fullData[0].y;
      return {zmin:trace.zmin,zmax:trace.zmax,correct:trace.z.every((row,i)=>row.every((z,j)=>{
        if (!(baseline[j] > 0)) return z === null;
        const expected = Math.max(-1,Math.min(1,Math.log2(trace.customdata[i][j]/baseline[j])));
        return Math.abs(z-expected)<1e-8;
      }))};
    });
    assert.deepEqual(relative, {zmin:-1,zmax:1,correct:true});
    await page.waitForFunction(() => {
      const p = document.getElementById('browser_allsamp-myPlotlyPlot');
      return ['mb_top_summary', 'mb_bottom_gene'].every(id =>
        Math.abs(document.getElementById('browser_allsamp-' + id)._fullLayout.xaxis._length - p._fullLayout.xaxis._length) < 2);
    });
    await page.evaluate(async () => {
      const p = document.getElementById('browser_allsamp-myPlotlyPlot');
      await Plotly.relayout(p, {'xaxis.range': [1, 3], 'yaxis.range': [1, 2]});
      p.emit('plotly_doubleclick');
    });
    await page.waitForFunction(() => {
      const p = document.getElementById('browser_allsamp-myPlotlyPlot');
      const side = document.getElementById('browser_allsamp-d');
      return p._fullLayout.xaxis.range[0] === 0.5 && p._fullLayout.xaxis.range[1] === 7.5 &&
        p._fullLayout.yaxis.range[1] === 5.5 && side?._fullLayout.yaxis.range[0] === 5.5;
    });
    await page.screenshot({path: '/tmp/megabrowser-collapsed-translons-desktop.png', fullPage: true});
    await page.setViewportSize({width: 390, height: 844});
    await page.waitForTimeout(1000);
    assert.equal(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth + 1), true);
    assert.equal(await page.locator('.mega-region-legend').isVisible(), true);
    assert.equal(await page.evaluate(() => [...document.querySelectorAll('.mega-region-legend span')].every(e =>
      e.getBoundingClientRect().right <= innerWidth)), true);
    await page.screenshot({path: '/tmp/megabrowser-collapsed-translons-mobile.png', fullPage: true});
    await page.setViewportSize({width: 1440, height: 1100});
    await page.locator('#browser_allsamp-collapsed_translons').uncheck();
    await page.waitForFunction(() => document.getElementById('browser_allsamp-myPlotlyPlot')?._fullData?.[0]?.z[0].length > 7 &&
      !document.documentElement.classList.contains('shiny-busy'));
    await page.evaluate(() => document.getElementById('browser_allsamp-d').emit('plotly_clickannotation', {fullAnnotation: {text: '1'}}));
    await page.waitForFunction(() => Shiny.shinyapp.$inputValues['browser_allsamp-mb_tabs'] === 'Factor enrichment' &&
      Shiny.shinyapp.$inputValues['browser_allsamp-factor_tabs'] === 'Result table');
    await page.locator('a[data-value="Translon enrichment"]').click();
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
