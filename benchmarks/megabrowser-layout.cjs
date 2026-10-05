const {chromium} = require('playwright');
const assert = require('node:assert/strict');

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
    await page.locator('#browser_allsamp-collapsed_clusters').check();
    await page.locator('#browser_allsamp-go').click();
    await page.waitForFunction(() => document.getElementById('browser_allsamp-myPlotlyPlot')?._fullData?.[0]?.z.length === 5);
    const assertHeatmapAlignment = async () => {
      await page.waitForFunction(() => {
        const plot = document.getElementById('browser_allsamp-myPlotlyPlot');
        const sidebar = document.querySelector('.mega-plot-sidebar .js-plotly-plot');
        if (!plot?._fullLayout || !sidebar?._fullLayout) return false;
        const extent = graph => {
          const axis = graph._fullLayout.yaxis;
          const top = graph.getBoundingClientRect().top + axis._offset;
          return [top, top + axis._length];
        };
        return extent(plot).every((value, i) => Math.abs(value - extent(sidebar)[i]) < 1);
      });
      assert.deepEqual(await page.evaluate(() => {
        const axis = document.getElementById('browser_allsamp-myPlotlyPlot')._fullLayout.xaxis;
        return [axis.showticklabels, axis.ticks, axis.automargin];
      }), [false, '', false]);
    };
    await assertHeatmapAlignment();
    await page.locator('#browser_allsamp-myPlotlyPlot .nsewdrag').dblclick();
    await page.waitForTimeout(200);
    await page.waitForFunction(() => {
      const axis = document.getElementById('browser_allsamp-d')?._fullLayout?.yaxis;
      return axis?.range[0] === 5.5 && axis.range[1] === 0.5;
    });
    await page.evaluate(() => Plotly.relayout('browser_allsamp-myPlotlyPlot', {'yaxis.range': [1.25, 3.75]}));
    await page.waitForFunction(() => {
      const range = document.getElementById('browser_allsamp-d')._fullLayout.yaxis.range;
      return range[0] === 4.75 && range[1] === 2.25;
    });
    await page.locator('#browser_allsamp-myPlotlyPlot .nsewdrag').dblclick();
    await page.waitForTimeout(200);
    await page.waitForFunction(() => {
      const range = document.getElementById('browser_allsamp-d')._fullLayout.yaxis.range;
      return range[0] === 5.5 && range[1] === 0.5;
    });
    await page.locator('#browser_allsamp-collapsed_clusters').uncheck();
    await page.waitForFunction(() => document.getElementById('browser_allsamp-myPlotlyPlot')?._fullData?.[0]?.z.length > 5);
    await assertHeatmapAlignment();
    await page.locator('#browser_allsamp-collapsed_clusters').check();
    await page.waitForFunction(() => document.getElementById('browser_allsamp-myPlotlyPlot')?._fullData?.[0]?.z.length === 5);
    await page.locator('.mega-focus-control').check();
    await page.waitForFunction(() => document.getElementById('browser_allsamp-myPlotlyPlot')._fullLayout.height === 700);
    assert.equal(await page.locator('.mega-summary-track').isVisible(), false);
    assert.equal(await page.locator('.mega-gene-track').isVisible(), false);
    await page.locator('.mega-focus-control').uncheck();
    await page.evaluate(() => {
      const slider = document.querySelector('.mega-height-control');
      slider.value = 900; slider.dispatchEvent(new Event('input', {bubbles: true}));
    });
    await page.waitForFunction(() => document.querySelector('.mega-plot-grid').getBoundingClientRect().height === 900);
    await page.waitForFunction(() => document.getElementById('browser_allsamp-myPlotlyPlot')._fullLayout.height === 675);
    await page.locator('.mega-sidebar-control').uncheck();
    await page.waitForFunction(() => getComputedStyle(document.querySelector('.mega-plot-sidebar')).display === 'none');
    await page.evaluate(() => Plotly.relayout('browser_allsamp-myPlotlyPlot', {'xaxis.range': [100, 200], 'yaxis.range': [1.5, 3.5]}));
    await page.locator('.mega-reset-view').click();
    await page.waitForFunction(() => {
      const plot = document.getElementById('browser_allsamp-myPlotlyPlot');
      return plot._fullLayout.yaxis.range[0] === 0.5 && plot._fullLayout.yaxis.range[1] === 5.5;
    });
    await page.locator('.mega-sidebar-control').check();
    await page.screenshot({path: '/tmp/megabrowser-layout-desktop.png'});
    await page.locator('a[data-value="Enrichment"]').click();
    await page.waitForFunction(() => document.getElementById('browser_allsamp-e')?._fullData?.length);
    await page.evaluate(() => document.getElementById('browser_allsamp-enrichment_metadata').selectize.setValue('CELL_LINE'));
    await page.waitForFunction(() => !document.documentElement.classList.contains('shiny-busy'));
    await page.locator('#browser_allsamp-mb_tabs a[data-value="Heatmap"]').click();
    await page.setViewportSize({width: 390, height: 844});
    await page.waitForTimeout(750);
    await page.screenshot({path: '/tmp/megabrowser-layout-mobile.png', fullPage: true});
    const width = await page.evaluate(() => ({body: document.documentElement.scrollWidth, viewport: innerWidth}));
    assert.ok(width.body <= width.viewport + 1, JSON.stringify(width));
    await page.locator('#browser_allsamp-toggle_settings').click();
    await page.locator('#browser_allsamp-floating_settings').waitFor({state: 'visible'});
    const panel = await page.locator('#browser_allsamp-floating_settings').boundingBox();
    assert.ok(panel.x >= 0 && panel.x + panel.width <= 390, JSON.stringify(panel));
    await page.locator('#browser_allsamp-toggle_settings').click();
    await page.setViewportSize({width: 1440, height: 1100});
    await page.evaluate(() => {
      const radio = document.querySelector('input[name="browser_allsamp-plotType"][value="ggplot2"]');
      radio.checked = true; jQuery(radio).trigger('change');
    });
    await page.locator('#browser_allsamp-go').click();
    await page.waitForFunction(() => document.getElementById('browser_allsamp-mb_subplot_static')?._fullLayout?.images?.length > 0);
    assert.ok(await page.locator('.mega-reset-view').isDisabled());
    assert.ok(await page.locator('.mega-focus-control').isDisabled());
    await page.screenshot({path: '/tmp/megabrowser-layout-static.png'});
    assert.deepEqual(errors, []);
    assert.deepEqual(await page.locator('.shiny-output-error').allTextContents(), []);
    console.log('Height, sidebar, reset, enrichment, desktop/mobile, settings and static checks passed');
  } finally {await browser.close();}
})().catch(error => {console.error(error); process.exitCode = 1;});
