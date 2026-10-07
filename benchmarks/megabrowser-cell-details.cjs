const {chromium} = require('playwright');
const assert = require('node:assert/strict');

(async () => {
  const browser = await chromium.launch({executablePath:'/opt/google/chrome/chrome', args:['--no-sandbox']});
  try {
    const page = await browser.newPage({viewport:{width:1440,height:1100}});
    page.setDefaultTimeout(60000);
    const errors = [];
    page.on('console', message => {
      if (message.type() === 'error' && /\[shiny\]|TypeError/.test(message.text())) errors.push(message.text());
    });
    page.on('pageerror', error => errors.push(error.message));
    await page.goto((process.argv[2] || 'http://127.0.0.1:7826/') + '#MegaBrowser');
    await page.waitForFunction(() => window.Shiny?.shinyapp?.$inputValues['browser_allsamp-gene'] === 'AMD1-ENSG00000123505' &&
      Shiny.shinyapp.$inputValues['browser_allsamp-tx'] === 'ENST00000451850' &&
      !document.documentElement.classList.contains('shiny-busy')).catch(async error => {
        console.log('Initialization:', await page.evaluate(() => ({
          gene:Shiny.shinyapp.$inputValues['browser_allsamp-gene'],
          tx:Shiny.shinyapp.$inputValues['browser_allsamp-tx'],
          busy:document.documentElement.className,
          errors:[...document.querySelectorAll('.shiny-output-error')].map(e=>e.textContent)})));
        throw error;
      });
    await page.evaluate(() => {
      const gene = document.getElementById('browser_allsamp-gene').selectize;
      gene.addOption({value:'ATF4-ENSG00000128272',label:'ATF4-ENSG00000128272'});
      gene.setValue('ATF4-ENSG00000128272');
    });
    await page.waitForFunction(() => Shiny.shinyapp.$inputValues['browser_allsamp-gene'] === 'ATF4-ENSG00000128272' &&
      document.getElementById('browser_allsamp-tx').value && !document.documentElement.classList.contains('shiny-busy'));
    await page.evaluate(() => {
      document.getElementById('browser_allsamp-tx').selectize.setValue('ENST00000674920');
      for (const id of ['add_translon','add_translons_transcode']) {
        const input = document.getElementById('browser_allsamp-'+id);
        input.checked = true; jQuery(input).trigger('change');
      }
    });
    await page.waitForFunction(() => Shiny.shinyapp.$inputValues['browser_allsamp-tx'] === 'ENST00000674920');
    await page.locator('#browser_allsamp-collapsed_translons').check();
    await page.locator('#browser_allsamp-collapsed_clusters').check();
    await page.locator('#browser_allsamp-go').click();
    await page.waitForFunction(() => document.getElementById('browser_allsamp-myPlotlyPlot')?._fullData?.[0]?.z.length === 5 &&
      document.getElementById('browser_allsamp-mb_bottom_gene')?._fullData?.[0]?.text[4] === 'clean CDS' && !document.documentElement.classList.contains('shiny-busy'));
    assert.deepEqual(await page.evaluate(() => document.getElementById('browser_allsamp-mb_bottom_gene')._fullData[0].text), ['uORF1','uORF2','uORF3','uORF4','clean CDS']);
    const coordinates = await page.evaluate(() => {
      const plot = document.getElementById('browser_allsamp-myPlotlyPlot');
      const rect = plot.getBoundingClientRect();
      return {x:rect.left+plot._fullLayout.xaxis._offset+plot._fullLayout.xaxis.d2p(4),
        y:rect.top+plot._fullLayout.yaxis._offset+plot._fullLayout.yaxis.d2p(2)};
    });
    const start = performance.now();
    await page.mouse.click(coordinates.x, coordinates.y);
    await page.waitForSelector('#browser_allsamp-cell_summary');
    await page.waitForFunction(() => document.getElementById('browser_allsamp-cell_distribution')?._fullData?.[0]?.type === 'box');
    console.log('Cell inspector ms:', performance.now()-start);
    assert.match(await page.locator('.modal-title').innerText(), /uORF4/);
    assert.equal(await page.locator('#browser_allsamp-cell_reference').inputValue(), 'R6');
    await page.evaluate(() => document.getElementById('browser_allsamp-cell_metric').selectize.setValue('Ratio'));
    await page.waitForFunction(() => {
      const title = document.getElementById('browser_allsamp-cell_distribution')?._fullLayout?.yaxis.title;
      return (typeof title === 'string' ? title : title?.text) === 'log2(region / reference density)';
    });
    await page.locator('a[data-value="Ratio by metadata"]').click();
    await page.evaluate(() => document.getElementById('browser_allsamp-cell_metadata').selectize.setValue('CONDITION'));
    await page.waitForFunction(() => document.querySelectorAll('#browser_allsamp-cell_ratio_metadata tbody tr').length > 0 && !document.documentElement.classList.contains('shiny-busy'));
    assert.match(await page.locator('#browser_allsamp-cell_ratio_metadata').innerText(), /Largest study/i);
    assert.match(await page.locator('#browser_allsamp-cell_ratio_metadata').innerText(), /Median_log2FC/i);
    await page.screenshot({path:'/tmp/ribocrypt-cell-details-desktop.png',fullPage:true});
    const download = page.waitForEvent('download');
    await page.locator('#browser_allsamp-cell_download').click();
    await (await download).saveAs('/tmp/ribocrypt-cell-libraries.csv');
    await page.setViewportSize({width:390,height:844});
    await page.evaluate(() => {document.getElementById('shiny-modal').scrollTop=0;window.scrollTo(0,0);});
    await page.waitForTimeout(750);
    await page.screenshot({path:'/tmp/ribocrypt-cell-details-mobile.png'});
    await page.locator('.modal-footer button').filter({hasText:'Close'}).click();
    await page.setViewportSize({width:1440,height:1100});
    const values = await page.evaluate(() => JSON.stringify(document.getElementById('browser_allsamp-myPlotlyPlot')._fullData[0].z));
    await page.evaluate(() => {
      const input = document.getElementById('browser_allsamp-region_min_coverage');
      input.value = '1000000000000'; input.dispatchEvent(new Event('input', {bubbles:true}));
      input.dispatchEvent(new Event('change', {bubbles:true}));
    });
    await page.waitForFunction(() => document.getElementById('browser_allsamp-myPlotlyPlot')._fullData.length === 2 &&
      document.getElementById('browser_allsamp-myPlotlyPlot')._fullData[1].x.length === 25, null, {timeout:10000}).catch(async error => {
        console.log('Support state:', await page.evaluate(() => ({
          inputs:Object.fromEntries(Object.entries(Shiny.shinyapp.$inputValues).filter(([key])=>key.includes('region_min'))),
          enabled:Shiny.shinyapp.$inputValues['browser_allsamp-region_support'],
          traces:document.getElementById('browser_allsamp-myPlotlyPlot')._fullData.map(t=>({type:t.type,n:t.x?.length})),
          data:document.getElementById('browser_allsamp-myPlotlyPlot').data.map(t=>({type:t.type,x:t.x,y:t.y})),
          errors:[...document.querySelectorAll('.shiny-output-error')].map(e=>e.textContent)})));
        throw error;
      });
    assert.equal(await page.evaluate(() => JSON.stringify(document.getElementById('browser_allsamp-myPlotlyPlot')._fullData[0].z)), values);
    await page.evaluate(() => {const input=document.getElementById('browser_allsamp-region_support');input.checked=false;jQuery(input).trigger('change');});
    await page.waitForFunction(() => document.getElementById('browser_allsamp-myPlotlyPlot')._fullData.length === 1);
    assert.deepEqual(errors, []);
    assert.deepEqual(await page.locator('.shiny-output-error').allTextContents(), []);
    const next = await browser.newPage();
    next.setDefaultTimeout(60000);
    await next.goto((process.argv[2] || 'http://127.0.0.1:7826/') + '#MegaBrowser');
    await next.waitForFunction(() => document.getElementById('browser_allsamp-gene')?.selectize &&
      document.getElementById('browser_allsamp-tx')?.selectize && !document.documentElement.classList.contains('shiny-busy'));
    await next.close();
    console.log('ATF4 cell inspection, clean-CDS ratios, metadata terms, export and support controls passed');
  } finally {await browser.close();}
})().catch(error => {console.error(error);process.exitCode=1;});
