const {chromium} = require('playwright');
const assert = require('node:assert/strict');
const base = process.argv[2] || 'http://127.0.0.1:7821/';
async function settle(page) {
  await page.waitForTimeout(600);
  await page.waitForFunction(() => !document.documentElement.classList.contains('shiny-busy'));
}
async function snapshot(page, id) {
  return page.evaluate(id => {
    const el = document.getElementById(id);
    return {x: el._fullLayout.xaxis.range, y: el._fullLayout.yaxis.range,
      bottom: el._fullLayout.yaxis2.range, traces: el.data.length};
  }, id);
}
(async () => {
  const browser = await chromium.launch({executablePath: '/opt/google/chrome/chrome',
    headless: true, args: ['--no-sandbox']});
  try {
    for (const obs of [false, true]) {
      const context = await browser.newContext({viewport: {width: 1440, height: 1000}});
      const page = await context.newPage(), errors = [];
      page.on('pageerror', error => errors.push(error.message));
      let copied;
      page.on('websocket', socket => socket.on('framereceived', event => {
        try {copied = JSON.parse(event.payload).custom?.['ribocrypt-copy-url']?.text || copied;} catch {}
      }));
      const prefix = obs ? 'observatory-browser_obs-' : 'browser-';
      const id = prefix + (obs ? 'browser_plot' : 'c');
      await page.goto(base + (obs ? '#Observatory' : '#browser'));
      if (obs) {
        await page.locator('a[data-value="Browse"]').click();
        await page.waitForFunction(p => document.getElementById(p + 'tx')?.value, prefix, {timeout: 120000});
        await page.locator('#' + prefix + 'go').click();
      }
      await page.waitForFunction(id => !!document.getElementById(id)?._localYCleanup, id, {timeout: 120000});
      await settle(page);
      const initial = await snapshot(page, id);
      assert.equal(await page.locator('#' + prefix + 'local_y_max').isChecked(), true);
      const window = await page.evaluate(id => {
        const el = document.getElementById(id), full = el._fullLayout.xaxis.range;
        const points = el.data.filter(t => (t.yaxis || 'y') === 'y' && t.visible !== 'legendonly')
          .flatMap(t => Array.from(t.x || [], (x, i) => [x, t.y[i]]))
          .filter(([x, y]) => Number.isFinite(x) && Number.isFinite(y));
        let best = null;
        const width = Math.max(20, Math.floor((full[1] - full[0]) / 20));
        for (let low = Math.ceil(full[0]); low + width < full[1]; low += width) {
          const peak = Math.max(0, ...points.filter(([x]) => x >= low && x <= low + width).map(p => p[1]));
          if (peak > 0 && (!best || peak < best.peak)) best = {range: [low, low + width], peak};
        }
        return best;
      }, id);
      assert.ok(window);
      await page.evaluate(({id, range}) => Plotly.relayout(document.getElementById(id), {'xaxis.range': range}), {id, range: window.range});
      await settle(page);
      const zoomed = await snapshot(page, id);
      assert.ok(zoomed.y[1] < initial.y[1], JSON.stringify({initial, zoomed, window}));
      assert.ok(zoomed.y[1] >= window.peak);
      assert.deepEqual(zoomed.bottom, initial.bottom);
      // Exercise the live setting through the actual checkbox, without generating a new plot.
      await page.locator('#' + prefix + 'toggle_settings').click();
      await page.locator('#' + prefix + 'floating_settings a[data-value="Settings"]').click();
      await page.locator('#' + prefix + 'local_y_max').uncheck();
      await settle(page);
      assert.deepEqual((await snapshot(page, id)).y, initial.y);
      await page.locator('#' + prefix + 'clip_button').click();
      for (let attempt = 0; !copied && attempt < 100; attempt++) await page.waitForTimeout(100);
      assert.ok(copied);
      const shared = await context.newPage();
      await shared.goto(copied);
      await shared.waitForFunction(id => !!document.getElementById(id)?._localYCleanup, id, {timeout: 120000});
      await settle(shared);
      assert.equal(await shared.locator('#' + prefix + 'local_y_max').isChecked(), false);
      const sharedInitial = await snapshot(shared, id);
      await shared.evaluate(({id, range}) => Plotly.relayout(document.getElementById(id), {'xaxis.range': range}), {id, range: window.range});
      await settle(shared);
      assert.deepEqual((await snapshot(shared, id)).y, sharedInitial.y);
      await shared.close();
      await page.locator('#' + prefix + 'local_y_max').check();
      await settle(page);
      assert.deepEqual((await snapshot(page, id)).y, zoomed.y);
      await page.locator('#' + prefix + 'toggle_settings').click();
      await page.screenshot({path: '/tmp/local-y-' + (obs ? 'observatory' : 'browser') + '.png'});
      await page.evaluate(id => Plotly.relayout(document.getElementById(id), {'xaxis.autorange': true}), id);
      await settle(page);
      assert.deepEqual((await snapshot(page, id)).y, initial.y);
      assert.deepEqual(errors, []);
      assert.deepEqual(await page.locator('.shiny-output-error').allTextContents(), []);
      console.log(JSON.stringify({obs, initial: initial.y, zoomed: zoomed.y, window, reset: true, liveToggle: true, url: true, errors}));
      await context.close();
    }
  } finally {await browser.close();}
})().catch(error => {console.error(error); process.exitCode = 1;});
