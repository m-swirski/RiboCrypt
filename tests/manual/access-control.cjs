const http = require('http');
// Synthetic local gateway only. Cookie identities here are NOT authentication.
const assert = require('assert/strict');
const { chromium } = require('playwright');
const { execFileSync } = require('child_process');

const proxy = http.createServer((req, res) => {
  const headers = trustedHeaders(req);
  const upstream = http.request({ host: '127.0.0.1', port: 7837, path: req.url, method: req.method, headers }, reply => {
    res.writeHead(reply.statusCode, reply.headers); reply.pipe(res);
  });
  upstream.on('error', () => { res.writeHead(502); res.end(); });
  req.pipe(upstream);
});
function trustedHeaders(req) {
  const user = /rc_test_user=(alice|bob)/.exec(req.headers.cookie || '')?.[1] || '';
  return { ...req.headers, 'ribocrypt-gateway-secret': 't'.repeat(64), 'ribocrypt-subject': user };
}
proxy.on('upgrade', (req, socket, head) => {
  const upstream = http.request({ host: '127.0.0.1', port: 7837, path: req.url, headers: trustedHeaders(req) });
  upstream.on('upgrade', (res, upstreamSocket, upstreamHead) => {
    socket.write(`HTTP/1.1 101 Switching Protocols\r\n${Object.entries(res.headers).map(([k,v]) => `${k}: ${v}`).join('\r\n')}\r\n\r\n`);
    if (head.length) upstreamSocket.write(head);
    if (upstreamHead.length) socket.write(upstreamHead);
    upstreamSocket.pipe(socket); socket.pipe(upstreamSocket);
  });
  upstream.on('error', () => socket.destroy()); upstream.end();
});

(async () => {
  await new Promise(resolve => proxy.listen(7838, '127.0.0.1', resolve));
  const browser = await chromium.launch({ executablePath: '/opt/google/chrome/chrome', headless: true, args: ['--no-sandbox'] });
  try {
    for (const user of ['', 'alice', 'bob']) {
      const context = await browser.newContext({ viewport: { width: 1440, height: 1000 } });
      if (user) await context.addCookies([{ name: 'rc_test_user', value: user, url: 'http://127.0.0.1:7838' }]);
      const page = await context.newPage();
      const errors = [];
      page.on('pageerror', e => errors.push(e.message));
      await page.goto('http://127.0.0.1:7838/', { waitUntil: 'domcontentloaded' });
      await page.waitForFunction(() => window.Shiny?.shinyapp && document.querySelector('#browser-dff')?.selectize, null, { timeout: 90000 });
      await page.waitForFunction(() => document.querySelector('#browser-dff')?.selectize.getValue() === 'human_all_merged_l50', null, { timeout: 60000 });
      await page.waitForFunction(() => typeof document.querySelector('#browser-dff')?.selectize.settings.load === 'function', null, { timeout: 60000 });
      await page.evaluate(() => { const s = document.querySelector('#browser-dff').selectize; s.loadedSearches = {}; s.open(); s.onSearchChange(''); });
      await page.waitForFunction(() => Object.keys(document.querySelector('#browser-dff').selectize.options).includes('human_all_merged_l50'), null, { timeout: 60000 });
      const options = await page.evaluate(() => Object.keys(document.querySelector('#browser-dff').selectize.options));
      console.log(JSON.stringify({ user, options }));
      assert(options.includes('human_all_merged_l50'));
      assert.equal(options.includes('alice_private'), user === 'alice');
      assert.equal(options.includes('bob_private'), user === 'bob');
      await page.waitForFunction(() => document.querySelectorAll('#browser-c .main-svg').length > 0 || document.querySelectorAll('.js-plotly-plot .main-svg').length > 0, null, { timeout: 90000 });
      if (user) {
        await page.evaluate(selected => document.getElementById('browser-dff').selectize.setValue(selected), `${user}_private`);
        await page.waitForTimeout(1000);
        await page.evaluate(() => {
          window.rcPlotUpdates = 0;
          $(document).on('shiny:value.rcPlotTest', event => { if (event.name === 'browser-c') window.rcPlotUpdates++; });
        });
        await page.locator('#browser-go').click();
        await page.waitForFunction(() => window.rcPlotUpdates > 0, null, { timeout: 30000 });
        assert.equal(await page.locator('.shiny-output-error').count(), 0);
      }
      await page.locator('#browser-toggle_settings').click();
      await page.waitForFunction(() => !document.querySelector('#browser-floating_settings').classList.contains('hidden'), null, { timeout: 30000 });
      await page.evaluate(() => {
        const pane = document.getElementById('browser-download_coverage').closest('.tab-pane');
        document.querySelector(`#browser-floating_settings a[href="#${pane.id}"]`).click();
      });
      await page.waitForFunction(() => document.querySelector('#browser-download_coverage')?.getAttribute('href')?.includes('/download/'));
      const download = new URL(await page.locator('#browser-download_coverage').getAttribute('href'), page.url()).href;
      const exported = await context.request.get(download);
      assert.equal(exported.status(), user === 'alice' ? 403 : 200);
      if (user !== 'alice') assert((await exported.text()).startsWith('position,'));
      const other = await browser.newContext();
      await other.addCookies([{ name: 'rc_test_user', value: user === 'bob' ? 'alice' : 'bob', url: 'http://127.0.0.1:7838' }]);
      assert.equal((await other.request.get(download)).status(), 403);
      await other.close();
      await page.locator('#browser-toggle_settings').click();
      await page.screenshot({ path: `/tmp/ribocrypt-account/${user || 'anonymous'}.png`, fullPage: true });
      await page.evaluate(() => document.querySelector('#navbarID a[data-value="Samples"]').click());
      await page.waitForFunction(() => window.jQuery?.fn.dataTable?.isDataTable(document.querySelector('#sample_info-sample_info table')), null, { timeout: 30000 });
      await page.evaluate(() => $('#sample_info-sample_info table').DataTable().search('ONLY_VISIBLE_TO_').draw());
      await page.waitForFunction(expected => $('#sample_info-sample_info table').DataTable().page.info().recordsDisplay === expected, user ? 1 : 0, { timeout: 30000 });
      const metadata = await page.evaluate(() => JSON.stringify($('#sample_info-sample_info table').DataTable().rows().data().toArray()));
      if (user) assert(metadata.includes(`ONLY_VISIBLE_TO_${user.toUpperCase()}`));
      assert(!metadata.includes(`ONLY_VISIBLE_TO_${user === 'alice' ? 'BOB' : 'ALICE'}`));
      await page.evaluate(() => document.querySelector('#navbarID a[data-value="Observatory"]').click());
      await page.waitForFunction(() => {
        const table = document.querySelector('#observatory-selector-libraries_data_table table[id]');
        return table && $.fn.dataTable.isDataTable(table) && $(table).DataTable().page.info()?.recordsTotal > 0;
      }, null, { timeout: 60000 });
      const libraryCount = await page.evaluate(() => $('#observatory-selector-libraries_data_table table[id]').DataTable().page.info().recordsTotal);
      assert.equal(libraryCount, 3857);
      await page.waitForFunction(() => document.querySelector('#observatory-selector-libraries_umap_plot .main-svg'), null, { timeout: 30000 });
      assert.deepEqual(errors, []);
      console.log(JSON.stringify({ user: user || 'anonymous', options, downloadStatus: exported.status(), metadata, libraryCount, errors }));
      if (user === 'bob') {
        execFileSync('sqlite3', ['/tmp/ribocrypt-account/data/access.sqlite', "DELETE FROM grants WHERE workspace='bob' AND dataset='bob_private';"]);
        await page.waitForEvent('domcontentloaded', { timeout: 15000 });
        await page.waitForFunction(() => typeof document.querySelector('#browser-dff')?.selectize.settings.load === 'function', null, { timeout: 60000 });
        await page.evaluate(() => { const s = document.querySelector('#browser-dff').selectize; s.loadedSearches = {}; s.onSearchChange(''); });
        await page.waitForFunction(() => Object.keys(document.querySelector('#browser-dff').selectize.options).includes('human_all_merged_l50'), null, { timeout: 30000 });
        assert(!(await page.evaluate(() => Object.keys(document.querySelector('#browser-dff').selectize.options))).includes('bob_private'));
        console.log('Active-session revocation: browser reloaded and private dataset removed.');
        execFileSync('sqlite3', ['/tmp/ribocrypt-account/data/access.sqlite', "INSERT INTO grants VALUES ('bob','bob_private',1);"]);
      }
      await context.close();
    }
  } finally { await browser.close(); proxy.close(); }
})().catch(e => { console.error(e); process.exitCode = 1; proxy.close(); });
