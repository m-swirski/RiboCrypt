// Loopback-only regression test; identities are synthetic, not a login provider.
const assert = require('assert/strict');
const fs = require('fs');
const path = require('path');
const { spawn, execFileSync } = require('child_process');
const { chromium } = require('playwright');
const root = fs.mkdtempSync('/tmp/ribocrypt-report-');
const origin = 'http://127.0.0.1:7849';
const headers = { 'RiboCrypt-Gateway-Secret': 'x'.repeat(64), 'RiboCrypt-Subject': 'alice' };
const log = fs.openSync(path.join(root, 'server.log'), 'w', 0o600);
const server = spawn('env', ['-u', 'LC_ALL', 'R', '--vanilla', '-q', '-f', 'tests/manual/fastq-report-launch.R'], {
  env: { ...process.env, RIBOCRYPT_REPORT_TEST_ROOT: root }, stdio: ['ignore', log, log] });
let browser;
(async () => {
  try {
    const deadline = Date.now() + 120000;
    while (true) {
      try { if ((await fetch(origin)).ok) break; } catch (_) {}
      assert(Date.now() < deadline, fs.readFileSync(path.join(root, 'server.log'), 'utf8'));
      await new Promise(r => setTimeout(r, 500));
    }
    browser = await chromium.launch({ executablePath: '/opt/google/chrome/chrome', headless: true,
      args: ['--no-sandbox'] });
    const page = await browser.newPage({ extraHTTPHeaders: headers });
    await page.context().addCookies([{ name: 'rc_report_test', value: 'alice', url: origin }]);
    await page.goto(origin);
    const iframe = page.locator('#report-frame');
    await iframe.waitFor();
    const url = new URL(await iframe.getAttribute('src'), origin);
    const response = await fetch(url, { headers });
    assert.equal(response.status, 200);
    assert.match(await response.text(), /SELECTED_REPORT/);
    assert.equal(response.headers.get('cache-control'), 'no-store');
    await page.frameLocator('#report-frame').locator('body[data-script="works"]').waitFor();
    for (const user of [null, 'bob']) {
      const altered = { 'RiboCrypt-Gateway-Secret': headers['RiboCrypt-Gateway-Secret'] };
      if (user) altered['RiboCrypt-Subject'] = user;
      assert.equal((await fetch(url, { headers: altered })).status, 403);
    }
    assert.equal((await fetch(url)).status, 403);
    for (const suffix of ['/tmpuser/sample.html', '/tmpuser/sibling.txt'])
      assert.equal((await fetch(origin + suffix, { headers })).status, 404);
    url.searchParams.set('file', 'sibling.txt');
    assert.doesNotMatch(await (await fetch(url, { headers })).text(), /PRIVATE_SIBLING/);
    execFileSync('sqlite3', [path.join(root, 'reports.sqlite'), 'DELETE FROM grants;']);
    assert.equal((await fetch(url, { headers })).status, 403);
    console.log('PASS real FASTQ selection/session endpoint, read-only access, sandbox scripts, identity checks, no static siblings, revocation');
  } catch (error) {
    console.error(fs.readFileSync(path.join(root, 'server.log'), 'utf8'));
    throw error;
  } finally {
    if (browser) await browser.close();
    const exited = new Promise(resolve => server.once('exit', resolve));
    server.kill('SIGINT');
    const timer = setTimeout(() => server.kill('SIGTERM'), 5000);
    await exited; clearTimeout(timer); fs.closeSync(log);
    fs.rmSync(root, { recursive: true, force: true });
  }
})().catch(error => { console.error(error); process.exitCode = 1; });
