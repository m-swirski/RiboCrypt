// Disposable loopback-only OIDC integration test; never deploy this harness.
const fs = require('fs');
const path = require('path');
const net = require('net');
const http = require('http');
const crypto = require('crypto');
const assert = require('assert/strict');
const { execFileSync, spawn } = require('child_process');
const { chromium } = require('playwright');
const { fetch, Agent, setGlobalDispatcher } = require('undici');
require('dns').setDefaultResultOrder('ipv4first');
const root = process.env.RIBOCRYPT_AUTH_TEST_ROOT;
assert(root && path.isAbsolute(root), 'Set RIBOCRYPT_AUTH_TEST_ROOT on a drive with sufficient space.');
const runroot = process.env.RIBOCRYPT_AUTH_TEST_RUNROOT || `/run/user/${process.getuid()}/ribocrypt-auth-overlay`;
const fixtures = path.join(root, 'fixtures');
fs.mkdirSync(fixtures, { recursive: true, mode: 0o700 });
const storage = process.env.RIBOCRYPT_AUTH_TEST_STORAGE || path.join(root, 'storage');
const podArgs = ['--root', storage, '--runroot', runroot, '--storage-driver', 'overlay',
  '--storage-opt', 'overlay.mount_program=/usr/bin/fuse-overlayfs'];
const containers = [];
const secret = () => crypto.randomBytes(32).toString('hex');
const adminPassword = secret(), dbPassword = secret(), clientSecret = secret(), gatewaySecret = secret();
const issuer = 'https://localhost:8444/realms/ribocrypt';
const origin = 'https://localhost:8443';
const password = `Test-${secret()}`;
const messages = [];
let browser, backend;
function write(name, value) { fs.writeFileSync(path.join(fixtures, name), value, { mode: 0o600 }); }
function pod(...args) { return execFileSync('podman', [...podArgs, ...args], { encoding: 'utf8', maxBuffer: 16 * 1024 * 1024, stdio: ['ignore', 'pipe', 'pipe'] }); }
async function wait(check, label, timeout = 180000) {
  const end = Date.now() + timeout;
  let lastError;
  while (Date.now() < end) { try { const result = await check(); if (result) return result; } catch (error) { lastError = error; } await new Promise(r => setTimeout(r, 500)); }
  throw new Error(`Timed out: ${label}; ${lastError?.cause || lastError || 'not ready'}`);
}
async function api(url, options = {}) {
  const response = await fetch(url, options);
  assert(response.ok, `${url}: ${response.status} ${await response.clone().text()}`);
  return response.status === 204 || !response.headers.get('content-type')?.includes('json') ? response : response.json();
}
async function ready(url, options = {}) {
  const response = await fetch(url, { ...options, signal: AbortSignal.timeout(5000) });
  await response.arrayBuffer();
  return response.ok;
}
// Test SMTP sink binds only to loopback and sends no external mail.
const smtp = net.createServer(socket => {
  socket.setEncoding('utf8'); socket.write('220 localhost test SMTP\r\n');
  let buffer = '', data = false, message = [];
  socket.on('data', chunk => {
    buffer += chunk;
    while (buffer.includes('\r\n')) {
      const index = buffer.indexOf('\r\n'), line = buffer.slice(0, index); buffer = buffer.slice(index + 2);
      if (data) {
        if (line !== '.') { message.push(line.replace(/^\.\./, '.')); continue; }
        messages.push(message.join('\r\n')); data = false; message = []; socket.write('250 received\r\n');
      } else if (/^EHLO|^HELO/i.test(line)) socket.write('250 localhost\r\n');
      else if (/^DATA/i.test(line)) { data = true; socket.write('354 send message\r\n'); }
      else if (/^QUIT/i.test(line)) socket.end('221 goodbye\r\n');
      else socket.write('250 OK\r\n');
    }
  });
});
const probe = http.createServer((req, res) => {
  res.setHeader('content-type', 'application/json');
  res.end(JSON.stringify({ subject: req.headers['ribocrypt-subject'] || '',
    gateway: req.headers['ribocrypt-gateway-secret'] === gatewaySecret,
    authorization: req.headers.authorization || '' }));
});
async function emailLink(since) {
  const raw = await wait(() => messages.slice(since).find(x => x.includes('localhost:8444')), 'email delivery');
  const decoded = raw.replace(/=\r\n/g, '').replace(/=([0-9A-F]{2})/gi, (_, hex) => String.fromCharCode(parseInt(hex, 16)));
  const link = decoded.match(/https:\/\/localhost:8444\/[^\s"<>]+/);
  assert(link, 'Email contains action link'); return link[0].replace(/&amp;/g, '&');
}
async function login(page, email, pass = password) {
  await page.goto(`${origin}/oauth2/start?rd=/`);
  await page.locator('#username').fill(email);
  await page.locator('#password').fill(pass);
  await page.locator('#kc-login').click();
}
async function run() {
  await new Promise(r => smtp.listen(11025, '127.0.0.1', r));
  await new Promise(r => probe.listen(17837, '127.0.0.1', r));
  execFileSync('openssl', ['req', '-x509', '-newkey', 'rsa:2048', '-nodes', '-days', '2',
    '-subj', '/CN=localhost', '-addext', 'subjectAltName=DNS:localhost,IP:127.0.0.1',
    '-keyout', path.join(fixtures, 'key.pem'), '-out', path.join(fixtures, 'cert.pem')], { stdio: 'ignore' });
  fs.chmodSync(path.join(fixtures, 'key.pem'), 0o600);
  const realm = JSON.parse(fs.readFileSync('deploy/auth/realm.json', 'utf8'));
  realm.smtpServer = { host: '127.0.0.1', port: '11025', from: 'test@ribocrypt.invalid', ssl: 'false', starttls: 'false', auth: 'false' };
  write('realm.json', JSON.stringify(realm));
  // CA trust is scoped to this child/test process, never installed system-wide.
  require('https').globalAgent.options.ca = fs.readFileSync(path.join(fixtures, 'cert.pem'));
  setGlobalDispatcher(new Agent({ connect: { ca: fs.readFileSync(path.join(fixtures, 'cert.pem')) } }));
  let shared = fs.readFileSync('deploy/auth/shiny-proxy.conf.example', 'utf8')
    .replace('127.0.0.1:3838', '127.0.0.1:17837').replaceAll('/etc/nginx/ribocrypt/', '/fixtures/');
  write('shiny-proxy.conf', shared);
  write('gateway-secret.conf', `proxy_set_header RiboCrypt-Gateway-Secret "${gatewaySecret}";\n`);
  const app = fs.readFileSync('deploy/auth/nginx-app.conf.example', 'utf8').replaceAll('/etc/nginx/ribocrypt/', '/fixtures/').replaceAll('127.0.0.1:4180', '127.0.0.1:14180');
  const loginConfig = fs.readFileSync('deploy/auth/nginx-login.conf.example', 'utf8')
    .replaceAll('127.0.0.1:8180', '127.0.0.1:18180')
    .replaceAll('$host', '$http_host').replace('X-Forwarded-Port 443', 'X-Forwarded-Port 8444');
  write('nginx.conf', `pid /tmp/ribocrypt-nginx.pid; events {} http {
    access_log off; error_log /dev/stderr warn;
    map $http_upgrade $ribocrypt_connection_upgrade { default upgrade; '' close; }
    server { listen 127.0.0.1:8443 ssl; server_name localhost;
      ssl_certificate /fixtures/cert.pem; ssl_certificate_key /fixtures/key.pem;
      add_header X-RiboCrypt-Test-Security preserved always;
      ${app}
    }
    server { listen 127.0.0.1:8444 ssl; server_name localhost;
      ssl_certificate /fixtures/cert.pem; ssl_certificate_key /fixtures/key.pem;
      ${loginConfig}
    }
  }`);
  write('oauth2-proxy.cfg', fs.readFileSync('deploy/auth/oauth2-proxy.cfg', 'utf8')
    .replace('0.0.0.0:4180', '127.0.0.1:14180') + `\nprovider_ca_files = ["/fixtures/cert.pem"]\n`);
  containers.push('rc-auth-postgres');
  pod('run', '-d', '--name', 'rc-auth-postgres', '--network', 'host', '-e', 'POSTGRES_DB=keycloak', '-e', 'POSTGRES_USER=keycloak', '-e', `POSTGRES_PASSWORD=${dbPassword}`,
    'docker.io/library/postgres:17', '-p', '15432', '-c', 'listen_addresses=127.0.0.1');
  await wait(() => pod('exec', 'rc-auth-postgres', 'pg_isready', '-p', '15432').includes('accepting'), 'PostgreSQL');
  containers.push('rc-auth-keycloak');
  pod('run', '-d', '--name', 'rc-auth-keycloak', '--user', '0', '--network', 'host', '-v', `${fixtures}/realm.json:/opt/keycloak/data/import/realm.json:ro`,
    '-e', 'KC_DB=postgres', '-e', 'KC_DB_URL=jdbc:postgresql://127.0.0.1:15432/keycloak', '-e', 'KC_DB_USERNAME=keycloak', '-e', `KC_DB_PASSWORD=${dbPassword}`,
    '-e', 'KC_HOSTNAME=https://localhost:8444', '-e', 'KC_HTTP_ENABLED=true', '-e', 'KC_HTTP_HOST=127.0.0.1', '-e', 'KC_HTTP_MANAGEMENT_HOST=127.0.0.1', '-e', 'KC_HTTP_PORT=18180', '-e', 'KC_PROXY_HEADERS=xforwarded',
    '-e', 'KC_BOOTSTRAP_ADMIN_USERNAME=admin', '-e', `KC_BOOTSTRAP_ADMIN_PASSWORD=${adminPassword}`,
    '-e', 'RIBOCRYPT_APP_HOST=localhost:8443', '-e', `RIBOCRYPT_OIDC_CLIENT_SECRET=${clientSecret}`,
    'quay.io/keycloak/keycloak:26.8.0', 'start', '--import-realm', '--cache=local');
  containers.push('rc-auth-nginx');
  pod('run', '-d', '--name', 'rc-auth-nginx', '--network', 'host', '-v', `${fixtures}:/fixtures:ro`, 'docker.io/library/nginx:1.28', 'nginx', '-c', '/fixtures/nginx.conf', '-g', 'daemon off;');
  console.log('Waiting for production-mode Keycloak/PostgreSQL startup');
  await wait(() => ready(`${issuer}/.well-known/openid-configuration`), 'OIDC discovery');
  pod('exec', 'rc-auth-nginx', 'nginx', '-t', '-c', '/fixtures/nginx.conf');
  for (const route of ['/admin', '/admin/', '/admin/master/console/', '/realms/master/', '/%61dmin/'])
    assert.equal((await fetch(`https://localhost:8444${route}`)).status, 404, route);
  assert.equal((await fetch(origin)).headers.get('x-ribocrypt-test-security'), 'preserved');
  containers.push('rc-auth-gateway');
  pod('run', '-d', '--name', 'rc-auth-gateway', '--user', '0', '--network', 'host', '-v', `${fixtures}:/fixtures:ro`,
    '-e', `OAUTH2_PROXY_OIDC_ISSUER_URL=${issuer}`, '-e', `OAUTH2_PROXY_REDIRECT_URL=${origin}/oauth2/callback`, '-e', `OAUTH2_PROXY_CLIENT_SECRET=${clientSecret}`,
    '-e', `OAUTH2_PROXY_COOKIE_SECRET=${crypto.randomBytes(32).toString('base64url')}`,
    '-e', `OAUTH2_PROXY_BACKEND_LOGOUT_URL=${issuer}/protocol/openid-connect/logout?id_token_hint={id_token}`,
    'quay.io/oauth2-proxy/oauth2-proxy:v7.15.0', '--config=/fixtures/oauth2-proxy.cfg');
  await wait(() => ready('http://127.0.0.1:14180/ping'), 'OAuth2 Proxy', 30000);
  const anonymous = await api(origin, { headers: { 'RiboCrypt-Subject': 'forged', 'RiboCrypt-Gateway-Secret': 'forged', Authorization: 'Bearer forged' } });
  assert.deepEqual(anonymous, { subject: '', gateway: true, authorization: '' });
  assert.equal((await fetch(`${origin}/_ribocrypt_auth`)).status, 404);
  console.log('PASS anonymous access, spoofed identity stripped, internal auth endpoint inaccessible');
  browser = await chromium.launch({ executablePath: '/opt/google/chrome/chrome', headless: true, args: ['--no-sandbox'] });
  const context = await browser.newContext({ ignoreHTTPSErrors: true, viewport: { width: 1440, height: 1000 } });
  const page = await context.newPage();
  await page.goto(`${origin}/oauth2/start?rd=/`);
  await page.getByRole('link', { name: 'Register', exact: true }).click();
  write('registration.html', await page.content());
  await page.locator('#email').fill('alice@ribocrypt.invalid');
  if (await page.locator('#firstName').count()) await page.locator('#firstName').fill('Alice');
  if (await page.locator('#lastName').count()) await page.locator('#lastName').fill('Test');
  const passwordBeforeVerification = await page.locator('#password').count();
  if (passwordBeforeVerification) {
    await page.locator('#password').fill(password);
    await page.locator('#password-confirm').fill(password);
  }
  const before = messages.length;
  await page.getByRole('button', { name: 'Register', exact: true }).click();
  await page.getByText('You need to verify your email address to activate your account.').waitFor();
  assert.equal((await context.request.get(`${origin}/oauth2/auth`)).status(), 401);
  await page.goto(await emailLink(before));
  if (!passwordBeforeVerification) {
    await page.locator('#password-new').fill('short');
    await page.locator('#password-confirm').fill('short');
    await page.getByRole('button', { name: 'Submit', exact: true }).click();
    await page.getByText(/minimum length 12/i).waitFor();
    await page.locator('#password-new').fill(password);
    await page.locator('#password-confirm').fill(password);
    await page.getByRole('button', { name: 'Submit', exact: true }).click();
  }
  await page.waitForURL(`${origin}/`, { timeout: 60000 });
  const identity = JSON.parse(await page.locator('body').innerText());
  assert.equal((await context.request.get(origin)).headers()['x-ribocrypt-test-security'], 'preserved');
  assert(identity.subject && identity.subject !== 'alice@ribocrypt.invalid'); assert(identity.gateway);
  const subject = identity.subject;
  const cookies = await context.cookies(origin);
  const sessionCookie = cookies.find(c => c.name === '__Host-ribocrypt');
  assert(sessionCookie?.secure && sessionCookie.httpOnly && sessionCookie.path === '/' && sessionCookie.sameSite === 'Lax');
  const tampered = await browser.newContext({ ignoreHTTPSErrors: true });
  await tampered.addCookies([{ ...sessionCookie, value: `${sessionCookie.value}forged` }]);
  assert.equal((await tampered.request.get(`${origin}/oauth2/auth`)).status(), 401);
  const forgedCallback = await tampered.request.get(`${origin}/oauth2/callback?code=forged&state=forged`);
  assert([400, 403, 500].includes(forgedCallback.status()));
  assert.equal((await tampered.request.get(`${origin}/oauth2/auth`)).status(), 401);
  await tampered.close();
  console.log('PASS account registration, verification mail/action, PKCE login, immutable subject and secure cookie');
  await page.goto(`${origin}/ribocrypt/logout`);
  await page.waitForURL(`${origin}/`);
  assert.equal(JSON.parse(await page.locator('body').innerText()).subject, '');
  await page.goto(`${origin}/oauth2/start?rd=/`);
  await page.locator('#username').waitFor();
  console.log('PASS logout clears gateway cookie and Keycloak SSO');
  await page.getByRole('link', { name: 'Forgot Password?', exact: true }).click();
  await page.locator('#username').fill('alice@ribocrypt.invalid');
  const resetBefore = messages.length;
  await page.getByRole('button', { name: 'Submit', exact: true }).click();
  await page.goto(await emailLink(resetBefore));
  const newPassword = `New-${secret()}`;
  await page.locator('#password-new').fill(newPassword);
  await page.locator('#password-confirm').fill(newPassword);
  await page.getByRole('button', { name: 'Submit', exact: true }).click();
  // A reset may authenticate immediately; explicitly log out before checking credentials.
  await page.goto(`${origin}/ribocrypt/logout`);
  await login(page, 'alice@ribocrypt.invalid', password);
  await page.getByText('Invalid username or password.', { exact: true }).waitFor();
  await page.locator('#password').fill(newPassword); await page.locator('#kc-login').click();
  await page.waitForURL(`${origin}/`);
  assert.equal(JSON.parse(await page.locator('body').innerText()).subject, subject);
  console.log('PASS password reset mail/action; old password denied, new password accepted');
  const backup = pod('exec', 'rc-auth-postgres', 'pg_dump', '-U', 'keycloak', '-p', '15432', '--no-owner', '--no-privileges', 'keycloak');
  write('identity-backup.sql', backup);
  pod('exec', 'rc-auth-postgres', 'createdb', '-U', 'keycloak', '-p', '15432', 'keycloak_restore');
  execFileSync('podman', [...podArgs, 'exec', '-i', 'rc-auth-postgres', 'psql', '-U', 'keycloak', '-p', '15432', '-v', 'ON_ERROR_STOP=1', '-d', 'keycloak_restore'],
    { input: backup, encoding: 'utf8', maxBuffer: 16 * 1024 * 1024 });
  const restored = pod('exec', 'rc-auth-postgres', 'psql', '-U', 'keycloak', '-p', '15432', '-d', 'keycloak_restore', '-At', '-c',
    "SELECT id || '|' || email_verified FROM user_entity WHERE email='alice@ribocrypt.invalid';").trim();
  assert.equal(restored, `${subject}|true`);
  pod('restart', 'rc-auth-keycloak');
  await wait(() => ready(`${issuer}/.well-known/openid-configuration`), 'Keycloak restart');
  await page.goto(`${origin}/ribocrypt/logout`);
  await login(page, 'alice@ribocrypt.invalid', newPassword);
  await page.waitForURL(`${origin}/`);
  assert.equal(JSON.parse(await page.locator('body').innerText()).subject, subject);
  console.log('PASS PostgreSQL backup/restore and provider restart preserve verified account and credentials');
  const log = fs.openSync(path.join(fixtures, 'app.log'), 'w', 0o600);
  backend = spawn('env', ['-u', 'LC_ALL', 'R', '--vanilla', '-q', '-f', 'tests/manual/access-control-launch.R'], {
    env: { ...process.env, RIBOCRYPT_TEST_ISSUER: issuer, RIBOCRYPT_TEST_ALICE_SUBJECT: subject, RIBOCRYPT_TEST_GATEWAY_SECRET: gatewaySecret }, stdio: ['ignore', log, log]
  });
  await wait(() => ready('http://127.0.0.1:7837', { headers: { 'RiboCrypt-Gateway-Secret': gatewaySecret } }), 'RiboCrypt startup');
  write('shiny-proxy.conf', shared.replace('127.0.0.1:17837', '127.0.0.1:7837'));
  pod('exec', 'rc-auth-nginx', 'nginx', '-s', 'reload', '-c', '/fixtures/nginx.conf');
  await wait(async () => {
    const response = await context.request.get(origin);
    return response.status() === 200 && (await response.text()).includes('id="authorized_app"');
  }, 'NGINX app reload', 15000);
  await page.goto(`${origin}/?rc_auth_test=${Date.now()}`);
  await page.waitForFunction(() => typeof document.querySelector('#browser-dff')?.selectize.settings.load === 'function', null, { timeout: 90000 });
  await page.evaluate(() => { const s = document.querySelector('#browser-dff').selectize; s.loadedSearches = {}; s.onSearchChange(''); });
  await page.waitForFunction(() => Object.keys(document.querySelector('#browser-dff').selectize.options).includes('alice_private'));
  assert(!(await page.evaluate(() => Object.keys(document.querySelector('#browser-dff').selectize.options))).includes('bob_private'));
  await page.waitForFunction(() => document.querySelector('#browser-c .main-svg'), null, { timeout: 90000 });
  assert.equal((await fetch('http://127.0.0.1:7837')).status, 403);
  await page.evaluate(() => {
    window.rcPlotUpdates = 0;
    $(document).on('shiny:value.rcAuthTest', event => { if (event.name === 'browser-c') window.rcPlotUpdates++; });
    document.getElementById('browser-dff').selectize.setValue('alice_private');
  });
  await page.waitForTimeout(1000);
  await page.locator('#browser-go').click();
  await page.waitForFunction(() => window.rcPlotUpdates > 0, null, { timeout: 30000 });
  assert.equal(await page.locator('.shiny-output-error').count(), 0);
  await page.screenshot({ path: path.join(fixtures, 'authenticated-app.png'), fullPage: true });
  console.log('PASS real OIDC identity reaches RiboCrypt over WebSocket; private catalog isolated, ATF4 rendered, untrusted backend request denied');
  await page.locator('#browser-toggle_settings').click();
  await page.waitForFunction(() => !document.querySelector('#browser-floating_settings').classList.contains('hidden'));
  await page.evaluate(() => {
    const pane = document.getElementById('browser-download_coverage').closest('.tab-pane');
    document.querySelector(`#browser-floating_settings a[href="#${pane.id}"]`).click();
  });
  await page.waitForFunction(() => document.querySelector('#browser-download_coverage')?.getAttribute('href')?.includes('/download/'));
  const download = new URL(await page.locator('#browser-download_coverage').getAttribute('href'), origin).href;
  assert.equal((await context.request.get(download)).status(), 403);
  const anonymousContext = await browser.newContext({ ignoreHTTPSErrors: true });
  assert.equal((await anonymousContext.request.get(download)).status(), 403);
  await anonymousContext.close();
  const reloaded = page.waitForEvent('domcontentloaded', { timeout: 45000 });
  const logoutTab = await context.newPage();
  await logoutTab.goto(`${origin}/ribocrypt/logout`);
  await reloaded;
  await page.waitForFunction(() => typeof document.querySelector('#browser-dff')?.selectize.settings.load === 'function', null, { timeout: 60000 });
  await page.evaluate(() => { const s = document.querySelector('#browser-dff').selectize; s.loadedSearches = {}; s.onSearchChange(''); });
  await page.waitForFunction(() => Object.keys(document.querySelector('#browser-dff').selectize.options).includes('human_all_merged_l50'));
  assert(!(await page.evaluate(() => Object.keys(document.querySelector('#browser-dff').selectize.options))).includes('alice_private'));
  console.log('PASS read-only CSV denial, anonymous session URL denial, logout in another tab clears active private session');
  pod('stop', 'rc-auth-gateway');
  assert.equal((await context.request.get(origin)).status(), 500);
  console.log('PASS unavailable authentication gateway fails closed');
  await context.close();
  console.log('ALL LOCAL AUTHENTICATION CHECKS PASSED');
}
run().catch(async error => { console.error(error); process.exitCode = 1;
  if (browser) {
    for (const context of browser.contexts()) for (const page of context.pages()) {
      try { write('failed-page.html', await page.content()); await page.screenshot({ path: path.join(fixtures, 'failed-page.png') }); } catch {}
    }
  }
  for (const name of containers) { try { write(`${name}.log`, pod('logs', name)); } catch {} }
}).finally(async () => {
  if (browser) await browser.close();
  if (backend) { backend.kill('SIGINT'); await Promise.race([new Promise(r => backend.once('exit', r)), new Promise(r => setTimeout(r, 5000))]); if (backend.exitCode === null) backend.kill('SIGTERM'); }
  for (const name of containers.reverse()) { try { pod('rm', '-f', '-v', name); } catch {} }
  smtp.close(); probe.close();
});
