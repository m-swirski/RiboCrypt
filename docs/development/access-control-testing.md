# Access-control verification

The account integration is opt-in; production Keycloak/Google login must be
verified separately after the operator configures DNS, TLS, SMTP and credentials.
No production accounts or secrets were created during local testing.

## Automated checks

Run from the repository root using the source-loaded package:

```sh
env -u LC_ALL R --vanilla -q -e 'devtools::load_all("."); testthat::test_local(".")'
```

The new registry/session tests cover:

- Anonymous public access and no automatic private grants for new accounts.
- Separate issuers, users/workspaces and read/download capabilities.
- Persistence, immutable dataset mappings, SQL parameterization and revocation
  visible from a second database connection.
- Gateway-header spoofing, forged experiment paths, cross-identity session URLs
  and locked Shiny HTTP-method wrapping.
- Real ORFik experiment files, session-local environments and `bindCache`
  isolation with identical cache keys in two sessions.
- Run-filtered metadata/UMAP, authorized collection columns, path containment,
  denied shared reference artifacts and revoked cache retrieval.
- Collection normalization, summary coverage and clustering excluding library
  columns outside the selected experiment.
- Secret-redacted configuration printing and safe account redirect routes.
- FASTQ selected-file responses, sandbox/no-store headers, symlink containment
  and permission revocation without enabling download capability.

For the actual FASTQ selection helper and real Shiny HTTP routing, run:

```sh
NODE_PATH=/tmp/ribocrypt-megabrowser/node_modules node tests/manual/fastq-report.cjs
```

This disposable loopback fixture uses synthetic identities, starts a source-loaded
R app on port 7849, and removes its server/temporary files afterwards. It verifies
script-enabled sandbox rendering in Chrome, read-only access, cross-user/logout
denial, revoked grants, ignored file-query parameters and absence of globally
published report/sibling URLs. It does not replace production login tests.

## Real-browser harness

The developer-only files under `tests/manual` are excluded from the R build.
They use the local human merged experiment, collection and metadata launcher to
exercise the actual app with ATF4. The launcher defaults to
`~/Desktop/run_ribocrypt.R`; `RIBOCRYPT_RUN_SCRIPT` can select another launcher.
It stops parsing at the first `RiboCrypt_app()` call, then enables authentication
on that app. `RIBOCRYPT_REPO` can override the source checkout path.

```sh
env -u LC_ALL R --vanilla -q -f tests/manual/access-control-launch.R
# In a second terminal, with Playwright available through NODE_PATH:
node tests/manual/access-control.cjs
```

Requirements: Chrome at `/opt/google/chrome/chrome`, Playwright, SQLite CLI,
local `human_all_merged_l50` and `all_samples-Homo_sapiens`, plus the metadata
used by the launcher. The backend uses port 7837 and the test gateway 7838.
The harness writes synthetic experiments/database/screenshots beneath
`/tmp/ribocrypt-account`. It deletes/recreates that test database on startup.

**The gateway is not a login provider.** It trusts deliberately synthetic
cookie identities (`alice`, `bob`) on loopback to exercise app authorization.
Never deploy it or its fixed fixture secret, and never use real private data
with it. Private fixtures are renamed copies of already-public coverage, with
distinct synthetic Run IDs and metadata markers.

Checks include:

- Anonymous startup has only public choices; Alice and Bob each have only their
  own private experiment in addition to public choices.
- ATF4 loads automatically and each user can render the assigned private fixture.
- Actual CSV responses: public/anonymous 200, Alice read-only 403, Bob permitted
  200. Another identity cannot reuse the session's download URL.
- The actual server-side Samples DT exposes Alice's synthetic metadata only to
  Alice, Bob's only to Bob, and neither to anonymous visitors.
- Observatory renders the UMAP and all 3,857 public libraries for each identity.
- Revoking Bob's grant while the page is open clears/reloads the page and removes
  the private choice. The fixture grant is restored after the check.
- No JavaScript errors in the checked sessions; desktop screenshots are retained
  for inspecting the rendered plot and controls.

## Local results: 2026-10-07

The source-loaded full test suite passed 3,404 assertions with no failures,
warnings, errors or skipped tests. The Chrome harness passed for anonymous,
Alice and Bob sessions, including private fixture plots, metadata isolation,
CSV permissions, Observatory rendering and live grant revocation. These results
validate app authorization through the synthetic gateway, not production OIDC
registration, email verification or Google login.

## Deployment validation

### Real local OIDC stack

`tests/manual/auth-stack.cjs` exercises production-mode Keycloak backed by
PostgreSQL, OAuth2 Proxy and NGINX using the deployment templates. It generates
disposable secrets, a short-lived localhost TLS certificate and a loopback SMTP
sink. It never installs a system CA, changes system NGINX or sends external mail.
Rootless Podman and `fuse-overlayfs` are required; Docker/sudo are not needed.

Use dedicated container storage on a Linux filesystem with sufficient space.
NTFS copy-based storage is impractically slow. On this workstation the system
disk was almost full, so a dedicated `/dev/shm` directory was used for overlay
storage after checking available RAM; fixtures/logs were kept on the S drive.
Do not use RAM-backed storage on a memory-constrained production server.

Pull `quay.io/keycloak/keycloak:26.8.0`, `docker.io/library/postgres:17`,
`quay.io/oauth2-proxy/oauth2-proxy:v7.15.0` and `docker.io/library/nginx:1.28`
into the dedicated Podman store before running. These are local test image
versions, not a substitute for choosing maintained production pins.

```sh
env RIBOCRYPT_AUTH_TEST_ROOT=/media/roler/S/ribocrypt-auth-verification \
    RIBOCRYPT_AUTH_TEST_STORAGE=/dev/shm/ribocrypt-auth-storage \
    NODE_PATH=/tmp/ribocrypt-megabrowser/node_modules \
    node tests/manual/auth-stack.cjs
```

The script needs Playwright, `undici`, OpenSSL and the local RiboCrypt fixtures
described above. It reserves loopback ports 8443, 8444, 11025, 14180, 15432,
17837, 18180 and 7837. Do not run concurrently with the synthetic gateway
harness. It removes its named containers/volumes after success or failure;
test images and protected fixtures remain for inspection. Secrets, mail action
links and provider logs must never be committed.

Local verification on 2026-10-07 passed through the real provider/gateway:

- Production-mode Keycloak with PostgreSQL and the realm import's environment
  placeholders; NGINX configuration validation and localhost HTTPS.
- Registration, email verification, rejection of unverified authentication,
  minimum password length, PKCE login and immutable subject forwarding.
- Secure/HTTP-only/SameSite cookies, tampered-cookie rejection, failed forged
  callbacks and stripping browser-supplied identity/authorization headers.
- Provider and gateway logout; password reset with old-password rejection.
- PostgreSQL dump restored to a separate database; provider restart preserves
  the account, verified-email state, subject and working password.
- Actual RiboCrypt WebSocket startup and ATF4 rendering, isolated private
  catalog and private-fixture plotting, read-only CSV denial and anonymous
  reuse of session URLs denied.
- Logout in another tab clears/reloads an already-open authorized session.
- Gateway failure returns an error rather than silently granting access.

The full source-loaded R test suite was rerun: 3,404 passing assertions, zero
failures, warnings, errors or skips. A real cookie-key format pitfall found by
the harness is now documented: the proxy requires **URL-safe** base64 encoding
of 32 random bytes, not arbitrary ordinary base64. The realm template now
explicitly sets a 12-character minimum password policy.

Still unverified: Google brokering (needs Google credentials), real SMTP delivery
and domain reputation, publicly trusted production TLS/DNS, production firewall
and service-user permissions, MFA/rate-limit policy, and server-specific backup
operations. The loopback SMTP sink and self-signed test certificate do not
validate those deployment concerns. Host system packages/services, trust store,
firewall and existing application sessions were not modified.

The final harness passed on two consecutive full runs with independently
generated secrets; the last run also plotted the private fixture. Test containers
and RAM-backed images were removed after verification. Protected test artifacts
and the earlier NTFS image cache remain under the dedicated S-drive directory.

### Review fixes verified: 2026-10-07

After the security review, the source-loaded full suite passed **3,415 assertions**
with zero failures, errors, warnings or skipped tests. The dedicated FASTQ Chrome
harness passed on the final code: selected report/scripts render for a read-only
account, session URLs reject absent/different identities and revoked grants,
query parameters cannot select a sibling, and `/tmpuser` URLs return 404.
Unit coverage additionally rejects report symlinks outside the report directory.

The complete real OIDC harness also passed with the actual restricted login-host
template and NGINX 1.28: public admin/master routes (including an encoded admin
path) return 404, while realm discovery, registration, email verification, login,
logout and password recovery continue to work. An existing server-level security
header is retained on both anonymous and authenticated app responses. Database
backup/restore, provider restart, real ATF4/private-fixture plotting, read-only
download denial, cross-tab logout and unavailable-gateway denial also passed.

The cookie-header directive now lives at server scope. Operators must explicitly
carry forward any inherited `http`-level security headers on older NGINX versions;
this test cannot validate the production server's existing custom configuration.
Restart all existing Shiny workers to discard historical static report mappings.

### Production validation

`docker compose -f deploy/auth/compose.yaml config --quiet` validates the compose
template when its required environment variables are supplied. This does not
start Keycloak or validate the entire OIDC login flow. Host NGINX is not installed
on this development machine and Docker daemon access is unavailable; the local
OIDC harness uses rootless containers instead. Run
`nginx -t` and live registration/verification/Google/logout tests on the actual
server before enabling private production datasets.

See [deployment instructions](../deployment/authentication.md) and the
[saved roadmap](accounts-and-analysis-platform.md).
