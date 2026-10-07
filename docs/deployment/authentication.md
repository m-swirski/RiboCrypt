# Authentication deployment

This design targets the existing personal NGINX server. Authenticated Shiny
routing is opt-in; public operation remains unchanged. Deployment templates are
in [deploy/auth](../../deploy/auth). Adapt them to the existing server rather
than replacing its working NGINX configuration blindly.

## Identity provider

Run Keycloak in production mode with PostgreSQL, persistent volumes, backups,
a fixed HTTPS hostname and restricted administrative access. Do not use
`start-dev` in production. Create a `ribocrypt` realm and enable self-registration,
email verification and password reset with a configured SMTP service. Google
login is an optional identity broker in the same realm; application identity is
still the Keycloak issuer/subject pair. Do not automatically link accounts based
on an unverified email address.

Create a confidential OIDC client for OAuth2 Proxy with exact callback URLs,
for example `https://ribocrypt.org/oauth2/callback`. Avoid wildcard redirects.
Store client credentials outside Git. Use the provider's account interface for
password changes and recovery; no password enters RiboCrypt or its database.

Use a maintained, pinned OAuth2 Proxy release. Configure generic OIDC discovery,
the realm issuer, secure HTTP-only cookies and `set_xauthrequest`. Confirm the
user claim sent to NGINX is the immutable `sub`, not an email or display name.
The generic OIDC provider in OAuth2 Proxy 7.15 uses `sub` for its User field;
verify this against the pinned release before granting data. The supplied
configuration uses PKCE, nonce checking and verified email, without insecure
issuer/email bypasses.

## NGINX boundary

1. Keep Shiny, OAuth2 Proxy and Keycloak backend ports inaccessible externally.
   Bind local services to loopback or use private networking with firewall rules.
2. Route `/oauth2/` to OAuth2 Proxy; use its auth endpoint in an internal
   `auth_request` location. A successful auth response supplies the subject;
   a 401 means anonymous public access, not an authenticated private request.
   Backend failures must not be mistaken for successful authentication.
3. Overwrite `RiboCrypt-Subject` from the auth subrequest result, or clear it for
   anonymous users. Never forward a browser-supplied identity header.
4. Overwrite `RiboCrypt-Gateway-Secret` with an independently generated secret
   (at least 32 random bytes, encoded as hex). Inject from a protected NGINX
   include and supply the same value through the Shiny service environment.
   Never log or send this secret to the browser. Clear bearer/token headers not
   needed by the application. Restrict configuration-file and environment access.
5. Apply the same authentication to WebSocket handshakes and download requests,
   not only the initial HTML response. Preserve standard Shiny upgrade handling.
6. Configure logout to clear both the proxy cookie and the provider session,
   with an exact allowlisted return URL. Login redirects must also be restricted
   to this site. Add gateway rate limits; require MFA for administrators.

Use the fixed issuer from server configuration, never an issuer supplied by a
request. Request headers map to `HTTP_RIBOCRYPT_SUBJECT` and
`HTTP_RIBOCRYPT_GATEWAY_SECRET` in httpuv. `ribocrypt_access_identity()` rejects
missing/incorrect gateway credentials even for anonymous requests. It does not
validate tokens itself and cannot compensate for an exposed backend or proxy
configuration that forwards forged headers.

## Database operations

Install DBI and RSQLite only if using this feature. Place the registry outside
the package and web roots, in a protected directory (0700) and database file
(0600). Directory protection also covers SQLite journal files. Database writes
are administrator operations, never public HTTP handlers. New accounts need
explicit workspace membership and grants; registration alone adds no access.
Public datasets permit both reading and downloading; private grants may permit
reading without export. IDs identify immutable dataset versions. To replace
storage, register a new version instead of silently retargeting an existing ID.

Accounts are recorded after trusted gateway authentication, by issuer/subject
with first/last-seen timestamps. No password is stored. Google account linkage
must be configured in Keycloak, not by changing registry identities. Audit
trails and comprehensive account deletion remain future administrative work.

## Install and enable

1. Choose maintained, pinned image versions in a protected `.env` based on
   `deploy/auth/.env.example`. Use PostgreSQL 17 for the supplied volume path;
   other major versions can require a different data mount. Generate independent
   database, administrator, client, cookie and gateway secrets. Protect the
   directory and `.env`; do not paste real secrets into chat or commit them.
   The cookie secret must encode exactly 32 random bytes as URL-safe base64:
   `openssl rand -base64 32 | tr '+/' '-_' | tr -d '\n'`. Ordinary base64
   containing `+` or `/` is not accepted by the pinned proxy. Gateway/client
   secrets are independent and do not use this cookie-key format.
2. Create `login.ribocrypt.org` (or another chosen login hostname), provision its
   TLS certificate and adapt `nginx-login.conf.example` to that server block.
   The public proxy allows only `/realms/ribocrypt/` and `/resources/`; adjust
   the realm allowlist if you rename it. `/admin/` and the master realm are not
   public. Use a separate VPN/SSH-accessible private administrative ingress,
   configured with Keycloak's administrative hostname and trusted proxy headers.
   Backend ports remain loopback-only. Configure SMTP in the imported realm
   before enabling registration: email verification is required by default.
   The realm template sets a 12-character minimum password policy. With the
   tested Keycloak version, registration verifies email before prompting for a
   password. For an existing realm, explicitly apply policy changes in Keycloak;
   restarting with an updated import file will not overwrite the realm.
3. Start PostgreSQL/Keycloak with `docker compose up -d postgres keycloak` from
   `deploy/auth`, then configure/verify the realm in Keycloak. The initial realm
   import supports environment-variable placeholders and is skipped on later
   starts if the realm already exists. Changing the JSON will not overwrite live
   configuration. Configure Google brokering separately using Google's console
   and Keycloak's exact broker callback URL.
4. Start OAuth2 Proxy after the HTTPS issuer is reachable. Adapt the NGINX app
   and proxy snippets, add the upgrade map in the `http` context, and install the
   protected gateway-secret include. Run `nginx -t` before a deliberate reload.
   Keep the refresh-cookie `add_header` at server scope alongside the site's
   security headers, not inside a location. On NGINX before 1.29.3, adding any
   header at a scope suppresses all inherited headers: explicitly repeat any
   existing `http`-level security headers at server scope. Verify actual response
   headers for anonymous and logged-in requests, not just configuration syntax.
5. Create the access registry outside the repository and explicitly register
   public datasets as well as private datasets/grants. Public collections must
   contain only public libraries; visibility is at experiment granularity, not
   inferred from the fact that files share a reference directory. Run IDs must
   be present for collection/metadata filtering. Grant shared annotation access
   only after checking the reference-wide files contain no unauthorized evidence.
6. Source-load the package and add the configuration to the existing launcher:

```r
devtools::load_all(".")
access <- ribocrypt_access_control(
  database = "/srv/ribocrypt/access.sqlite",
  issuer = "https://login.ribocrypt.org/realms/ribocrypt",
  gateway_secret = Sys.getenv("RIBOCRYPT_GATEWAY_SECRET")
)
app <- RiboCrypt_app(metadata = metadata, browser_options = browser_options,
                    access_control = access)
shiny::runApp(app, host = "127.0.0.1", port = 3838)
```

`all_exp` is not used to grant access in authenticated mode: the registry is
authoritative. An explicit `all_exp_meta` can classify which authorized datasets
are collections; otherwise names starting with `all_samples-` are used. New
grants appear after reloading. Read-only private grants conservatively disable
server downloads for the whole session. Viewing itself necessarily reveals data
that a user can retain, irrespective of export controls.

The session HTTP guard wraps Shiny's `handleRequest` method because Shiny routes
DT/download endpoints before the app handler. Browser testing must be repeated
when upgrading Shiny. Backend session requests also validate the gateway secret
and exact identity, including after logout or switching accounts.

FASTQ HTML reports use session-scoped data-object endpoints, not global static
directories. Each request rechecks read permission; only the selected file is
served, with `no-store` and an opaque-origin script-enabled sandbox. Sibling files
and symlinks outside the report directory are not published. Read-only grants
can view reports without enabling server downloads. Restart existing app workers
when deploying this change: an old worker may retain its prior `tmpuser` mapping.

NGINX header placement follows its documented
[inheritance rules](https://nginx.org/en/docs/http/ngx_http_headers_module.html#add_header).

Use a short proxy session lifetime/refresh policy. Shiny authenticates the
WebSocket handshake, not each incoming frame; an existing socket can outlive a
provider cookie. The browser checks an identity-bound session HTTP endpoint every
30 seconds and when a tab becomes visible, so expired cookies/account switches
clear and reload that page. Registry revocations are checked each second and
every guarded read. Production logout must clear both provider/proxy state and
reload/close the originating page, and should be tested with multiple already-open
tabs. Background-tab timer throttling and offline browsers can delay checks;
data already displayed cannot be recalled.

## Acceptance before rollout

- Existing anonymous browsing and URLs still work, without private catalog names.
- Two different users see only their authorized libraries and metadata.
- Forged dataset IDs, experiment names, paths and identity headers fail closed.
- Caches, reference-wide FST files and downloads do not leak another workspace.
- Revoking membership or a grant clears/terminates active private sessions.
- Restarting services preserves grants, and backup/restore is exercised.
- Registration, Google login, recovery and logout work in real browsers.

Authoritative deployment references:
[Keycloak containers](https://www.keycloak.org/server/containers),
[Keycloak reverse proxy](https://www.keycloak.org/server/reverseproxy),
[OAuth2 Proxy configuration](https://oauth2-proxy.github.io/oauth2-proxy/configuration/overview/).
