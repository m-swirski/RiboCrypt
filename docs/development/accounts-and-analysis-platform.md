# Accounts and analysis platform

## Deployment decision

RiboCrypt runs behind NGINX on a personal server. Use Keycloak for ordinary
registration and optional Google identity brokering, and OAuth2 Proxy for OIDC
sessions at the gateway. Do not implement passwords or accept Google tokens in
Shiny. Account keys are `(issuer, subject)`, never an email address.

Keep the biological service independent of this identity provider. Anonymous
visitors retain public access; registration alone grants no private datasets.
Administrators assign workspace membership and dataset permissions.

## Roadmap

1. Identity and administrator-assigned datasets. Persist accounts, workspaces,
   membership and dataset grants. Resolve authorization on the server, covering
   experiment loading, forged inputs, URLs, metadata, downloads and caches.
   Require isolation tests before enabling private data in Shiny.
2. Saved projects: genes, sample selections, regions, analysis parameters and
   permission-aware sharing. Record dataset versions for reproducibility.
3. A common analysis service for Shiny and AI clients, with bounded background
   workers separate from interactive sessions. Expose a versioned API; MCP can
   provide an additional interface, not a second authorization system.
4. Paid AI access using service accounts, scoped capabilities and budget caps.
   Plan, validate, quote, reserve, execute, deliver and settle each job. Support
   idempotency, cancellation and a provider-independent usage ledger. Start with
   prepaid credits rather than unlimited autonomous jobs.

Workspaces are the sharing boundary. Preserve existing ORFik/FST storage; a
small relational database holds identity and grants, later jobs/accounting.
Public caches may be shared. Private caches initially belong to a session;
workspace caches later require authorization-safe keys and dataset versions.
Biological results should expose study-aware evidence, uncertainty and an
inspectable provenance manifest rather than unsupported conclusions.

## Implementation status

Implemented:

- Optional SQLite registry with immutable dataset IDs and canonical experiment
  directory/name mappings; explicit public visibility, deny by default.
- Administrator functions for workspaces, memberships and read/download grants.
- Gateway identity verification using a deployment secret and fixed OIDC issuer.
- An account record (issuer/subject, first/last seen) is created only after a
  gateway-authenticated session; passwords and email verification remain in Keycloak.
- Authorization queries and an explicit resolver that checks current grants on
  every call. Revocation takes effect on the next resolver call.
- Tests for anonymous access, separate users/workspaces, overlapping experiment
  names, spoofed identity, revocation, persistence and capability separation.

`RiboCrypt_app(access_control = ribocrypt_access_control(...))` now enables
gateway-authenticated sessions. Without this argument, legacy public behavior
is unchanged and all data passed to that app must still be public.

Authenticated operation:

- The startup HTML contains no experiment catalog. The authorized app is built
  in the verified WebSocket session and its server starts after controls bind.
- The registry, not global ORFik configuration or `all_exp`, is the source of
  dataset visibility and canonical experiment paths. Authorized experiment names
  must be unique within a session. Two users may have same-named experiments in
  different directories, provided neither has both versions assigned.
- All package experiment reads use the guarded resolver. Legacy URLs retain
  experiment names but can address only the session's authorized catalog.
- Metadata and UMAP rows are filtered by experiment Run IDs. Collection coverage
  is restricted to authorized library columns before statistics/normalization.
- Reference-wide translon/protein artifacts cannot reliably be filtered by Run.
  They are denied by default; administrators must explicitly approve
  `reference_annotations = TRUE` when registering that dataset. Approval means
  all reference-wide annotations are suitable for that dataset's readers.
- Experiments and plot caches are session-local. Cache retrieval/storage and
  session HTTP endpoints recheck current grants; another identity cannot use a
  session's DT/download URL. Revocation closes the session on the next poll
  (default one second), clears its cache and reloads the browser.
- Server exports require download permission for **every dataset in the current
  session catalog**, conservatively preventing mixed-public/private exports.
  This is not DRM: users can retain data already displayed. An administrator
  should not grant viewing if copying displayed data is unacceptable.

Production activation still requires DNS/TLS, Keycloak/SMTP/Google configuration,
secrets, NGINX integration and live identity-provider acceptance testing. No
production services or credentials are created by loading this package.

## Administrator example

Run with the source-loaded package; install DBI and RSQLite for this optional
feature. Keep the database outside the repository, with restricted filesystem
permissions and backups. Public registration is explicit, never inferred from
the experiments installed on disk.

```r
devtools::load_all(".")
db <- ribocrypt_access_db("/srv/ribocrypt/access.sqlite")
ribocrypt_access_workspace(db, "lab-a", "Lab A")
ribocrypt_access_member(db, "https://login.example.org/realms/ribocrypt",
                       "oidc-subject-from-provider", "lab-a")
ribocrypt_access_dataset(db, "custom-study-v1", "my_experiment",
                        "/srv/ribocrypt/private/experiments")
ribocrypt_access_grant(db, "lab-a", "custom-study-v1", download = TRUE)
DBI::dbDisconnect(db)
```

These functions are administrator APIs, not unauthenticated web endpoints.
Protect database write access; only trusted administrators may change grants.
See [gateway deployment](../deployment/authentication.md).
