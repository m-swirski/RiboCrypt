# Observatory URL review and regression checks

Date: 2026-09-24. Branch: `perf/atf4-startup`.
Base commit: `8c1e8bb` (deferred page initialization).
This is a correctness review, not a new performance measurement.

## Findings and changes

- A delayed URL initialization event invalidated the Observatory module after
  selection restoration. Its reset observer replaced the saved cohort and label
  with All merged. Reset now follows explicit Generate Plot clicks or actual
  dataset changes, not every module invalidation.
- Restoring group 1 did not necessarily change the active ID, so selection
  restoration could be missed. Initialization now has an explicit completion
  signal. Restoring only group 2 also starts with an empty store rather than
  retaining an unrequested default group 1.
- Initial Plotly highlight messages could precede JavaScript handler
  registration. Plotly and DT now announce widget readiness and receive a
  targeted selection resync. This does not rebuild the full table.
- Get URL omitted frame subsets, frame colors, log scales, conservation,
  mapability, summary settings and translon settings. Capture, restoration and
  readiness now share one control registry. Observatory has no unique-alignment
  input, so that Browser-only setting is not included.
- Autostart previously waited for gene/transcript and any nonempty selection,
  but could run with default controls or the wrong cohort. It now waits for all
  saved controls and exact saved group membership, and fires once.
- Missing/invalid transcript defaults were resolved for input initialization
  but not the saved state used by autostart. Both now use the resolved defaults.
- A different collection in the URL now loads its annotation data before gene
  and transcript defaults are resolved, rather than starting with the default
  collection's annotation.
- URL payloads now validate version, scalar controls, group IDs and run indices
  before reaching reactive code. Legacy uncompressed base64 JSON remains
  supported alongside compressed base64url JSON. Unknown control fields are
  ignored. Invalid payloads fall back to normal startup.
- Shared links on other HTTPS hosts now retain their protocol, optional port
  and application path. Explicit hash arguments to `getPageFromURL()` are used.

## Verification

R tests use source code via `devtools::load_all(".")`, never an installed
RiboCrypt. Added coverage includes delayed input acknowledgements, active group
1, active group 2, group-2-only state, all-control update messages, false/empty
settings, cohort readiness, malformed payloads, legacy links, annotation
fallbacks, cross-collection initialization and HTTPS paths.
The final full suite passed 2,731 assertions with no failures, warnings or skips.

Real Chrome checks used a separately launched app at `127.0.0.1:7820`, the local
`human_all_merged_l50` Browser default, `all_samples-Homo_sapiens` collection,
and metadata from `~/Desktop/run_ribocrypt.R`. No user's server was replaced.
The test restored two groups containing `DRR277023`, with group 2 active and
named Saved group. DT showed one library, UMAP highlighted one point and Browse
rendered 19 traces without a Go click. A URL generated through the real Get URL
button restored the cohort, labels, active group and every exported Browse
control in a second session. Remove subset restored the full table; Subset to
page reduced it to 15 libraries. Clicking a row did not subset immediately;
Subset to selected then reduced the table to that one library. A `go = FALSE`
link did not render a Browse plot. The final Chrome run passed all assertions,
captured a screenshot, and reported no page errors.

The browser checks are reproducible with `observatory-url-smoke.cjs`:

```sh
npm install --prefix /tmp/observatory-url-browser playwright@1.51.1
NODE_PATH=/tmp/observatory-url-browser/node_modules \
  node benchmarks/observatory-url-smoke.cjs http://127.0.0.1:7820/
env -u LC_ALL R --vanilla -q -e 'devtools::load_all("."); devtools::test()'
```

Start the app from this checkout using `devtools::load_all(".")`, with the real
human data above. The script defaults to `/opt/google/chrome/chrome`; set
`CHROME` for another Chrome executable. The screenshot is written to
`/tmp/observatory-url-smoke.png`. Timing depends on local data/cache state;
the smoke test deliberately does not assert a startup latency threshold.

## Limits and remaining design boundaries

- A yeast URL resolved to `YAL069W` / `YAL069W_mRNA`, but its plot could not be
  verified: the local yeast experiment points to unavailable collection FST
  tables. The cross-collection initializer has a focused unit test; production
  yeast plotting still needs a check on a server with those tables.
- Shared URLs represent materialized cohorts, not saved DT search expressions,
  pagination, clicked-row editing state or UMAP camera position. Get URL writes
  the same cohort into both compact selection maps. Existing codec support for
  separate plot/table maps is retained, but differing maps without saved filter
  expressions are not covered by the end-to-end checks.
- State is restored at session/module initialization. Pasting another state
  into the fragment of an already initialized Observatory is not a live import
  mechanism; open the link in a new page or reload.
- Missing/renamed runs or collections are not migrated by these changes.
  Autostart requires the requested cohort, rather than silently plotting a
  different one. Invalid payloads currently fall back without a detailed
  user-facing validation message.
- Full non-default phyloP/mapability data access is not exercised by this smoke
  test; their saved controls are covered by restoration tests.

No master commit, merge, push or production deployment is part of this change.
