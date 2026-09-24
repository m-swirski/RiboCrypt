# Warm ATF4 startup: cache the static UI response

Date: 2026-09-24. Branch: `perf/atf4-startup`.
Control commit: `e214846` (including the Observatory URL fixes).

## Decision

Keep this iteration. In an alternating comparison, warm median completion fell
from **2.700 s to 2.182 s**, a **0.518 s / 19.2% improvement**. The visible ATF4
plot arrived 0.503 s sooner. This is an incremental improvement over the earlier
deferred-server work, not a comparison against master.

No plot calculations, scientific settings, cache keys, dropdown data or hidden
page behavior were changed. Further browser-side Plotly work is deferred; its
cost is now larger than the remaining server interval and deserves a separate
profile before changing rendering behavior.

## Profile and implementation

Fresh measurements confirmed repeated static HTML rendering was worth targeting.
A direct warm invocation of the original Shiny UI HTTP handler took 0.373 s.
Three subsequent profiled calls sampled 1.30 s in total, with 95.4% under
`renderTags()`, including tag serialization and dependency resolution. Rprof
used an effective 10 ms interval; the requested 1 ms interval was not supported
on this machine. Sampling totals are not additive function wall times.

`cache_static_app_ui()` wraps the app's existing HTTP handler. It retains the
first successful plain `GET /` response and reuses it for later plain requests.
It keeps the complete Shiny response, including headers, singleton registration
and dependency declarations. It does not rebuild HTML through private rendering
functions or cache any session's reactive state, input values or plots.

Query-string requests, non-root paths and non-GET methods bypass this cache.
This preserves the original path for gene/library links, bookmarking and
showcase query parameters. Fragment-only links share the static document but
still restore their own state after connection, as before.

The cache is local to the app object and fills only after the first ordinary
request. It does not improve process startup or first-user HTML rendering.
It retains one response of roughly the original 190 kB decoded HTML size;
the initial payload and client asset requirements are unchanged.

**Constraint:** this wrapper is appropriate only because RiboCrypt constructs a
request-independent static UI. If UI construction becomes user-, request- or
runtime-configuration-dependent, remove or redesign the wrapper. Shiny's normal
session HTTP routes remain outside the cached root request. App-object handler
compatibility and full response identity are covered by tests.

## Environment and method

Same computer/data configuration as the
[hardware baseline](2026-09-22-atf4-startup-local.md#hardware-and-software-baseline):
i7-8850H, 6 cores/12 threads, 31 GiB OS-visible RAM, Ubuntu 24.04.2.
R remains development revision r86148 (2024-03-18); Shiny 1.11.1,
htmltools 0.5.8.1, Plotly R 4.12.0. Benchmark app processes loaded both RiboCrypt
and local ORFik source with `devtools::load_all()`. ORFik source commit remains
`e43cf261bb85a5698a71e8349a1974bd3fbc3d44`.

The control server was loaded before editing package code, on port 7820. A
separate R process loaded the candidate on 7821. Neither replaced a user's
server or changed master. Both used the launcher's experiment/metadata setup,
`human_all_merged_l50`, ATF4, RFP, transcript ENST00000674920, translons enabled,
and the human collection available for deferred pages.

Each server first received six visits in fresh Chrome contexts to warm app
caches and R compilation. Then a fresh Chrome process alternated candidate and
control visits, using a new empty-cache browser context for every page, at
1440 x 1000. Six rounds were collected; round 0 was treated as browser warmup
for both servers. The remaining five rounds are the primary comparison. The
candidate's round-0 completion was an unusually slow 5.816 s (control 2.696 s),
mostly between plot receipt and visibility; it is recorded here rather than
hidden. This is not a production-tail latency benchmark.

Completion requires `browser-c._fullData`, a main SVG, then clearing Shiny's
busy class. Every trial asserted ATF4, the expected transcript, 21 traces and no
page/output errors. The measurement harness polls at 25 ms, captures navigation
response end, Shiny connection and plot-value receipt, and waits an unmeasured
second before closing each context. Screenshots confirmed coverage and sequence
tracks. Browser timing instrumentation differs slightly from the September 22
report, so use this run's paired control for the claimed improvement.

No Rprof or R tests ran during the final alternating comparison. Other system
activity, CPU frequency and thermal state were not controlled. Visits were
sequential, not concurrent, and candidate-first ordering was not randomized.

## Results

All values below are seconds from navigation start.

| Warm round | Control visible | Candidate visible | Control complete | Candidate complete |
|---|---:|---:|---:|---:|
| 1 | 2.574 | 2.203 | 2.622 | 2.229 |
| 2 | 2.742 | 2.214 | 2.773 | 2.266 |
| 3 | 2.652 | 2.140 | 2.700 | 2.182 |
| 4 | 2.770 | 2.150 | 2.793 | 2.167 |
| 5 | 2.645 | 2.121 | 2.690 | 2.167 |
| Median | **2.652** | **2.150** | **2.700** | **2.182** |

| Median milestone | Control | Candidate |
|---|---:|---:|
| HTML response complete | 0.545 | 0.005 |
| Shiny connected | 1.154 | 0.621 |
| ATF4 plot value received | 1.440 | 0.938 |
| Plot visible | 2.652 | 2.150 |
| Visible and Shiny idle | 2.700 | 2.182 |

Milestone medians are summarized independently, not additive stage samples.
Almost all of the end-to-end gain is explained by avoiding repeated HTML
rendering. The exploratory six-visit runs gave a similar warm completion change
(2.871 s to 2.311 s). Their first users completed at 5.634 s and 5.690 s:
there is no evidence of a first-user improvement.

## Verification and reproduction

- Full R suite: **2,755 assertions passed**, no failures, warnings or skips.
- Added tests compare the complete cached response to Shiny's original response,
  including singleton and dependency registration. They cover query/method/path
  bypass, per-app isolation and retry after unsuccessful responses.
- Real Chrome Observatory checks passed for saved cohorts, DT/UMAP highlighting,
  all exported settings, Get URL round trip, subset/reset actions and `go=FALSE`.
- ATF4 startup was verified on both servers with screenshots and 21-trace checks.
- Both a full Browser gene/transcript URL and a short `dff`, `library`, `go=true`
  URL automatically plotted ATF4 on the candidate. Their query requests still
  passed through the uncached handler.

```sh
env -u LC_ALL R --vanilla -q -e 'devtools::load_all("."); devtools::test()'
NODE_PATH=/tmp/observatory-url-browser/node_modules RESULTS=/tmp/warm-alternating.json \
  node benchmarks/warm-session-init.cjs http://127.0.0.1:7821/,http://127.0.0.1:7820/ 6
NODE_PATH=/tmp/observatory-url-browser/node_modules \
  node benchmarks/observatory-url-smoke.cjs http://127.0.0.1:7821/
```

The scripts use Playwright 1.51.1 and `/opt/google/chrome/chrome` (override with
`CHROME`). Start the two app versions using the same source-loaded/data settings
documented in the earlier baseline, and warm both before an alternating run.
`RESULTS` saves individual timings as JSON; screenshots go to
`/tmp/warm-init-<port>.png`. Local exploratory results and profiling logs remain
under `/tmp/warm-*` and are not required to read the retained tables above.
