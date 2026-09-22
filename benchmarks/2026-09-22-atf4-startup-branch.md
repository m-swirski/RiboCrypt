# ATF4 startup: deferred server initialization

Date: 2026-09-22. Branch: `perf/atf4-startup`.
Base: `fdea24ce4972382928ee8f8530c236e783e07d1e`.

## Plan and implementation

The [original profile](2026-09-22-atf4-startup-local.md) showed that warm users
already reuse the ATF4 plot cache. The first candidates were unnecessary hidden
server initialization and repeated work in initial gene selection.

1. **Defer hidden server groups until first entry.** Browser starts immediately.
   MegaBrowser and Observatory start together on first entry to either page,
   preserving their shared collection state. Analysis, Metadata and tutorial
   start on first entry to their respective tabs. Each group starts once per
   session; leaving a tab does not destroy its state. The existing navbar/URL
   router controls activation, including direct hash links.
2. **Remove redundant gene-label deduplication.** Gene existence and default
   selection checks use the label column directly. Duplicated transcript labels
   do not affect membership or the first fallback gene. Dropdown choice creation
   still deduplicates labels as before.
3. **Leave static UI serialization for a separate pass.** Its approximately
   0.6-second warm HTTP response remains a candidate, but pre-rendering must
   preserve Shiny dependency and singleton registration. This branch keeps that
   behavior and the initial markup unchanged.

The changes are in `R/app_startup.R`, `R/RiboCrypt_app.R` and
`R/Select_panels_update.R`. Experiment reads, plot calculations, Plotly traces,
DT configuration and plot cache keys were not changed.

## Environment and method

Same hardware, R version, ORFik checkout and data configuration as the
[local baseline](2026-09-22-atf4-startup-local.md#hardware-and-software-baseline):
Intel i7-8850H, 6 cores/12 threads, 31.03 GiB OS-visible RAM, Ubuntu 24.04.2,
R development revision r86148, Chrome 146 headless, Playwright 1.51.1.

Both runs used `devtools::load_all()` and the same ATF4 experiment, transcript,
RFP selection, metadata and enabled collection pages. Each used a fresh R process,
followed by four sequential browser contexts with empty browser caches and a
1440 x 1000 viewport. The first context primed the app cache; the remaining
three measured warm-user loads. Visible-plot and Shiny-idle criteria match the
original report. Function tracing and 10 ms Rprof sampling remained enabled.

The final branch benchmark ran on port 7820. A fresh control benchmark ran
afterward on port 7821 from a detached worktree at the base commit, without
switching or modifying `master`. The two benchmark page-load sequences ran
sequentially; another app process could remain idle. A smoke-test page was still
waiting for a table action when the control sequence began. First-user times
exclude R/package loading and construction of the app object.

## Results

| Session | Fresh master: plot visible | Branch: plot visible | Fresh master: idle | Branch: idle |
|---|---:|---:|---:|---:|
| First user | 5.864 s | 5.695 s | 5.986 s | 5.882 s |
| Second user | 4.268 s | 3.054 s | 4.556 s | 3.281 s |
| Third user | 3.346 s | 2.857 s | 3.632 s | 3.080 s |
| Fourth user | 3.373 s | 2.694 s | 3.628 s | 2.917 s |
| Warm median | 3.373 s | 2.857 s | 3.632 s | 3.080 s |

- Warm median completion improved by **0.552 s (15.2%)** versus the fresh control.
- Warm median visible plot improved by **0.516 s (15.3%)**.
- The second user specifically improved by **1.275 s (28.0%)** to idle.
- Versus the original 3.841-second warm median, completion improved by
  **0.761 s (19.8%)**. The fresh control is the preferred comparison.
- First-user completion improved by only 0.104 s in these runs. Do not claim
  a substantial cold-start improvement from this sample.

The initial, intermediate branch experiment measured a 2.926-second warm median.
That was before the final gene validation change and is not the final result;
the difference also illustrates normal run-to-run variation.

## Where time changed

Fourth-user sequential intervals:

| Stage | Fresh master | Branch |
|---|---:|---:|
| Navigation to completed HTML response | 0.599 s | 0.594 s |
| HTML response to Shiny connection | 0.597 s | 0.572 s |
| Shiny connection to received plot value | 1.057 s | 0.341 s |
| Received plot value to visible plot | 1.120 s | 1.188 s |
| Visible plot to observed idle | 0.255 s | 0.223 s |
| Total | 3.628 s | 2.917 s |

Wall-time traces inside the server interval:

| Server registration | Fresh master | Branch |
|---|---:|---:|
| First user | 0.935 s | 0.209 s |
| Second user | 1.412 s | 0.470 s |
| Third user | 0.657 s | 0.105 s |
| Fourth user | 0.579 s | 0.119 s |

The fourth-user registration reduction was **79.4%**. Traces confirm the hidden
server groups did not initialize during the branch's default Browser loads.
Warm sessions still skipped the controller, sequence-panel and coverage-panel
calculations through the existing app cache. Browser-side rendering remains
approximately 1.2 seconds and was not optimized here.

App-object construction before accepting users was 2.924 s for the fresh control
and 2.988 s for the branch. There is no meaningful pre-server improvement in
these measurements.

## Verification

- Full suite: `env -u LC_ALL R --vanilla -q -e 'devtools::load_all("."); devtools::test()'`.
  **2,671 assertions passed, zero failures, warnings or skips.**
- Added tests cover deferred activation, one-time initialization, independent
  sessions, module reactivity/state after leaving and returning, and gene
  selection with duplicated labels, missing values and empty input tables.
- Chrome exercised Observatory page subset, single-row subset, removal back to
  3,857 libraries and subset preservation when leaving and returning.
- Chrome opened MegaBrowser, searched the server-side Codon gene selector,
  opened Samples and tutorial, and returned to the original Browser plot.
- Both full and gene/transcript-omitting Browser URLs automatically plotted ATF4.
- Direct hash entry worked for Observatory, MegaBrowser, Codon dwell time,
  Studies and Predicted Translons. The latter verified selector initialization;
  its explicit Search action was not part of this check.
- All four benchmark sessions displayed the same ATF4 transcript and 21 traces.
  The completed browser smoke suite recorded no JavaScript or non-silent Shiny errors.

### Existing Observatory URL restoration issue

An additional compressed `obs_state` URL test did not restore a subset on either
the branch or the unchanged base commit. The test specified human libraries,
`view = "umap"`, tissue/cell-line coloring, two groups with run `DRR277023`, and
active group `2` labeled `URL restored group`. In both versions group 2 became
active, but its label reset to `All merged` and the DT displayed all 3,857
libraries instead of the saved single-library subset. The same failure occurred
with only group 2 in the saved state.

This is an observed existing issue, not a passing URL-restoration test or a
regression introduced by deferred initialization. It remains unfixed in this
performance branch. The separate check used a zlib-compressed, base64url-encoded
JSON state matching the current URL decoder. Logs and screenshots are in
`url-state.log` and `url-state-7820.png` / `url-state-7821.png` in the local
profiling directories. Investigate this separately before demonstrating saved
Observatory subset links.

## Tradeoffs and remaining work

First entry to a deferred group now pays its initialization cost. This is
deferred work, not removal of those features or a promise that every tab is faster.
Related pages still initialize as groups to preserve their existing shared state.
Their static UI remains in the initial page, so this change does not shrink the
initial asset payload.

There are only three warm observations per run, not a production percentile or
concurrency benchmark. Local network, tracing, browser polling and system load
affect the measurements. The control followed the branch rather than using a
randomized order. Gene-validation savings were not isolated from lazy startup
in the final comparison. Full analyses, exports and external services were not
exercised by the browser smoke test.

Local evidence is in `/tmp/ribocrypt-branch-profile/` and
`/tmp/ribocrypt-master-recheck/`: browser JSON, screenshots, function timing logs,
profiling scripts and Rprof data. Those temporary files are not part of the
repository. The tables above preserve the results for future comparisons.
