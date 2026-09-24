# Remaining warm startup costs

Profiled master `80cd002` on 2026-09-24 after merging the small-panel SVG and
sequence-startup changes. No application code changed during this investigation.

## Setup

Same local machine, R and real-data launcher as the
[sequence startup comparison](2026-09-24-sequence-startup.md#method).
RiboCrypt and local ORFik loaded with `devtools::load_all()`;
`env -u LC_ALL R --vanilla -q -e` used to launch the app on localhost:7820.
Default human merged ATF4, ENST00000674920, RFP, translons enabled.
Fresh headless Chrome contexts at 1440 x 1000, warm app caches, empty browser
caches. This does not measure a returning user's warm browser cache or WAN load.

## Uninstrumented load times

Existing `warm-session-init.cjs`, six sequential contexts. Visit 0 primes the
app and is excluded from warm medians; visit 1 is retained despite being slower.
All sessions asserted the expected gene, transcript, 21 traces and no page errors.

| Visit | Visible (ms) | Visible and idle (ms) |
|---|---:|---:|
| 0, cold app | 4814.2 | 4832.7 |
| 1 | 2195.7 | 2246.4 |
| 2 | 1802.0 | 1851.2 |
| 3 | 1907.5 | 1939.1 |
| 4 | 1820.3 | 1862.8 |
| 5 | 1930.5 | 1978.9 |
| Warm median | 1907.5 | 1939.1 |

Median sequential intervals (independent medians need not sum to total):

| Interval | ms |
|---|---:|
| Navigation to HTML response end | 5.5 |
| HTML response end to Shiny connected | 588.4 |
| Shiny connected to plot value received | 330.9 |
| Plot value received to visible | 931.2 |
| Visible to observed idle | 48.4 |

The connection-to-value interval includes client/server exchanges, not just R CPU.

## Ranked candidates

1. **Earlier Plotly dependency loading.** Four separate instrumented visits
   showed plot-receipt-to-`newPlot` gaps of 701.0, 709.9, 800.6 and 388.9 ms.
   The last visit enabled CPU profiling and is not directly comparable.
   Plotly's main JS request starts late, around 0.99-1.08 s after navigation.
   It transfers 1,136,953 bytes and takes 126-132 ms locally. Installed
   htmlwidgets calls `Shiny.renderDependencies(data.deps)` before widget rendering;
   sampled browser stacks include synchronous jQuery loading/evaluation.
   Test earlier fetching first, then dependency registration if necessary.
   Download time is only part of the delay; do not promise removal of the whole
   0.7-0.8 s. Early execution could delay Shiny connection, so compare total load.
2. **Defer hidden-page UI construction/binding in the browser.** Server modules
   are already deferred, but all tab UIs are still included by `RiboCrypt_app()`.
   A separate DOM check found 69 initialized Selectize controls: 10 in Browser
   and 59 elsewhere, with about 3,100 DOM elements overall. The CPU-profiled
   session sampled about 157 ms under Selectize stacks, including descendants
   and input updates; not all of this is hidden-page or removable work.
   Keep this separate from dependency work. URL restoration, first-tab visits,
   returning-tab state and dependency initialization need regression coverage.
3. **Construct final sequence/axis state before `newPlot`.** The two remaining
   callback mutations cost about 74-83 ms synchronously in instrumented visits.
   Prepare the placeholder and axis style in R to avoid these redraws, retaining
   correct narrow-URL initialization and all browser variants. This is a smaller
   candidate; its savings are not established and cannot simply be added to the
   dependency estimate. Do not remove the required state changes themselves.

Initial `newPlot` synchronous work was 161-184 ms; its start-to-visible interval
was 269-309 ms including callback updates and browser scheduling. There was one
initial `newPlot`, not repeated whole-plot construction. The inspected R path
still caches controller, bottom panel, browser plot and rendered output.
Main gene, experiment and library updates already use server-side Selectize.
Further server micro-optimization is lower priority than the measured client work.

## Reproduction

Run the existing warm harness for six visits, then the Plotly profiler separately
with `VISITS=4`, against the same source-loaded app. Raw local outputs were
`/tmp/master-remaining-warm.json`, `/tmp/master-remaining-profile.log`,
`/tmp/plotly-render-profile.json` and `/tmp/plotly-render.cpuprofile`.
Profiler files are overwritten by later runs. No experimental trace rewriting
was enabled. The owned profiling server was stopped afterward.
