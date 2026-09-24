# Earlier Plotly dependency loading

Date: 2026-09-24. Branch: `perf/early-plotly-loading`.
Control: master `80cd002`, including SVG small panels and optimized sequence startup.

## Decision

Keep normal early dependency registration. Warm median page-ready time fell from
**1.950 to 1.642 seconds**, an improvement of **308 ms / 15.8%**. Median time to
the visible plot fell from **1.922 to 1.594 seconds**, about **328 ms / 17.1%**.
A second experiment using a deferred main script was slower, at **1.800 seconds**
page-ready, and is not included in the implementation.

This is an incremental comparison against current master, not the original app.
These local measurements are not a guarantee for remote users or slower devices.

## Implementation

For apps configured with `plot_on_start=TRUE`, `RiboCrypt_app()` includes
`plotly_startup_dependencies()` in the initial UI. This helper obtains the normal
dependency list from `plotly::plot_ly()$dependencies`, without building or rendering
a dummy plot. Shiny/htmltools handle dependency registration, version resolution
and deduplication. There are no hard-coded bundle paths, custom JavaScript loaders,
async races, altered Plotly bundles or plot-data changes.

The full existing Plotly bundle is used, including WebGL support. Dynamic widget
rendering sees the dependencies as already registered and does not load them again.
Apps with `plot_on_start=FALSE` retain lazy loading. This condition uses the app
configuration, not individual URL `go` overrides; URLs still control whether a
plot is generated in that session. The additional `xml2` Suggests dependency is
used only by tests inspecting the generated HTML.

## Setup and method

Same hardware/software and real-data configuration as the
[current-master profile](2026-09-24-remaining-startup.md): i7-8850H, 6 cores /
12 threads, 31 GiB RAM, Ubuntu 24.04.2, R development r86148 (2024-03-18),
Shiny 1.11.1, Plotly R 4.12.0 / Plotly.js 2.25.2, htmltools 0.5.8.1.
RiboCrypt and local ORFik were source-loaded with `devtools::load_all()`;
R commands used `env -u LC_ALL R --vanilla -q -e`.

Three separate app processes used identical launcher/data settings: human merged
ATF4, ENST00000674920, RFP, translons enabled and human collection available.
Control on 7820 loaded master before code edits; candidate on 7821 loaded the
new helper. The experimental process on 7822 used the same helper but replaced
only the `plotly-main` dependency's script with
`list(src = original_script, defer = TRUE)` in memory before constructing the app.
No deferred-script override is shipped.

Each app received three primer visits. Then a single headless Chrome process
made seven rounds in control/early/deferred order, using a fresh empty-cache
1440 x 1000 browser context each time. Round 0 is reported but excluded as browser
warmup. No R tests or other browser jobs ran during the measured comparison.
The existing warm-session harness checked ATF4, transcript, 21 final traces and
absence of JavaScript/Shiny errors for every visit.

## Results

Navigation to visible plot and Shiny idle, milliseconds:

| Round | Control | Early normal script | Early deferred script |
|---|---:|---:|---:|
| 0, excluded | 2000.0 | 1660.0 | 1812.0 |
| 1 | 1955.0 | 1641.3 | 1824.6 |
| 2 | 1939.6 | 1657.3 | 1849.7 |
| 3 | 2061.1 | 1642.6 | 1802.7 |
| 4 | 1892.4 | 1701.3 | 1795.3 |
| 5 | 1945.1 | 1615.7 | 1796.8 |
| 6 | 1984.8 | 1610.7 | 1749.4 |
| Median, 1-6 | 1950.1 | 1642.0 | 1799.8 |

Median navigation-to-Shiny-connection increased from 612 to 1043 ms for the
retained candidate: earlier script execution does delay connection. However,
plot-receipt-to-visible decreased from 951 to 243 ms, more than offsetting it.
Connection-to-plot-receipt medians were 344 versus 303 ms, despite unchanged
server computation. Independent interval medians do not sum exactly to the
overall median, and small samples/fixed ordering limit causal precision.

## Regression checks

- Full R suite: 2,793 assertions passed, with no failures, warnings or skips.
- New unit tests check original dependency identity, disabled automatic plotting,
  deduplication and version registration in Shiny HTML, preserved jQuery selection,
  cached responses and query-string responses.
- New browser smoke checks confirm Plotly is ready at Shiny connection, exactly
  one main-bundle request/script exists, and WebGL sequence zoom works. Both fresh
  and populated browser HTTP caches passed. The fresh page fetched the bundle at
  10.5 ms after navigation and finished at 295.2 ms, before plot receipt at 1401.3 ms.
- Existing Chrome checks passed for codon hover, sequence copy, zoom/reset,
  viewport resize, initial URL zoom, Observatory settings round trip and
  highlighting, subset/reset actions and `go=FALSE`. No page errors were reported.

## Reproduction

After starting and priming the three source-loaded variants:

```sh
RESULTS=/tmp/early-plotly-threeway.json \
  NODE_PATH=/tmp/observatory-url-browser/node_modules \
  node benchmarks/warm-session-init.cjs \
  http://127.0.0.1:7820/,http://127.0.0.1:7821/,http://127.0.0.1:7822/ 7
NODE_PATH=/tmp/observatory-url-browser/node_modules \
  node benchmarks/plotly-startup-smoke.cjs http://127.0.0.1:7821/
env -u LC_ALL R --vanilla -q -e 'devtools::load_all("."); testthat::test_local()'
```

Use a fresh source-loaded app to test the change through normal `RiboCrypt_app()`.
Existing processes retain cached UI and will not gain the new dependencies.
