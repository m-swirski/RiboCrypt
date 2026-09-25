# Observatory Select libraries startup

Measured 2026-09-24 against the user's existing server at
`http://127.0.0.1:4723/#Observatory`. No application code was changed, no R
profiling hooks were installed in that process, and the server was not restarted.
Checkout at investigation time: `7166831`; the running process's exact source
revision cannot be established from browser measurements.

## Setup

Same computer as the [hardware baseline](2026-09-22-atf4-startup-local.md#hardware-and-software-baseline).
Six independent empty-cache Chrome contexts, 1440 x 1000, sequential requests
to an already running app. Visit 0 is browser warmup, not a cold-server measurement.
Human collection `all_samples-Homo_sapiens`; 3,857 libraries/UMAP points,
158 UMAP traces, default DT page size 15, all libraries displayed, no subset.
Every visit checked the collection, populated table, visible UMAP and absence
of JavaScript/Shiny errors. A screenshot confirmed both rendered correctly.

The current local launcher calls `devtools::load_all()` and specifies
`plot_on_start=FALSE`, with `init_tab_focus="MegaBrowser"`; the URL opens the
Observatory Select libraries view. These launcher settings are context, not a
fingerprint of the already running R session. Unlike the earlier ATF4 benchmark,
this is not a newly launched, controlled app process with automatic Browser plotting.

## Timings

Milliseconds from navigation. DT/UMAP values are Shiny output arrival times;
DT draw is its first `draw.dt` event. Both-ready is the first polling observation
of populated visible widgets, not an exact paint timestamp. Ready also requires
Shiny's busy class to clear. Polling interval is 25 ms and cannot run while the
main thread is blocked. Each page stays open another 1.5 seconds for observation.

| Visit | Connected | DT value | UMAP value | DT draw | Both ready | Ready and idle |
|---|---:|---:|---:|---:|---:|---:|
| 0, excluded | 698.9 | 1880.3 | 2025.6 | 3469.9 | 3621.3 | 3708.6 |
| 1 | 615.6 | 1580.1 | 1722.1 | 3197.4 | 3420.3 | 3420.3 |
| 2 | 620.9 | 1636.2 | 1777.9 | 3250.7 | 3464.5 | 3464.5 |
| 3 | 645.3 | 1793.5 | 1942.3 | 3617.4 | 3855.5 | 3855.5 |
| 4 | 625.3 | 1666.1 | 1814.1 | 3384.5 | 3543.2 | 3645.2 |
| 5 | 649.4 | 1667.7 | 1813.0 | 3369.7 | 3531.5 | 3636.5 |
| Median, 1-5 | 625.3 | 1666.1 | 1813.0 | 3369.7 | 3531.5 | 3636.5 |

Warm page-ready median is **3.64 seconds**, with **3.42-3.86 seconds** observed.
HTML response end is only 6.4 ms median. Connection-to-UMAP-value takes 1.164 s
median; UMAP-value-to-observed-visible takes another 1.718 s. These are mixed
client/server intervals, not isolated R function costs.

## Where time goes

Three additional sessions instrumented public Plotly calls separately from the
main timings. Each made one UMAP `newPlot`, one `relayout` and one `react`.
No hidden Browser/MegaBrowser Plotly construction was observed. DT had exactly
one Shiny output value and one draw in every measured and instrumented session.

- **Late Plotly dependency loading:** 729-755 ms between UMAP value arrival and
  `newPlot`. The live page loads the 1,136,953-byte main bundle via late XHR,
  not an initial script tag. In uninstrumented visits it starts around 1.78-1.99 s
  and transfers in 138-195 ms. The rest of the gap includes dependency evaluation
  and widget setup. The newly merged early-loading helper only applies when
  `plot_on_start=TRUE`, consistent with the local launcher's FALSE setting.
- **UMAP rendering:** `newPlot` uses 371-394 ms synchronously for 158 traces and
  3,857 points. The startup relayout adds 28-34 ms and selection reset `react`
  adds 41-51 ms. Promise durations overlap and must not be added as separate costs.
- **DT waits beyond data availability:** first-page responses are about 10 kB.
  Three of five warm requests complete in 59-63 ms, yet first draw occurs much
  later, after Plotly's blocking work. Two requests instead take 1.29 and 1.46 s;
  browser observations cannot assign those delays exclusively to R computation.
  There is no evidence here of repeated table reconstruction or downloading all
  3,857 rows into the browser.
- **Collection initialization:** the pre-output interval is substantial. Source
  inspection shows `RiboCrypt_app()` initializes both MegaBrowser and Observatory
  on the first visit to either tab. Browser network traces confirm dropdown
  requests for hidden MegaBrowser and Observatory Browse controls. UMAP output
  is cached; DT construction is per session. Their precise R costs require a
  separate controlled source-loaded process, not attribution from wall time alone.

## Next experiments

1. Apply early Plotly loading to the Observatory entry path without requiring
   automatic Browser plotting. Measure total load; do not assume the entire
   dependency gap is removable.
2. Separate MegaBrowser initialization from Observatory, and consider deferring
   Observatory Browse until opened. Preserve saved URLs, selection state and
   first-use behavior with the existing browser regression suite.
3. Avoid the redundant dragmode relayout (the R plot already specifies select)
   and investigate skipping an already-clear initial selection reset. Preserve
   explicit subset/reset and restored selections; do not remove synchronization.

## Reproduction

```sh
NODE_PATH=/tmp/observatory-url-browser/node_modules \
  node benchmarks/observatory-startup.cjs 'http://127.0.0.1:4723/#Observatory' 6
PROFILE=1 RESULTS=/tmp/observatory-startup-profile.json \
  NODE_PATH=/tmp/observatory-url-browser/node_modules \
  node benchmarks/observatory-startup.cjs 'http://127.0.0.1:4723/#Observatory' 3
```

Default results: `/tmp/observatory-startup.json`; screenshot:
`/tmp/observatory-startup.png` in this original run (the updated harness now uses
`/tmp/observatory-startup-<port>.png`). Profiling is opt-in and is not used for the
primary timings. Browser contexts were closed afterward; the user's server remains running.
