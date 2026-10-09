# Warm Browser initial rendering, 2026-10-09

## Scope and implementation

Control: master `3c72c86`. Candidate: branch
`perf/warm-browser-init-2026-10-09`, based on that commit.

Profiling found two avoidable initial Plotly mutations: adding the wide-range
sequence-copy placeholder, then relaying final x-axis styling. Prepare those
states in the R plot instead. Narrow initial URL ranges still use the existing
nucleotide/amino-acid callback. The callback itself, coverage calculations,
WebGL selection, URL readiness, access controls and app caches are unchanged.

This is a modest, localized optimization, not a solution to the remaining
approximately one second before Shiny connects.

## Method

Two separately source-loaded apps, control on 7820 and candidate on 7821, using
`benchmarks/warm-session-launch.R` and the first app definition in
`~/Desktop/run_ribocrypt.R`. Both load local ORFik sources too. Browser is the
initial tab, with real `human_all_merged_l50` data, ATF4-ENSG00000128272,
ENST00000674920, RFP coverage and predicted translons. Each plot has 21 traces.

Prime each app before comparison. Use `warm-session-init.cjs`, fresh Chrome
contexts at 1440 x 1000, empty browser HTTP caches, shared warm R/app caches and
loopback networking. Seven measured pairs visit control first; six further
pairs visit candidate first. Exclude visit zero in each comparison. Timing
runs finished before regression tests were started.

Visible means Plotly full data and main SVG exist; idle additionally means
Shiny's busy class is absent. These are automated readiness milestones, not
instrumented monitor scan-out. Inspect screenshots separately for a correctly
drawn coverage and annotation plot. This is not a production/WAN or concurrency
benchmark; small differences remain sensitive to scheduling and CPU load.

## Results

Milliseconds, medians:

| Milestone | Control | Candidate | Reduction |
| --- | ---: | ---: | ---: |
| Plot visible, combined 13 visits/version | 1616.2 | 1512.9 | 103.3 (6.4%) |
| Shiny idle, combined | 1650.4 | 1560.4 | 90.0 (5.5%) |
| Plot receipt to visible, combined | 237.0 | 168.9 | 68.1 (28.7%) |
| Shiny connected, combined | 1036.3 | 1028.9 | 7.4 |
| Plot received, combined | 1378.6 | 1344.0 | 34.6 |
| HTML response end, combined | 6.6 | 6.5 | 0.1 |

Control-first visible medians: 1548.9 versus 1499.0 ms (3.2%).
Candidate-first visible medians: 1633.0 versus 1540.3 ms (5.7%).
Twelve of thirteen paired visits favor the candidate. Differences before plot
receipt are not caused directly by this rendering change; the strongest
mechanistic evidence is the removal of redraws and the receipt-to-visible gain.

Individual visible timings, control/candidate:

| Order | Pair | Control | Candidate |
| --- | ---: | ---: | ---: |
| Control first | 1 | 1527.3 | 1483.8 |
| Control first | 2 | 1548.9 | 1438.3 |
| Control first | 3 | 1662.1 | 1518.4 |
| Control first | 4 | 1554.9 | 1505.2 |
| Control first | 5 | 1538.7 | 1470.1 |
| Control first | 6 | 1647.6 | 1598.6 |
| Control first | 7 | 1546.9 | 1499.0 |
| Candidate first | 1 | 1616.2 | 1542.1 |
| Candidate first | 2 | 1555.8 | 1619.5 |
| Candidate first | 3 | 1714.4 | 1512.9 |
| Candidate first | 4 | 1756.0 | 1493.3 |
| Candidate first | 5 | 1636.6 | 1606.5 |
| Candidate first | 6 | 1629.3 | 1538.5 |

Exploratory cold-app visible times were 5020.9 ms control and 4949.4 ms
candidate. These single cold trials are excluded and do not establish a cold
startup improvement.

Instrumented `plotly-render-profile.cjs` confirms control calls `newPlot`,
`addTraces` and `relayout` on initialization. Candidate calls only `newPlot`.
Control `newPlot` completion was approximately 208-223 ms; the two candidate
profile visits were 182.0 and 170.6 ms. These instrumented trials are separate
from the timing comparisons; nested call durations must not be added together.

## Verification

- Full `testthat::test_local()` suite passed after `devtools::load_all()`.
- New focused tests cover placeholder coordinates, existing-trace preservation,
  range threshold, fractional bounds, narrow URL ranges and prepared axis style.
- Live `sequence-startup-smoke.cjs`: narrow URL nucleotide traces and full-range
  restoration passed.
- Live `plotly-startup-smoke.cjs`: single early Plotly bundle, fresh/cached
  browser visits and WebGL nucleotide zoom passed.
- Live `aa-svg-smoke.cjs`: hover, full sequence copy, zoom, reset and resize passed.
- Live `local-y-zoom-smoke.cjs`: Browser and Observatory local Y adjustment,
  reset, setting toggle and URL round-trip passed.
- Live `observatory-url-smoke.cjs`: settings/selection restoration, row/page
  subset, reset, `go=FALSE` and deferred manual plotting passed.
- Live `coverage-download-smoke.cjs`: displayed-track CSV downloads in Browser
  and Observatory passed.
- All measured visits had the expected gene/transcript, 21 traces and no
  browser errors. Candidate screenshot was inspected for coverage, annotations,
  sequence placeholder and correctly positioned axes.

## Hardware and software

- Ubuntu 24.04.2 LTS, x86_64.
- Intel Core i7-8850H, 6 cores / 12 threads, 2.60 GHz nominal, 4.30 GHz maximum.
- 31 GiB RAM reported by `free -h` (32 GB installed), 2 GiB swap.
- R Under development (unstable), 2024-03-18 r86148.
- shiny 1.11.1, plotly 4.12.0, htmlwidgets 1.6.4; installed ORFik metadata
  version 1.29.11, with local ORFik source used in app processes.
- Google Chrome 146.0.7680.164, headless; Playwright 1.51.1 on Node 18.19.1.

Raw local artifacts, tooling and R temporary files are on the data drive at
`/media/roler/S/ribocrypt-warm-2026-10-09/`: `comparison.json`,
`reverse-comparison.json`, `baseline.json`, `candidate-prime.json`,
`control-profile.json`, `candidate-profile.json` and app logs. They are not
package assets and are not committed. Launch/measurement scripts remain in
`benchmarks/` to reproduce on another machine.

## Next candidates

The largest remaining stage is navigation to Shiny connection (~1.03 s),
including dependencies, HTML parsing and input binding. Deferring hidden-tab
UI/selectize construction is the next meaningful investigation, but requires
explicit URL restoration, cross-tab navigation and access-control tests.
Connection to plot receipt (~0.3 s) also warrants profiling before changing
reactive readiness. The initial Plotly render itself is now ~0.17 s; do not
trade away large-gene WebGL behavior to chase a smaller startup gain.
