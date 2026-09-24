# Avoid redundant sequence startup redraws

Date: 2026-09-24. Branch: `perf/plotly-render-profile`.
Control: `89acf9c`, already using SVG for small codon panels.
Candidate: that commit plus the sequence callback change described here.

## Decision

Keep this localized change. Warm median navigation-to-visible improved from
**1.902 to 1.834 seconds** (68 ms / 3.6%). Navigation-to-visible-and-idle improved
from **1.941 to 1.888 seconds** (53 ms / 2.7%). This is a modest incremental
gain, not another large startup reduction or a production latency guarantee.

`inst/js/render_on_zoom.js` now initializes from the actual rendered range,
including URL zoom, rather than a synthetic zero-based full-gene range. It
does not restyle an already centered placeholder or reapply unchanged axis
properties. An unused sequence slice and its zoom-time console logging were
removed. Scientific data, R caches, renderer thresholds, copy behavior and
zoomed-letter rendering remain unchanged.

## Method

Same machine and software baseline as the
[warm UI report](2026-09-24-warm-ui-response.md#environment-and-method):
Intel i7-8850H (6 cores / 12 threads), 31 GiB OS-visible RAM, Ubuntu 24.04.2,
R development revision r86148 (2024-03-18), Shiny 1.11.1, Plotly R 4.12.0,
htmltools 0.5.8.1. Local ORFik source and RiboCrypt were loaded with
`devtools::load_all()`; R commands used `env -u LC_ALL R --vanilla -q -e`.

The control app on 7820 was started and primed before editing JavaScript.
Its script and plot caches retained the old callback. The candidate on 7821
was started after editing. Separate browser instrumentation subsequently
confirmed the old versus new calls; neither server replaced a user's app.
Both used `human_all_merged_l50`, ATF4-ENSG00000128272, ENST00000674920,
RFP and translons enabled, with the human collection available.

After priming both apps, one headless Chrome process made seven alternating
control/candidate rounds, with a fresh empty-cache 1440 x 1000 browser context
per navigation. Round 0 is reported but excluded as browser warmup. No other
browser tests, CPU profiling or R tests ran during this comparison.
Completion uses the existing `warm-session-init.cjs` checks: Plotly fullData
and main SVG, followed by Shiny's busy class clearing. Each visit also checked
ATF4, its transcript, 21 final traces and no JavaScript page errors.

## Timings

Milliseconds from navigation; every measured pair favored the candidate.

| Round | Control visible | Candidate visible | Control idle | Candidate idle |
|---|---:|---:|---:|---:|
| 0, excluded | 1884.0 | 2292.0 | 1948.3 | 2326.1 |
| 1 | 1961.5 | 1839.6 | 2057.8 | 1913.0 |
| 2 | 1878.0 | 1808.9 | 1904.6 | 1844.9 |
| 3 | 1900.6 | 1857.9 | 1933.5 | 1897.6 |
| 4 | 1888.0 | 1843.8 | 1933.0 | 1905.0 |
| 5 | 1969.1 | 1828.6 | 2014.6 | 1878.0 |
| 6 | 1904.0 | 1812.0 | 1948.2 | 1862.1 |
| Median, 1-6 | 1902.3 | 1834.1 | 1940.9 | 1887.8 |

Median plot-value-received-to-visible fell from 983.4 to 907.0 ms. Receipt
itself was similar (909.2 versus 919.3 ms median), consistent with a client-side
improvement. These are small local samples, not confidence intervals; localhost,
headless graphics, fixed ordering and browser compilation affect the result.

Separate instrumentation recorded five startup mutations in every control visit
(add, restyle, relayout, restyle, relayout), versus two in every candidate visit
(add, relayout). The remaining relayout changes only ticks and zeroline on ATF4.
Combined synchronous mutation time was 133.1 / 157.7 / 140.8 ms for control and
71.8 / 69.7 / 70.0 ms for candidate. The last visit in each set had CPU profiling
enabled; these instrumented timings are diagnostic, not the primary comparison.
Promises overlap and their durations must not be summed.

## Regression checks

- Full R suite: 2,777 assertions passed, no warnings, skips or failures.
- Seven Node tests evaluate the shipped callback, including mutation counts,
  narrow initial range, conditional axis styling, pan, zoom/reset, existing
  placeholder and copy behavior. Before the change, six failed; all now pass.
- Chrome: actual codon hover and sequence-copy click, zoom to letters, autorange
  reset to `[1, 2294]`, viewport resize and desktop screenshot inspection passed.
- A fresh URL with `zoom_range=1000:1050` immediately displayed six populated
  sequence traces, no placeholder and the correct range; full-range restoration
  passed. This is exercised by `sequence-startup-smoke.cjs`.
- Observatory URL round trip, restored settings and highlighting, reset/page
  subset, selected-row subset and `go=FALSE` all passed in Chrome.

The copy check deliberately preserves existing slice semantics; it does not
claim to fix sequence-coordinate conventions. Large-gene renderer policy is
unchanged. The resize check is not certification of the entire mobile layout.

## Reproduction

Start separate source-loaded control/candidate apps with the same data and
prime them before measuring. The control must retain its old script: `fetchJS`
and plots are cached per R process, so restarting it from edited sources would
invalidate the comparison.

```sh
RESULTS=/tmp/sequence-warm-comparison.json \
  NODE_PATH=/tmp/observatory-url-browser/node_modules \
  node benchmarks/warm-session-init.cjs http://127.0.0.1:7820/,http://127.0.0.1:7821/ 7
NODE_PATH=/tmp/observatory-url-browser/node_modules \
  node benchmarks/sequence-startup-smoke.cjs http://127.0.0.1:7821/
node --test tests/testthat/js/sequence-startup.cjs
```

Use the existing Plotly profiling and SVG/Observatory smoke harnesses separately
from timed runs. Reload the package and restart the app to test this change;
an existing app can still serve its cached old callback.
