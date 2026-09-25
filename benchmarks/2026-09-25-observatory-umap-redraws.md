# Skip redundant Observatory UMAP updates

Date: 2026-09-25. Branch: `perf/observatory-startup`.
Control: `b0430c9`, including the previous independent collection initialization.
Candidate: only the no-op guards in `inst/js/umap_plot_extension.js`.

## Result

Keep the small callback change. Warm median time to populated visible UMAP/DT
improved from **2.769 to 2.717 s** (52 ms / 1.9%). Readiness plus Shiny idle
improved from **2.847 to 2.778 s** (69 ms / 2.4%). Four of six measured pairs
favored the candidate; this is a modest local improvement, not a major reduction
or a production latency guarantee.

The callback no longer sets dragmode when Plotly already reports `select`.
An already-clear reset returns without `Plotly.react()`. Actual selected points,
including empty selection arrays, and drawn regions still reset. Both raw data
and Plotly's computed data/layout are checked to retain native box/lasso resets.
Event registration, readiness signals, restored selections and server state are
unchanged. No scientific data, trace types, DT code or R startup code changed.

## Measurement

Same machine, R and dependency environment as the
[hardware/software baseline](2026-09-22-atf4-startup-local.md#hardware-and-software-baseline).
Two source-loaded local processes: control 7820, candidate 7821. RiboCrypt and
local ORFik were loaded with `devtools::load_all()`; R was invoked through
`env -u LC_ALL R --vanilla -q -e`. Both used the local launcher metadata,
`plot_on_start=TRUE`, initial Observatory focus and `all_samples-Homo_sapiens`.
The control stayed running with its original cached callback. The candidate was
restarted after removing experimental changes and warmed by regression visits.

One Chrome process alternated candidate then control for seven rounds, each
in a fresh empty-cache 1440 x 1000 context. Round 0 is excluded as warmup.
No unit tests or other browser jobs ran during these timings. Each visit checked
3,857 libraries/points, 158 traces, 15 visible table rows, populated visible
widgets and no JavaScript/Shiny errors. Polling is every 25 ms and cannot run
during main-thread blocking; readiness is not an exact paint timestamp.

Milliseconds from navigation:

| Round | Control both ready | Candidate both ready | Control idle | Candidate idle |
|---|---:|---:|---:|---:|
| 0, excluded | 2806.6 | 3338.5 | 2891.3 | 3364.2 |
| 1 | 2756.1 | 2773.7 | 2756.1 | 2814.7 |
| 2 | 2786.9 | 2654.8 | 2864.4 | 2693.3 |
| 3 | 2814.8 | 2673.8 | 2899.3 | 2709.1 |
| 4 | 2781.0 | 2742.1 | 2857.2 | 2786.5 |
| 5 | 2710.0 | 2692.0 | 2828.1 | 2769.3 |
| 6 | 2754.8 | 2835.2 | 2837.2 | 2962.1 |
| Median, 1-6 | 2768.6 | 2717.1 | 2847.2 | 2777.9 |

All trials, including the slower final candidate visit, are retained. The small
sample and fixed pair order limit precision. Historical 3.64-second measurements
used different launcher settings and are not the control for this patch.

A separate instrumented pair confirmed the control's extra dragmode relayout
(28.8 ms synchronous work) and reset react (35.0 ms) disappeared. Initial
`newPlot` remained approximately 388 ms in both. A later htmlwidgets resize
remains (161 ms control, 200 ms candidate in that pair). These synchronous
durations are diagnostic samples; overlapping promise durations are not summed.

## Rejected experiments

- Deferring Observatory Browse until its tab is opened saved roughly 158 ms
  in an earlier comparison, but the expanded smoke test exposed an empty
  transcript when generating a plot soon after first opening Browse from a
  selected page. The production and unit-test changes were reverted.
- Reserving a root scrollbar gutter did not remove the resize in the profile;
  it was reverted, leaving page layout unchanged. The resize should be traced
  further before assuming that scrollbar appearance causes it.

The selector-URL regression test explicitly clicks Generate Plot after switching
to Browse. Existing semantics deliberately do not auto-plot `view='selector'`,
even when saved browser settings contain `go=TRUE`.

## Validation

- Full R suite: 2,806 assertions, zero failures, warnings or skips.
- Five Node tests evaluate the actual shipped UMAP callback: no-op startup,
  correction of other drag modes, nonempty/empty selectedpoint resets, native
  computed selections/regions and subsequent selection-event propagation.
- Real-app collection smoke passes in both tab-entry orders, including page
  subset, Generate Plot, subset retention, reset and retained MegaBrowser inputs.
- Real-app URL smoke passes saved groups/settings, highlighted selections,
  browser auto-plot, URL round trip, page/selected subsets, reset, `go=FALSE`
  and selector URL followed by manual plotting.
- Desktop screenshot confirms the populated UMAP and library table.
- ATF4 Browser smoke passes hover, sequence copying, zoom, full-range reset
  and desktop-to-mobile resizing without JavaScript errors.

## Reproduction

```sh
RESULTS=/tmp/obs-noop-comparison.json \
  NODE_PATH=/tmp/observatory-url-browser/node_modules \
  node benchmarks/observatory-startup.cjs \
  'http://127.0.0.1:7821/#Observatory,http://127.0.0.1:7820/#Observatory' 7
NODE_PATH=/tmp/observatory-url-browser/node_modules \
  node benchmarks/collection-startup-smoke.cjs http://127.0.0.1:7821/
NODE_PATH=/tmp/observatory-url-browser/node_modules \
  node benchmarks/observatory-url-smoke.cjs http://127.0.0.1:7821/
```

Set `PROFILE=1` for separate diagnostic visits. The harness now records relayout
arguments and call sites, so startup mode updates can be distinguished from
htmlwidgets resizes. Benchmark files are excluded from R builds.
