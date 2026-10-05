# MegaBrowser output cache and view toolbar

Date: 2026-10-05. Working tree on master, based on `ba01e05` plus the preceding
collapsed-clusters implementation. Nothing in this iteration is pushed.

## Changes

- Cache the rendered interactive heatmap output, including serialization,
  rather than only its R plot object. Keys include plot hash, enrichment mode,
  collapsed state and output namespace. Shiny's existing bounded app cache
  controls retention; no unbounded extra per-browser matrix store is added.
- Ignore size-only and y-only events when synchronizing peer x axes. Previously
  they sent empty JSON arrays to Plotly relayout, causing rejected promises.
- Expose a view toolbar: collapse, sidebar visibility, plot height (400-1200 px),
  reset zoom, actual library/display-row/position counts, and enrichment metadata.
- Keep view size/sidebar controls client-only. Coverage and clustering are not
  recomputed. Responsive margins give small screens more usable plotting width;
  the settings overlay fits the viewport and scrolls internally.
- Give enrichment its own tab rather than placing the chart below the heatmap.
  Shiny can suspend that hidden output until the tab is visited. Statistics and
  Enrichment share the visible metadata selector.

The full matrix remains the analysis source, including statistics and results.
No changes to coverage precision, normalization, filtering or clustering were
made. Static raster output remains uncached because it depends on client size.

## Warm AMD1 comparison

Hardware, software and input settings are the same as in the
[preceding report](2026-10-05-megabrowser-collapsed-clusters.md): human AMD1,
ENST00000451850, leader+CDS, kmer 1, maxNormalized, minimum count 100, five
clusters, 1,754 libraries and 483 positions. Chrome 146, loopback connection,
1440 x 1000, real source-loaded app.

A fresh three-pair baseline was run before this iteration. A separate app
process loaded the updated code and ran the same harness afterward. Each
session generated the full heatmap once before alternating collapsed/full.
Completion is expected row count plus Shiny idle plus 150 ms settling time.
These are warm display toggles, not cold app startup or initial gene loading.

| Display toggle (seconds) | Trial 1 | Trial 2 | Trial 3 | Median |
|---|---:|---:|---:|---:|
| Previous collapsed | 0.353 | 0.480 | 0.498 | 0.480 |
| Updated collapsed | 0.554 | 0.507 | 0.514 | 0.514 |
| Previous full | 1.736 | 1.753 | 1.767 | 1.753 |
| Updated full | 1.191 | 1.098 | 1.093 | 1.098 |

Full-view warm completion improves **37% (1.60x)**. Collapsed completion does
not improve meaningfully in this small sample; its output was already compact.
An earlier updated run measured full median 1.085 s and collapsed median 0.509 s.

The full payload is still approximately 8.55 MB; caching saves server-side JSON
work but does not eliminate network transfer or client rendering. Final median
heatmap `newPlot` time is 501 ms full and 29 ms collapsed. Further full-view
speedups need a separate payload/viewport-rendering design, not statistical
aggregation or indiscriminate numeric rounding.

Changing enrichment metadata to CELL_LINE took 111 ms in the final run and
caused no heatmap redraw. All 1,754 libraries remained in the result table.

## Verification

- Full R suite: **2,972 assertions**, no failures, warnings or skips.
- Reactive regression counters prove restored full output comes from cache and
  collapse/enrichment changes do not rerun coverage/clustering.
- Added size-only/y-only event tests: no empty peer update; sidebar y zoom works.
- Live Chrome smoke checks exercise height, sidebar, reset, enrichment,
  desktop/mobile widths, and settings-panel bounds without JavaScript errors.
- Static-layout coverage checks CSS-variable heights through the spinner, as
  as well as interactive/static renderer switching. Reset is disabled in the
  fixed-axis static view instead of altering the raster framing.
- Tutorial Rmd and served HTML are updated. `git diff --check` is clean.

Run the same toggle benchmark with `megabrowser-clusters.cjs`, and the layout
smoke checks with:

```bash
NODE_PATH=/path/to/node_modules node benchmarks/megabrowser-layout.cjs http://127.0.0.1:7823/
```

Temporary measurement outputs are `/tmp/megabrowser-browser.json` (updated) and
`/tmp/megabrowser-browser-baseline-iteration2.json` (baseline). The layout harness
saves desktop, mobile and static screenshots under `/tmp/megabrowser-layout-*`.
