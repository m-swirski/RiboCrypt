# MegaBrowser adaptive rendering and group focus

Date: 2026-10-05. Working tree on master, based on `ba01e05` plus the two
preceding MegaBrowser iterations. This compares against the immediately
preceding working version, not unmodified master. Changes are not pushed.

## Implementation

- Standard Plotly heatmap for matrices of up to 2,000,000 cells; WebGL retained
  above that limit for full views. Collapsed views always use the standard
  renderer. This conservative boundary is a guard, not a universal crossover
  benchmark for every computer or gene.
- `Visible groups` selects complete existing clusters, numeric bins or ordering
  categories. An empty selection restores all groups. Saved row orders are
  remapped, not recomputed; original group labels and member counts remain.
- Focus is applied before optional display collapse. Coverage, normalization,
  clustering, enrichment and result tables retain the original full matrix.
  Summary and gene annotation tracks remain full-collection references.
- `Heatmap only` uses the available height for the heatmap without new coverage
  requests. Disabled in the combined, fixed-axis static view. Preferences are
  scoped to the MegaBrowser workspace.
- Separate analysis/display and layout toolbar rows. Counts distinguish visible
  libraries, total libraries and displayed mean rows. Single focused groups keep
  their sidebar label and result-table link.
- Small plots use integer row ticks. Large plots leave `dtick` unspecified:
  serializing a NULL field as `[]` incorrectly forced unit ticks in Plotly and
  made large views expensive. This issue was caught and fixed before final timing.
- Empty/unset focus states share the same cache key, allowing exact cached
  restoration rather than rebuilding the full output after clearing focus.

No numeric rounding, downsampling or changes to scientific analysis were made.

## Setup

Same machine and software as the
[preceding report](2026-10-05-megabrowser-view-toolbar.md): i7-8850H, 6 cores /
12 threads, 31 GiB RAM; R devel 2024-03-18 r86148; Chrome 146.0.7680.164;
Shiny 1.11.1, plotly 4.12.0 and source-loaded ORFik 1.29.11.
RiboCrypt uses `devtools::load_all()` from this checkout.

Human AMD1, ENST00000451850, leader+CDS, kmer 1, maxNormalized, minimum count
100, five clusters, TISSUE/CELL_LINE ordering, Matrix color theme. Full view:
1,754 libraries x 483 positions (847,182 cells). Collapsed view: five rows.
Headless Chrome, 1440 x 1000, loopback HTTP/WebSocket; no WAN latency.

Generate full AMD1 once, then alternate collapsed/full three times within one
session. Warm app/browser caches; completion is expected row count, Shiny idle,
plus 150 ms settling time. These are display toggles, not cold startup or
first-time coverage computation. Baseline and updated runs use separate app
processes; clustering seeds are not fixed in the live app.

## Timings

| Warm toggle (seconds) | Trial 1 | Trial 2 | Trial 3 | Median |
|---|---:|---:|---:|---:|
| Previous full (WebGL) | 1.154 | 1.098 | 1.073 | 1.098 |
| Updated full | 0.865 | 0.868 | 0.876 | **0.868** |
| Previous collapsed | 0.344 | 0.488 | 0.492 | 0.488 |
| Updated collapsed | 0.421 | 0.264 | 0.235 | **0.264** |

Full refresh is **21% faster**; collapsed median is **46% faster** in this
comparison. An earlier updated run measured medians 0.852 / 0.248 seconds,
supporting the direction of the result. First-collapse cost differs from cached
repeats; this is not a guarantee that every operation completes in 0.26 seconds.

Final heatmap `newPlot` median: 301 ms full and 26 ms collapsed. Previously
full drawing was around 500 ms. Transport is essentially unchanged: about
8.55 MB full versus 119 KB collapsed. The improvement is primarily rendering,
not loss of coverage precision or fewer libraries in the analysis.

CELL_LINE enrichment switching took 145 ms in the final run, with no coverage
heatmap redraw. One live focus example selected 322 of 1,754 libraries in
521 ms, collapsed those to one row, then focused two nonadjacent groups and
restored the full view. The result table still contained all 1,754 libraries.
Focus timing uses a first-use display state and a different completion wait;
do not directly compare it to the warm toggle medians.

## Verification

- Full R suite: **3,010 assertions**, zero failures, warnings or skips.
- Exact single/nonadjacent group membership, saved-order remapping, mean values,
  immutable original tables, original labels and statistics/results invariants.
- Cache counters verify focus/collapse do not reload coverage or rerun clustering;
  clearing focus restores the original cached full output.
- Boundary tests retain WebGL above two million cells. This iteration does not
  benchmark a very large gene such as TTN; its large-view renderer is unchanged.
- Live Chrome: single/multi-group focus and restore, real full-library result
  table, heatmap-only mode, height/sidebar/reset, mobile/settings bounds and
  static-renderer switching. No page/Shiny output errors in passing checks.
- Tutorial Rmd and served HTML updated; `git diff --check` clean.

Reproduction:

```bash
NODE_PATH=/path/to/node_modules node benchmarks/megabrowser-clusters.cjs http://127.0.0.1:7823/
NODE_PATH=/path/to/node_modules node benchmarks/megabrowser-focus.cjs http://127.0.0.1:7823/
NODE_PATH=/path/to/node_modules node benchmarks/megabrowser-layout.cjs http://127.0.0.1:7823/
```

Temporary outputs: `/tmp/megabrowser-browser-baseline-iteration3.json`,
`/tmp/megabrowser-browser.json`, `/tmp/megabrowser-focus.json` and PNG screenshots
under `/tmp/megabrowser-focused.png` and `/tmp/megabrowser-layout-*`.
