# MegaBrowser collapsed clusters and independent enrichment

Date: 2026-10-05. Source: working-tree changes on master, based on `ba01e05`.
These measurements compare two display modes of the new implementation, not an
old-master versus new-master performance comparison.

## Scope and invariants

- `Collapsed clusters` defaults to FALSE. Each existing cluster/bin becomes one
  displayed mean profile, calculated after the usual normalization/binning.
- Coverage loading, filtering, clustering and statistics use the full matrix.
  Result tables retain individual libraries. Neither collapse nor a new
  enrichment metadata selection reruns coverage loading or clustering.
- Statistics metadata can be selected independently of heatmap ordering;
  categorical missing values are separate and numeric fields use five bins.
- Full display keeps Plotly WebGL. Compact display uses a standard heatmap,
  avoiding half-row gaps observed with WebGL on five-row matrices.
- Summary and annotation outputs are reused on collapse. The plot container is
  rebuilt only when the renderer type changes, not for each heatmap update.

## Machine and software

- Intel Core i7-8850H @ 2.60 GHz, 6 physical cores / 12 threads; 31 GiB RAM.
- Linux x86_64; local reference/collection data on mounted `/media/roler/S`.
- R Under development (unstable), 2024-03-18 r86148.
- Shiny 1.11.1, plotly 4.12.0, ORFik 1.29.11, data.table 1.18.99,
  ComplexHeatmap 2.18.0; Chrome 146.0.7680.164.
- RiboCrypt loaded from this checkout with `devtools::load_all()`; ORFik was
  source-loaded for the app using the local launcher configuration.
- Headless Chrome, 1440 x 1000 viewport, loopback HTTP/WebSocket, no WAN latency.

## Inputs and method

Human `all_samples-Homo_sapiens`, AMD1 (`ENSG00000123505`), transcript
`ENST00000451850`, leader+CDS, kmer 1, maxNormalized, minimum count 100,
five k-means clusters, metadata TISSUE and CELL_LINE, enrichment mode Clusters.
There are **1,754 passing libraries and 483 positions**. Collapsed rendering has
five rows; both modes analyze the identical full library set.

Backend: three repeated full load/normalize/cluster calculations, seed 42,
followed by building/serializing each display. File caches are warm; these are
not cold-disk or app-startup measurements. Backend uses White-Blue colors;
browser uses the launcher's Matrix colors, multiplier 3.

Browser: generate AMD1 once, then alternate collapsed/full three times in the
same session. Completion requires the expected heatmap row count and Shiny idle,
plus a fixed 150 ms settling wait. Includes input scheduling, transport, widget
work and browser rendering. It is not initial page-load time. Repeat trials
can reuse Shiny computation caches but still transfer/render the changed widget.

## Results

| Warm display toggle (seconds) | Trial 1 | Trial 2 | Trial 3 | Median |
|---|---:|---:|---:|---:|
| Collapsed | 0.359 | 0.569 | 0.528 | **0.528** |
| Full | 1.741 | 1.733 | 1.760 | **1.741** |

Median completion is **3.3x faster (70% less time)** with collapsed display.
This is not a guarantee of a sub-0.5-second end-to-end response.

| Backend stage (seconds) | Trial 1 | Trial 2 | Trial 3 |
|---|---:|---:|---:|
| Full load, normalize and cluster | 0.365 | 0.371 | 0.345 |
| Sample grouping metadata | 0.004 | 0.005 | 0.004 |
| Display-only collapse | 0.011 | 0.010 | 0.010 |
| CELL_LINE enrichment | 0.003 | 0.003 | 0.003 |
| Full Plotly build | 0.074 | 0.069 | 0.076 |
| Collapsed Plotly build | 0.008 | 0.008 | 0.007 |
| Full JSON serialization | 0.575 | 0.573 | 0.585 |
| Collapsed JSON serialization | 0.010 | 0.010 | 0.009 |

Full heatmap JSON is 8,403,958 bytes versus 107,271 bytes collapsed, about
**98.7% smaller**. Browser toggle WebSocket frames total 8,551,181 versus
119,022 bytes, including sidebar/output overhead. Plotly's measured heatmap
`newPlot` medians are 521 ms full and 31 ms collapsed. These timings overlap
end-to-end work and should not be added to it.

Changing Statistics from Grouping to CELL_LINE completed in **49 ms** in the
final browser run, with no heatmap Plotly calls. The result table still contained
1,754 libraries. An earlier independent run measured 119 ms for that switch;
one trial is not a service-level guarantee.

## HTML and profiling observations

Screenshots confirmed five continuous, aligned heatmap/sidebar rows. Browser
instrumentation showed only the heatmap and sidebar redrawing on collapse;
summary and gene annotation were not rebuilt. No browser errors or Shiny output
errors occurred. Large full-view long tasks remain during widget processing.

The backend Rprof sampling captured `num_to_char` (JSON numeric serialization),
`fstretrieve`, and matrix transformations. Forced GC from repeated `system.time`
calls accounted for 80.9% of sampled profile time and is excluded from its elapsed
measurements; **do not interpret that proportion as production GC overhead**.
Measured JSON serialization and browser rendering are the clearer bottlenecks.

Further promising work would reduce full-view numeric payload duplication or
render a resolution-limited viewport, but requires a separate accuracy/zoom
design. No rounding, downsampling or statistical changes were introduced here.

## Verification and reproduction

The final full suite passed **2,954 assertions**, with no failures, warnings or
skips. Live Chrome checks also verified collapsed zoom/double-click reset and
the five-row ComplexHeatmap raster. Focused regression tests cover exact cluster membership,
full-library contingency tables, immutable cached metadata, mean profiles,
missing/constant numeric fields, single-position matrices, renderer choice, and
reactive counters proving dropdown/collapse actions do not recompute coverage.
Tutorial Rmd and served HTML were updated.

```bash
env -u LC_ALL R --vanilla -q -e 'devtools::load_all("."); source("benchmarks/megabrowser-clusters.R")'
NODE_PATH=/path/to/node_modules node benchmarks/megabrowser-clusters.cjs http://127.0.0.1:7823/
```

The browser harness requires Playwright and `/opt/google/chrome/chrome`, and a
source-loaded app with the human collection. Outputs are
`/tmp/megabrowser-backend.json`, `/tmp/megabrowser-backend-rprof.out`,
`/tmp/megabrowser-browser.json`, and collapsed/statistics screenshots.
