# Warm ATF4 browser-side Plotly profile

Date: 2026-09-24. Application commit: `6921bc7`.
Profiling branch: `perf/plotly-render-profile`.
No R application code or shipped JavaScript was changed in this investigation.

## Findings

There is useful work left, but the previous approximately 1.2-second interval
between receiving a plot and seeing it was not all Plotly drawing.

1. **Dependencies load after the plot arrives.** The initial document includes
   the widget binding, but Plotly's main bundle is fetched when htmlwidgets calls
   `Shiny.renderDependencies(data.deps)`. In a CPU-profiled baseline visit the
   gap before `newPlot` was 418 ms: the main JS resource took 138 ms to fetch
   (1,136,953 transferred bytes), with subsequent script compilation/evaluation
   and module setup. Other instrumented warm visits had gaps around 0.7-0.8 s.
   The sampled stacks include synchronous jQuery script loading/evaluation.
   Starting the download earlier is worth testing, but merely moving blocking
   script evaluation earlier could delay Shiny connection instead of saving time.
2. **Only one initial plot is constructed.** No duplicate whole-widget
   `newPlot` or `react` occurred during the observed default ATF4 loads. The
   column-switch and x-range-clamp callbacks did not issue startup updates in
   this configuration. Removing those callbacks is not supported by this profile.
3. **The small codon/frame panel starts WebGL.** The initial plot contains 20
   traces, of which three `scattergl` traces belong to the bottom panel (`y4`).
   They have 2, 83 and 245 coordinate entries. Coverage is already SVG in this
   default lines view. `automateTicksAA()` in `R/plotly_helpers.R` converts the
   desktop panel to WebGL, despite its constructors producing SVG traces.
4. **The sequence callback issues five startup mutations.** In
   `inst/js/render_on_zoom.js`, initialization adds a placeholder, immediately
   restyles its unchanged center, and updates axis styling. A subsequent
   relayout moves the center from 1147 to 1147.5 and repeats the axis styling.
   The initial synthetic range starts at zero, while the actual range starts
   at one. These five calls consumed roughly 140-177 ms of synchronous work
   in the alternating run. Not all of that time is automatically removable:
   the placeholder and event bindings still need to exist.

## Browser-only renderer experiment

The profiling harness optionally changes only those three initial `y4` traces
from `scattergl` to `scatter`, in the browser's received object before
`Plotly.newPlot()`. It leaves the server, numerical arrays, layout and application
files unchanged. This is a diagnostic, not an implemented optimization.

Seven independent browser contexts were opened against the same warmed app:
one baseline browser warmup, then three alternating baseline/SVG pairs.
All visits retained ATF4, its expected transcript and 21 final traces, and
reported no JavaScript page errors. The final SVG visit also enabled the Chrome
CPU profiler; its dependency-loading gap was shorter, so end-to-end differences
must not be attributed entirely to the renderer change.

Milliseconds, rounded to one decimal:

| Pair / mode | Before newPlot after value | newPlot synchronous | Startup mutations synchronous | newPlot start to visible | Navigation to complete |
|---|---:|---:|---:|---:|---:|
| 1 baseline | 717.6 | 346.4 | 176.8 | 808.3 | 2606.4 |
| 1 SVG | 679.4 | 154.9 | 152.6 | 347.4 | 1946.0 |
| 2 baseline | 741.1 | 270.0 | 141.4 | 626.5 | 2378.5 |
| 2 SVG | 663.3 | 156.5 | 166.0 | 362.9 | 2003.1 |
| 3 baseline | 742.2 | 258.8 | 145.4 | 655.7 | 2300.2 |
| 3 SVG, CPU profiled | 409.4 | 186.4 | 164.8 | 402.2 | 1743.8 |

Median plot-start-to-visible decreased from **655.7 to 362.9 ms**, about
**293 ms / 45%**. Median navigation-to-complete was 2378.5 versus 1946.0 ms,
but the dependency-loading differences above make this only exploratory evidence
for end-to-end gains. Promise completion times overlap with other work and must
not be summed as independent durations. Synchronous nested `widget.renderValue`
and `newPlot` timings also overlap; the table lists only `newPlot` once.

The first baseline context spent about four seconds synchronously inside
`newPlot`, unlike subsequent contexts. Cold graphics/compiler state in headless
Chrome is a plausible explanation, not a proven diagnosis. This outlier is
excluded from the warm comparison and should not be generalized to real users.

## Visual caveat

The baseline used three canvases; the SVG experiment used none. The SVG
screenshot displayed start/stop-codon marker lines in the colored bottom bands
that were absent from the earlier baseline screenshot on this headless machine.
Thus **visual equivalence is not established**. This could reflect a WebGL
rendering/layering issue, but it requires checking in a normal browser before
claiming a general bug or shipping a renderer change.

Zoomed sequence letters, hover labels, copy-sequence behavior, genomic mode,
long genes, full annotation, extra BigWig tracks and export have not been
regression-tested for the experimental renderer. The existing app was not
changed, so the full R suite was not rerun for this profiling-only task.

## Recommended next iteration

1. Test keeping small codon/frame panels in SVG, including the visual discrepancy
   above. Retain WebGL where trace size and interaction measurements justify it.
   Do not replace all coverage or zoomed sequence traces indiscriminately.
2. Initialize the placeholder/axis state once and avoid no-op startup restyles.
   Use the actual initial range rather than an artificial zero-based range,
   preserving zoom, reset and copy-sequence behavior with browser tests.
3. Separately benchmark earlier Plotly dependency fetching. Check total visible
   load time and Shiny connection, not just whether the post-value interval shrinks.

The first two are concrete, localized candidates. The measured savings should
not be added together without an integrated benchmark.

## Reproduction

Same local source-loaded app/data setup as
[the warm UI profile](2026-09-24-warm-ui-response.md): human merged ATF4, local
Chrome headless, Playwright 1.51.1, fresh 1440 x 1000 browser contexts, localhost
port 7820. The app was started using `devtools::load_all(".")` and the local
ORFik source checkout. User servers and master were not modified.

```sh
NODE_PATH=/tmp/observatory-url-browser/node_modules \
  node benchmarks/plotly-render-profile.cjs http://127.0.0.1:7820/
VISITS=7 COMPARE_SVG=1 NODE_PATH=/tmp/observatory-url-browser/node_modules \
  node benchmarks/plotly-render-profile.cjs http://127.0.0.1:7820/
```

The harness instruments Plotly and htmlwidgets entry points before app scripts
run, records long tasks and resource timing, and CPU-profiles the final visit.
It writes `/tmp/plotly-render-profile.json`, `/tmp/plotly-render.cpuprofile`, and
`/tmp/plotly-render-profile.png`; subsequent runs overwrite those files.
`VISITS` defaults to five. Without `COMPARE_SVG`, every visit is unmodified.
Results are local and instrumented, not production percentiles; GPU configuration,
script compilation caches, machine load and browser differences remain relevant.
