# Developer benchmarks

Performance reports live here for developers and AI contributors. This directory
is kept in Git and excluded from R source bundles by `^benchmarks$` in
`.Rbuildignore`, so it is not installed with the package.

## Reports

- [2026-10-06: MegaBrowser translon enrichment](2026-10-06-megabrowser-translon-enrichment.md):
  lazy interval-pair ratios, clean CDS, real ATF4 duplicate/overlap checks and timings.

- [2026-10-05: MegaBrowser preparation and group inspection](2026-10-05-megabrowser-preparation.md):
  fewer matrix copies, deferred statistics, group-summary exports and backend timings.

- [2026-10-05: MegaBrowser adaptive rendering and group focus](2026-10-05-megabrowser-group-focus.md):
  canvas/WebGL selection, display-only group focus, heatmap-only layout and warm timings.
- [2026-10-05: MegaBrowser output cache and view toolbar](2026-10-05-megabrowser-view-toolbar.md):
  warm rendering comparison, responsive view controls, and no-op zoom-sync fixes.
- [2026-10-05: MegaBrowser collapsed clusters](2026-10-05-megabrowser-collapsed-clusters.md):
  display-only aggregation, independent enrichment selection, AMD1 timings and profiling.
- [2026-09-25: UMAP startup redraws](2026-09-25-observatory-umap-redraws.md):
  skip no-op Plotly updates, paired timings, regression coverage and rejected experiments.
- [2026-09-25: independent collection startup](2026-09-25-observatory-module-startup.md):
  defer the unused MegaBrowser/Observatory module, with alternating control timings.
- [2026-09-24: Observatory startup](2026-09-24-observatory-startup.md):
  live Select libraries UMAP/DT milestones, dependency loading and initialization candidates.
- [2026-09-24: earlier Plotly loading](2026-09-24-early-plotly-loading.md):
  normal versus deferred dependency loading and a three-way warm-session comparison.
- [2026-09-24: remaining startup costs](2026-09-24-remaining-startup.md):
  current master timings and ranked dependency-loading, hidden-UI and redraw candidates.
- [2026-09-24: sequence startup redraws](2026-09-24-sequence-startup.md):
  remove redundant Plotly mutations, with warm-session comparison and URL-zoom tests.
- [2026-09-24: Plotly rendering profile](2026-09-24-plotly-render-profile.md):
  dependency loading, startup redraws and an experimental SVG comparison.
- [2026-09-24: warm UI response cache](2026-09-24-warm-ui-response.md):
  next startup iteration, alternating control comparison and cache safety tests.
- [2026-09-24: Observatory URL review](2026-09-24-observatory-url-review.md):
  restoration fixes, browser round-trip checks and known validation limits.
- [2026-09-22: ATF4 startup performance branch](2026-09-22-atf4-startup-branch.md):
  implementation plan, fresh control comparison and regression checks.
- [2026-09-22: local ATF4 startup profile](2026-09-22-atf4-startup-local.md):
  first and subsequent user loads, stage timings, hardware/software baseline
  and instructions for comparison with the server.

## Adding a comparison

Use a dated, descriptive Markdown filename. Include the source commit, input
data/settings, hardware, R and relevant package versions, browser/network setup,
cache state, completion criteria, individual trial timings and limitations.
Keep process startup separate from per-session load. Distinguish warm app caches
from warm browser caches. Preserve prior reports rather than replacing baselines.

## Placement convention

R packages use several names for development performance material rather than
one mandatory directory. For example, [vroom's repository documentation](https://github.com/tidyverse/vroom)
points readers to benchmark code in `bench/`. The relevant packaging rule is
explicit exclusion of development-only files through
[`.Rbuildignore`, as described in R Packages](https://r-pkgs.org/structure#sec-rbuildignore).
Here, `benchmarks/` groups human-readable measurements outside `inst/` and tests.
