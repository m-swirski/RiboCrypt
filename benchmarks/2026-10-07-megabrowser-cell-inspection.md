# MegaBrowser region-cell inspection

## Scope

Inspired by the clean-CDS ratio and evidence-audit approach in
[dominant_cell_states](https://github.com/Roleren/dominant_cell_states).
Region/reference distributions now show log2 ratios, and metadata quartiles
are log2-transformed raw-ratio quartiles. Raw ratios and Log2FC are both exported.
Zero numerators have score -Inf and remain counted/included in rank tests but
are omitted from finite box plots; no pseudocount is added.

Only summaries of the already-loaded region density matrix are reused here:
no offline models, replicate pooling or additional FST reads.

Clicking a collapsed-translon cell opens raw-density and reference-ratio
distributions, selected-library records/CSV, and an on-demand metadata-term
table. Metadata comparisons use all analysis libraries, not collapsed group
means. Zero numerators stay valid; non-positive reference densities are undefined.
Term-versus-rest Wilcoxon tests, rank-biserial effects and BH adjustment within
the selected field are exploratory, not study-adjusted evidence. Study counts
and largest-study shares help identify concentration, without correcting it.

Support requires positive coverage sums meeting an editable threshold; defaults
are 10 per library and three supporting libraries per collapsed group. Ratio
support requires both numerator and reference support and does not filter tests.
Coverage sums are not independent read counts. Grey x overlays leave heatmap
values, palettes, cluster membership and full-matrix analyses unchanged.

## Real-data Check

Local source-loaded app, human collection, ATF4 ENST00000674920, T and TC enabled:
3,159 analysis libraries, five clusters and five displayed regions. Labels are
uORF1-uORF4 and clean CDS; source aliases remain available. Clicking uORF4 in
cluster 2 selected 196 libraries in one run and 234 in the final run; cluster
labels/membership are not assumed stable across separate analyses. Both CSVs
were parsed with data.table and checked against density/reference-density
ratios and the 973-nt clean-CDS coverage sum.

Observed click-to-first-distribution rendering ranged from approximately 0.83
to 1.89 seconds in local development runs. These observations are not an SLA or
a controlled comparison. Collection loading and initial clustering are excluded.

The modal is constructed under the app's Bootstrap theme without changing the
session theme. This avoids inactive dynamic tabs and replacement of themed
Selectize resource paths. Support traces are appended to built Plotly data so
they also work with the app's prebuilt heatmap template.
The layout listener also tolerates document-targeted Shiny updates, which would
otherwise interrupt client rendering before the new plot could be displayed.

Final validation after the log2-score update: 3,323 assertions passed with no warnings or skipped tests.
The real-data Playwright harness passed, including unchanged heatmap values
after support updates and working Selectize controls in a subsequent session.
The tutorial HTML was rebuilt from the updated R Markdown source.

## Reproduction

```bash
env -u LC_ALL R --vanilla -q -e 'devtools::load_all("."); testthat::test_local(".")'
NODE_PATH=/path/to/node_modules node benchmarks/megabrowser-cell-details.cjs http://127.0.0.1:7826/
```

The Playwright harness uses ATF4, exercises actual heatmap clicks, reference
switching, metadata tabs, CSV download, desktop/mobile views, support controls,
unchanged heatmap values and initialization of another browser session.
Temporary artifacts are `/tmp/ribocrypt-cell-details-desktop.png`,
`/tmp/ribocrypt-cell-details-mobile.png` and `/tmp/ribocrypt-cell-libraries.csv`.
