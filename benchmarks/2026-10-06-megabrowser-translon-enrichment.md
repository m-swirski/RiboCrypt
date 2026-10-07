# MegaBrowser translon enrichment

Date: 2026-10-06. Working tree based on `fd18afd` (improved megabrowser).
This report records the initial implementation and the subsequent user-region update.

## Analysis contract

- Uses the same annotation controller and gene-model panel as the annotation track.
  Coordinates are transcript/display-relative and retain exon gaps and strand mapping.
- Deduplicates complete interval sets, retaining T/TC/CDS source aliases. Identical
  outer bounds alone do not identify a duplicate.
- For each CDS, upstream translons starting before that CDS contribute to the union
  subtracted to produce `clean_cds`. Internal/downstream translons are not treated as
  upstream ORFs. Entirely removed CDS regions remain identifiable, with undefined ratios.
- Uses raw unbinned collection coverage, restricted to the original count-filtered
  libraries. Heatmap logarithms, normalization, collapse and display focus are not used
  to compute these ratios. Coverage density is the mean over the region's actual bases;
  each library's ratio is numerator density / denominator density.
- Every unordered region pair is represented once. Pair orientation favors translons
  in the numerator and CDS/clean CDS references in the denominator.
- No pseudocount. Zero denominators give undefined ratios; positive denominators and
  zero numerators give valid zero ratios. Log2 boxes omit zeros/undefined values; the
  statistics table counts them explicitly.
- Cluster summaries use all individual libraries: median, quartiles and proportion
  above one. Exploratory Wilcoxon group-versus-rest tests and rank-biserial effects
  use finite ratios, with BH adjustment across all pairs/groups. Clusters use the same
  coverage and libraries are not independent biological replicates; these are not
  replicate-aware differential-translation tests or independent confirmation.
- Calculation is gated on opening the Translon enrichment tab, before the cache.
  A new analysis while another tab is open does not start computation. Cached results
  can be reused on returning to the tab. Changing collapse/focus does not alter results.
- CSV exports include source labels, all library/pair ratios, statistics and region
  definitions. The selected pair changes only the visible chart/table.

## Real ATF4 case

Human all-samples collection, ATF4 `ENST00000674920`, leader+CDS, T and TC enabled,
minimum counts 100, maxNormalized heatmap, kmer 1, TISSUE/CELL_LINE metadata, five
k-means groups. Backend reproducibility script fixes the clustering seed to 42;
the live app does not fix this seed. 3,159 libraries and 2,228 positions.

| Region | Labels | Length (nt) | Display intervals |
|---|---|---:|---|
| CDS | ENST00000674920 | 1,056 | 1173-2228 |
| Translon | T1, TC1 | 6 | 915-920 |
| Translon | T2, TC2 | 12 | 977-988 |
| Translon | T3, TC3 | 27 | 1046-1072 |
| Translon | T4, TC4 | 180 | 1076-1255 |
| clean_cds | clean_cds (ENST00000674920) | 973 | 1256-2228 |

Four T/TC duplicate pairs become four unique translons. The overlapping upstream
ORF removes exactly 83 nt from the CDS. Six regions yield 15 pairs, 47,385 individual
ratios and 75 pair/group statistics rows. Some ratios are undefined or zero; these
are retained in the downloads and counts, not silently converted to pseudocount ratios.

## Measurements

Intel i7-8850H 2.60 GHz, six cores / 12 threads, 31 GiB RAM. R devel 2024-03-18
r86148; source-loaded RiboCrypt 1.15.2 and local ORFik 1.29.13. Chrome on loopback,
1440 x 1100 desktop and 390 x 844 mobile. No WAN latency. Shiny app source-loaded
with `devtools::load_all()`; installed RiboCrypt is not used.

One initial backend observation: 1.468 s. Final backend observation including CSV
source-label columns: **1.687 s**. Timing covers annotation extraction, raw coverage
reload, densities, all pair ratios and statistics, not initial coverage clustering.
These are observations, not a multi-trial performance claim.

Live first opening in the updated app process: **3.110 s** from tab click to the
ratio chart. An earlier implementation run observed 2.904 s. This includes calculation,
server/client updates and Plotly rendering. No guarantee is made for larger genes or
more regions; pair counts grow quadratically. Ordinary heatmap startup remains lazy
with respect to this analysis.

A second final browser run observed 3.508 s. Cache residency is not assumed across
these new browser sessions; these observations should not be interpreted as a warm
cache guarantee or a controlled performance comparison.

## Verification

### Tab grouping and user-region update

The top-level tabs are Heatmap, Factor enrichment and Translon enrichment.
Factor enrichment contains Enrichment, Group summary, Statistics and Result table.
The Regions **+** modal supports editable UR1/UR2 defaults and one-based inclusive
coordinates. User regions join sequential region IDs, all pairs and CSV exports,
without modifying annotation-derived clean CDS; they reset when the analysis changes.
Raw coverage is reused when adding regions rather than reloaded for each addition.

The updated focused suite passed **65 assertions**, and the full R regression
suite passed. Chrome verified the tab hierarchy, an invalid `40:39` interval,
an edited label for `40:80`, the next UR2 default, and a single-nucleotide `40`
region. Eight regions produced 28 pairs and **88,452** exported library ratios.
No JavaScript or Shiny output errors were observed. One first-tab observation was
**1.720 s**; this is not a controlled comparison with the earlier observations.

### Initial implementation checks

Unit tests cover duplicate aliases, different exon structures, union subtraction,
non-overlapping and fully removed CDS, negative-strand annotation coordinates,
exact raw-coverage means rather than heatmap log ratios, undefined/zero ratios,
rank effects/BH adjustment, tab gating and cache behavior across analysis changes.

Chrome checks the real ATF4 region definitions, five group-statistics rows, stable
statistics when collapsing/focusing groups, all-library CSV export, desktop/mobile
rendering and absence of JavaScript/Shiny output errors. Tutorial Rmd and served
HTML include methods and statistical limitations.

Final full suite: **3,143 assertions**, zero failures, warnings or skips. Focused
translon tests: **39 assertions**, all passing. `git diff --check` is clean.

### Initial collapsed-translon display (superseded below)

An independent, default-off **Collapse translons** control displays raw
per-nucleotide region densities in equal-width columns. It includes unique
CDS/translon and user-defined regions, excluding synthetic clean CDS. The
existing clusters and individual-library analysis remain unchanged. Combined
with cluster collapse, values are arithmetic means of member-library densities.
Summary and annotation tracks share region axes and fixed horizontal margins;
hover labels retain region aliases and coordinates. Raw coverage is shared
with enrichment and loaded only when either analysis/display needs it.

The new focused tests passed **44 assertions**, covering exact exon-aware means,
single-base regions, both collapse modes, focus, preserved clustering, template
hover labels, static rendering objects, reset ranges, validation, lazy loading
and raw-coverage reuse. The full R suite passed. Chrome checked real ATF4 with
3,159 individual rows or five cluster rows and seven region columns (including
two user regions), track alignment, double-click reset, restoration of the
positional heatmap, and desktop/mobile rendering without horizontal overflow.
No JavaScript or Shiny output errors were observed. These checks do not provide
a controlled speed comparison for the new display.

### Transcript ordering and region-relative cluster scores

The current display orders regions by transcript start/end, retains enrichment
region IDs, and replaces a CDS with its own clean CDS where available. A fully
removed clean CDS remains undefined and uses its parent position for ordering;
without a clean CDS the original CDS is retained. User-defined regions join
the same coordinate order. Enrichment comparisons still retain both CDS types.

With both collapse options enabled, colours show log2(cluster mean raw density /
region mean density across all original analysis libraries). The baseline is
library-weighted, not equally weighted across clusters, and remains fixed when
focusing groups. The user-selected colour theme and multiplier map -1/0/+1; colour
values saturate at 0.5x and 2x while hover retains raw means and exact fold
changes. Zero cluster coverage takes the lower colour bound. Zero/undefined
reference coverage stays undefined; no pseudocount is used. This descriptive
score does not correct for library depth or study effects and is not a test.

Focused display tests passed **66 assertions**, including equal colour scores
for equal relative changes in strong and weak regions, clipping/zero handling,
raw hover values, weighted baselines, focus invariance and clean-CDS ordering.
The full R suite passed. Real ATF4 Chrome checks verified clean CDS last after
upstream regions and compared every displayed score with raw cluster density
and the summary-track baseline, in addition to reset and desktop/mobile checks.

Palette selection applies equally to scored and positional heatmaps, including
the existing colour multiplier. Tests compare both built-in palettes at multiple
multiplier settings against the ordinary renderer, verify identical score values,
and check matching static-renderer interpolation and mobile-legend swatches.

Reproduction from repo root:

```bash
env -u LC_ALL R --vanilla -q -e 'devtools::load_all("."); devtools::load_all("/home/roler/Desktop/ORFik"); source("benchmarks/megabrowser-translons.R")'
NODE_PATH=/path/to/node_modules node benchmarks/megabrowser-translons.cjs http://127.0.0.1:7823/
```

Temporary outputs: `/tmp/megabrowser-atf4-translons.rds`,
`/tmp/megabrowser-translon-ratios.csv`, `/tmp/megabrowser-translons-desktop.png`,
`/tmp/megabrowser-translons-mobile.png`.
