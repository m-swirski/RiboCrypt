# RiboCrypt

Interactive visualization and browsing of sequencing coverage, with tools for
ribosome profiling, RNA-seq and other NGS experiments. Use the hosted app at
[RiboCrypt.org](https://ribocrypt.org/) or run it locally with your own
[ORFik experiments](https://bioconductor.org/packages/ORFik).

![RiboCrypt coverage browser](inst/images/tutorial_fig1.png)

## Choose a workflow

| Page | Use it to |
|---|---|
| Browser | Inspect individual-library coverage, reading frames, sequence and annotations; configure Y-axis limits and download displayed coverage. |
| MegaBrowser | Compare precomputed coverage across a collection of libraries, cluster profiles and inspect metadata or coding-region ratios. |
| Observatory | Build library groups with linked UMAP/table selections and compare their merged coverage tracks. |
| Translons | Search predicted translated regions, inspect annotations and download study tables. |
| Analysis and Metadata | Explore metagene, codon and QC results, or find libraries and studies. |

Available organisms, annotations and analyses depend on the data configured for
the server or local app. New features in this checkout may not yet be available
in a Bioconductor release or the hosted app.

## MegaBrowser workflow

1. Choose a collection, gene/transcript and region. Set metadata ordering,
   clustering, normalization and minimum counts, then **Generate Plot**.
2. Inspect the heatmap. **Collapsed clusters** displays one mean profile per
   existing cluster/bin; **Visible groups** focuses the display without
   reclustering. Sidebar visibility, heatmap-only mode, plot height and reset
   zoom controls adjust the layout without reloading coverage.
   **Collapse translons** instead gives equal-width CDS/translon and user-region
   columns in transcript order, preferring clean CDS over its parent CDS.
   Combine it with **Collapsed clusters** for region-relative log2 fold-change
   colours that expose changes in weak uORFs, using the chosen palette and colour
   scale zoom; hover retains raw densities and
   exact fold changes. Neither option changes
   full-matrix clustering or enrichment.
3. Open **Factor enrichment**, then **Group summary** to inspect metadata modes, agreement percentages,
   numeric summaries and library membership for the visible groups. Export
   summaries or memberships as CSV.
4. Within **Factor enrichment**, use **Enrichment** and **Statistics** for metadata enrichment,
   or **Result table** for library records. Change
   **Enrichment metadata**, for example from `TISSUE` to `CELL_LINE`, without
   rebuilding the coverage heatmap.
5. For coding-region comparisons, enable the desired predicted translons (T/TC),
   generate the plot, then open **Translon enrichment**. This analysis starts
   only when its page is opened. It deduplicates identical CDS/translon interval
   sets, retains source aliases and creates `clean_cds` by excluding overlapping
   upstream-ORF bases. Inspect per-library raw coverage-density ratios, group
   distributions, median/IQR summaries and exploratory rank tests with BH
   adjustment. Download ratios, statistics and interval definitions as CSV.
   In **Regions**, use **+** to add a transcript interval with an editable default
   label (`UR1`, `UR2`, etc.). `40:80` includes both endpoints; `40` is equivalent
   to `40:40`. Invalid coordinates show a validation message. Custom regions
   join comparisons and exports, reset on a new analysis and do not alter clean CDS.

Click a collapsed-translon heatmap cell to inspect coverage distributions,
library metadata and log2 ratios against clean CDS, with raw ratios retained in CSV exports. **Ratio by metadata** provides
exploratory term-versus-rest comparisons, support counts and study representation
without loading coverage again. Configurable low-support markers do not alter
scores or analyses. See the tutorial for ratio validity and statistical caveats.

Clustering and analyses use individual count-filtered libraries, not collapsed
rows. Display focus affects the heatmap and Group summary, but not metadata or
translon enrichment. Translon ratios use raw unbinned coverage, not heatmap
logarithms; zero denominators remain undefined without a pseudocount. Rank tests
are exploratory because clusters use the same coverage and libraries are not
necessarily independent biological replicates.

For a worked ATF4 example and complete settings/methods, see the
[app tutorial source](inst/rmd/tutorial.rmd).
The small ORFik example experiment does not include the precomputed collections
and prediction files required by MegaBrowser.

## Installation

Install the released package through Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("RiboCrypt")
```

Use the [Bioconductor installation guide](https://bioconductor.org/install/) to
choose a compatible R/Bioconductor environment; avoid mixing release and
development dependencies. To install this repository's development version:

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("m-swirski/RiboCrypt")
```

## Run locally

Register your sequencing libraries and references as ORFik experiments, then:

```r
library(RiboCrypt)
RiboCrypt_app()
```

The app discovers registered experiments. MegaBrowser and Observatory additionally
need configured all-samples collections, metadata and precomputed coverage;
enabling a page alone does not create those files. Predicted translons require
the corresponding annotation files. See the vignettes for example setup and the
app tutorial for collection-specific requirements.

## Documentation

- [Package overview](vignettes/RiboCrypt_overview.rmd): ORFik example data,
  launching the app, choosing a workflow and plotting from R.
- [Complete app tutorial](inst/rmd/tutorial.rmd): shared source for the in-app
  tutorial and [app vignette](vignettes/RiboCrypt_app_tutorial.rmd).
- [Published Bioconductor documentation](https://bioconductor.org/packages/release/bioc/html/RiboCrypt.html):
  release manual and rendered vignettes, which may describe an older version.
- [Developer benchmarks](benchmarks/README.md): performance reports, scientific
  analysis contracts and reproduction scripts; not installed with the package.

Report bugs through [GitHub issues](https://github.com/m-swirski/RiboCrypt/issues).

## Development

The [accounts and analysis platform plan](docs/development/accounts-and-analysis-platform.md)
tracks the optional identity/dataset-grant implementation and rollout work.
Authentication is opt-in through `RiboCrypt_app(access_control = ...)`; see the
[authentication deployment checklist](docs/deployment/authentication.md) before
configuring an identity provider.

From the repository root, source-load this checkout before testing:

```r
devtools::load_all(".")
testthat::test_local()
```

Edit `inst/rmd/tutorial.rmd` for shared app instructions and regenerate its served
HTML with `rmarkdown::render("inst/rmd/tutorial.rmd", output_format = "html_document")`.
The app vignette includes that same source; keep the overview and README focused
on orientation rather than maintaining separate copies of the complete tutorial.
