# MegaBrowser frame view, 5 November 2026

## Meaning and behavior

The experimental paper Figure 2 B uses three colors, interpolating from white to the phase color; frame corresponds to position modulo three relative to each annotated CDS start. The original MegaBrowser checkbox passed `split_by_frame` to `multiSampleBinRows`, but that function ignored it and ordinary binning mixed all phases.

With frame view enabled, the position heatmap assigns red/green/blue to frames 0/1/2 relative to the **displayed region start**. For a k-mer width greater than one, the width is rounded up to a whole number of codons; each window emits up to three cells, each summing only the corresponding phase. For example, k=9 yields 20 windows -> 60 cells for 180 bases. The x coordinate is the displayed slot within the window, not an assertion that all coverage contributing to that cell is at the single x coordinate. The figure's globally aligned CDS-start frame and this displayed-start frame can differ, e.g. in leaders or genomic views. No bin crosses phase. Existing normalization is applied after binning. Each cell's hover retains the actual normalized value; the encoded color coordinate is only for display.

With **Collapse translons** enabled, each region is calculated from raw unbinned coverage and gets three cells, or one cell if frame is disabled. Each frame is anchored at that translon's start; a clean CDS uses its parent CDS start even if bases were removed. Each value is a mean per nucleotide in that phase. Short or discontiguous regions can have phases with no bases, represented as missing rather than as zero. With **Collapsed clusters** as well, each phase is averaged across the group's libraries, preserving the original group logic. Cluster-versus-all scores remain region/phase-specific and use their original clipped log2 fold changes; hue shows frame and intensity maps scores from -1 (white) to +1 (saturated), with exact mean and fold change available in hover. Ratios, metadata tests and the underlying individual-library table remain available. The frame-off path and pre-existing translon analysis retain their original behavior.

A single discontinuous Plotly colorscale encodes phase and intensity in one heatmap. Actual coverage is supplied separately as hover customdata. Collapsed region scores use fixed limits [-1,1]; ordinary heatmaps scale over the currently displayed matrix. The RGB palette matches the figure. The alternative ComplexHeatmap renderer uses the same encoded values and palette. Frame view necessarily carries three times as many binned positional cells as ordinary k-mer binning, plus hover data; large views should be evaluated for browser memory and latency as well as server preparation.

## Performance

`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript benchmarks/megabrowser-frame.R` checks exact encoding equivalence against the initial implementation and compares five consecutive runs on a shared server. Frame color encoding on a 1,434 × 1,000 matrix had a median **0.091 s -> 0.077 s**. On 2,000 × 1,500, median **0.236 s -> 0.150 s**. A few GC outliers occurred. The optimization avoids materializing an additional full vector of all finite cells and avoids a full-matrix `sweep`; it does not claim the same speedup for data loading or client rendering. Individual timings are in `/tmp/megaframe_bench_final.log` for the originating session; rerun the script for fresh data.

A local Shiny/Playwright smoke test used the public human all-samples collection, ATF4 ENST00000674920, k=9, 3,742 retained libraries and 5 clusters. Frame + cluster view showed 447 phase cells. Collapsing both showed 5 rows and three CDS phase cells; uncollapsing clusters showed 3,742 rows and three cells. Turning off frame and regenerating gave one CDS cell. The color ranges, source values, hover text, legend and absence of browser JS errors were checked. On this shared machine one observed click-to-render measurement for the first frame/cluster view was about 5.9 seconds; it includes app loading/plotly work and is not a controlled before/after benchmark.

To reproduce when local public data and an isolated Playwright/Chromium installation are available:

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript benchmarks/megabrowser-frame-app.R
NODE_PATH=<node_modules> PLAYWRIGHT_BROWSERS_PATH=<browser_directory> node benchmarks/megabrowser-frame-smoke.cjs http://127.0.0.1:7892/
```

The branch does not install ORFik or deploy the public app. Regression tests include normal/short bins, count conservation, true phase-specific region means, overlapping or discontiguous regions, clean-CDS parent-frame anchoring, missing phases, both collapse modes, hover and region cell selection. Existing MegaBrowser regression suites are run as well.

Raw encoder, regression and Shiny smoke logs are in `benchmarks/results/2026-11-05-megabrowser-frames/`.
