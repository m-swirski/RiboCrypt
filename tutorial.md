---
title: "Tutorial"
output: 
  BiocStyle::html_document:
    self_contained: true
---




# Introduction

Welcome to `RiboCrypt` 
`RiboCrypt` is an R package for interactive visualization in genomics. `RiboCrypt` works with any NGS-based method, but much emphasis is put on Ribo-seq data visualization. 

This tutorial will walk you through usage of the app. 

`RibCrypt` app currently supports creating interactive browser views for NGS tracks, using
ORFik, Ribocrypt and massiveNGSpipe as backend.


# Browser
The browser is the main coverage plot display page. 
It contains a click panel on the left side and display panels on the right.
It displays coverage of NGS data in either transcript coordinates (default),
or genomic coordinates (like IGV).
Each part will now be explained:

<img src="../images/tutorial_browser_browser.png" alt="plot of chunk unnamed-chunk-1" width="1000px" height="500px" style="display: block; margin: auto auto auto 0;" />

## Display panel (browser)
The display panel shows what can be specified to display, the possible
select boxes are:
1. Select an organism: Either select "ALL" to keep all experiments, or select a specific organism.
2. Select an experiment: The experiments contain study names combined with organism (some studies
are multi species, so sometimes one study have multiple experiments). 
Select which one you want. There also exist merged experiments 
(all samples merged for the organism, etc)
3. Select a gene: A gene can be selected currently using:
 - Transcript id (ENSEMBL)
 - Gene id (ENSEMBL) (To be implemented)
 - Gene symbol (hgnc, etc) (To be implemented)
4. Select libraries
Each experiment usually have multiple libraries. Select which one
to display, by default if you select multiple libraries they will be shown
under each other.
5. Select frames display type: 
  -lines (single line, most clear for middle distance (> 100 nt))
  -columns (single point bars, most clear for single nt resolution)
  -stacks (Area under curve, stacked, most clear for long distance (> 1000 nt))
  -area (Area under curve, with alpha (see-through), most clear for long distance (> 1000 nt))
6. K-mer length:
When looking at a large region (> 100nt), pure coverage can usually be hard to inspect.
Using K-mer length > 1 (9 is a good starting point to try), you can easily look at
patterns over larger regions. 

## Display panel (Navigate)
Here additional options are shown:
1. 5' extension (extend viewed window upstream, outside defined region)
2. 3' extension (extend viewed window downstream, outside defined region)
3. Genomic View (Activate/deactivate genomic view, giving splice information and
correct positions in genome, but a lot harder to understand)
4. Protein structures (If you click the annotation name of a transcript
in the plot panel it will display the alpha-fold protein colored by the ribo-seq data displayed
in the plot panel)

## Plot panel
From the options specified in the display panel, when you press "plot" the data will be displayed.
It contains the specific parts:

1. Ribo-seq data (top), the single or multi-track data is displayed on top. By default
Ribo-seq is displayed in 3 colors, where 
 - red is 0 frame, the start frame of reference transcript.
 - green is +1 frame
 - blue is +2 frame 
2. Sequence track (top middle), displayes DNA sequence when zoomed in (< 100nt)
3. Annotation track (middle), the annotation track displays the transcript annotation, together
with black bars that is displayed on top of the data track.
4. Frame track (bottom), the 3 frames displayed with given color bars:
 - white  (Start codons)
 - black  (Stop codons)
 - purple (Custom motifs)
 When zoomed in, the amino acid sequence is displaced within each frame


# Heatmap

This tab displays a heatmap of coverage per readlength at a specific
region (like start site of coding sequences) over all genes selected. 

<img src="../images/tutorial_heatmap_browser.png" alt="plot of chunk unnamed-chunk-2" width="1000px" height="500px" style="display: block; margin: auto auto auto 0;" />

## Display panel (heatmap)
The display panel shows what can be specified to display, the possible
select boxes are:
1. Select an organism: same as for browser
2. Select an experiment: same as for browser
3. Select a gene: A gene can be selected currently using:
 - all (all transcripts)
 - Transcript id (ENSEMBL)
 - Gene id (ENSEMBL) (To be implemented)
 - Gene symbol (hgnc, etc) (To be implemented)
4. Select libraries
Only 1 library can be selected currently in heatmap mode.
5. View region
Select one of:
 - Start codon
 - Stop codon
6. Normalization
Normalization mode for data display, select one of:
 - transcriptNormalized (each gene counts sum to 1)
 - zscore (zscore normalization, will give better overview if 1 readlength
 is extreme)
 - sum (raw sum of counts, is very sensitive to extreme peaks)
 - log10sum (log10 sum of counts, is less sensitive to extreme peaks)
7. Min Readlength
 The minimum readlength to display from library
8. Max readlength 
 The maximum readlength to display from library

## Display panel (Navigate)
Here additional options are shown:
1. 5' extension (extend viewed window upstream from point, default 30)
2. 3' extension (extend viewed window downstreamfrom point, default 30)

## Plot panel
From the options specified in the display panel, when you press "plot" the data will be displayed.
It contains the specific parts:

1. Y-axis: the read lengths
2. X-axis: the position (in nucleotides)
3. Color: Per x,y coordinate, color by the count score

# metadata

Metadata tab displays information about studies.

## Study accession number
 Here you input a study accession number in the form of either:
 - SRP
 - GEO (GSE)
 - PRJNA (PRJ....)
 - PRJID (Only numbers)

## Output
On top the abstract of the study is displayed, and on bottom a table
of all metadata found from the study is displayed.

# Additional information

## massiveNGSpipe

The processing pipeline used is massiveNGSpipe which wraps over multiple tools:
1.  Fastq files are download with ORFik download.sra
2.  Adapter is detected with either fastqc (sequence detection) and falls
    back to fastp auto detection. 
3.  Reads are then trimmed with fastp
4.  Read are collapsed (get the set of unique reads and put duplication count in read header)
5. Reads are aligned with the STAR aligner (using the wrapper in ORFik), that supports
contamination removal. 
6. When all samples of study are aligned, an ORFik experiment is created that
connects each sample to metadata (condition, inhibitor, fraction, replicate etc)
7. Bam files are then converted to ORFik ofst format
8. These ofst files are then pshifted
9. Faster formats are then created (bigwig and covRLE) for faster visualization

## Introduction to Ribo-seq

If you're not familiar with terms like "p-shifting" or "p-site offset", it's best to walk through ORFikOverview vignette, especially chapter 6 "RiboSeq footprints automatic shift detection and shifting"

https://bioconductor.org/packages/release/bioc/vignettes/ORFik/inst/doc/ORFikOverview.html#riboseq-footprints-automatic-shift-detection-and-shifting

## About

This app is created as a collaboration with:

- University of Warsaw, Poland
- University of Bergen, Norway

Main authors and contact:

- Michal Swirski (Warsaw), email: michal.swirski@uw.edu.pl
- Håkon Tjeldnes (Bergen), email: hauken_heyken@hotmail.com
