# Testing the Translation page

Implemented in the server fork after fetching origin over HTTPS (SSH binary unavailable). Starting master/origin commit: `3c72c86510b40ef6d661abc3ee6e5047a34ca12d`. No live RiboCrypt service was replaced or published.

In an RStudio session:

```r
devtools::load_all("~/livemount/forks/RiboCrypt")
data_dir <- "~/livemount/Bio_data/RiboCrypt_translation/current"
ribocrypt_translation_app(data_dir)
```

To include the page in the normal full app, set `RIBOCRYPT_TRANSLATION_DIR` to that directory before calling `RiboCrypt_app()`, or include `translation_data_dir` in `browser_options`. Open **Analysis → Translation**. The Translation server starts only after clicking Translation. Its initial Species overview reads small display sidecars; counts, gene matrices and GO annotations load only when their views need them. Restart the main launcher to pick up the new code and release.

The current human release (`20261015_human_v5`) uses the validated 20261015 matched-count run: **1,321 matched profiles, 10,468 atlas genes, 85 contexts**. Its default Run mapping is an author/metadata-supported proposal, **not proof of shared biological material**. Start with descriptive effects; review biological replication before inferential analysis. For exact-material paired inference, upload your reviewed mapping and check the material-identity control. Both RNA and Ribo technical Run sets must be complete. NA IDs stay unmatched. Independent biological replicates are preserved; multiple distinct assay preparations under one biosample are rejected for review.

The prepared count store contains canonical-gene CDS counts in both assays. Its differential results therefore use Ribo CDS / RNA CDS. The generic collection adapter below defaults to Ribo CDS / RNA mRNA; the result panel displays this distinction.

## Using your all-samples experiments and mapping

```r
ribo <- ORFik::read.experiment("all_samples-Homo_sapiens", validate = FALSE)
rna <- ORFik::read.experiment("YOUR_RNA_COLLECTION_NAME", validate = FALSE)
map <- read.csv("your_run_biosample_mapping.csv", colClasses = "character")
paired <- ribocrypt_translation_pairs(ribo, rna, map)
saveRDS(paired, "your_translation_pairs.rds")
ribocrypt_translation_app(data_dir, paired_counts_file = "your_translation_pairs.rds")
```

Use the actual registered collection names. The required mapping columns are `Run,biosample_id`. Optional `condition,study,replicate,match_level` columns supply reviewed biological metadata. Add `assay_library_id` to identify technical sequencing Runs of the same assay library: these counts are summed within assay/library. Without this field Runs are treated as separate libraries, and ambiguous pooling is rejected. RNA and Ribo are never summed together. For sequencing Runs from different preparations, preserve separate identities and provide the corresponding biological pairing.

The adapter uses count-table colData Run identities, checks genome/species/annotation compatibility and joins regional feature rows by ID. If old collection count tables lack Run colData, pass explicit Run-keyed `ribo_counts` and `rna_counts` matrices. The regional intersection is at the original feature resolution (often transcripts), not an inferred gene-level sum; supply your canonical/gene matrices when that is your intended analysis. A custom paired object must match the atlas species; its differential features may differ from the atlas. For the full app, set `browser_options["translation_pairs_file"]` to the RDS path. Uploaded mapping edits are disabled for this override: rebuild the paired object explicitly. Custom objects/metadata are not persisted to the public fit cache.

## Views and performance

- **Species overview.** Shows curated, count-ready and atlas pair totals, matched RNA–Ribo CLR R² distributions, context coverage and cohort exclusion counts. The current release contains 1,478 curated pairs, 1,430 count-ready pairs, 1,321 atlas profiles, 85 contexts and 163 studies.
- **What changed?** Select a study or multiple studies, biological context, biosamples, reference and comparison. Descriptive mode has no p/FDR. Inferential mode estimates separate assay depth factors and fits RNA, Ribo and joint DESeq2 models. The joint interaction is tested directly; effects are unshrunk. Exact-material mode uses a full-rank sample-blocked interaction; multi-study independent mode adds study to the model and rejects confounding. At least two reviewed biological units per group are required; three or more are preferable. Fits run through Shiny ExtendedTask/future in one background worker. Display thresholds reuse results.
- **Which genes move together?** Signed TE neighbours, RNA comparisons, centered context profiles and a query-centred network. Neighbour queries use a small vectorized calculation, avoiding the cost of importing the full precomputed graph to display 30 neighbours. A negative association is displayed explicitly. These are associations, not measured pathway activation.
- **What functions are suggested?** Gene/pathway exploration, known GO members, TE/RNA term performance and candidate evidence. Default candidates are the 90 rows across 73 genes retained across model sensitivity screens. Broader exploratory mode exposes 1,465 rows. Neither catalogue has candidate-level FDR or a certified novelty claim. Context filtering changes neighbour queries, not precomputed GO validation.

Matched-pair correlation QC is enabled by default at R² ≥ 0.188, the reference cutoff in the [author dataset](https://zenodo.org/records/10373032). Unscored pairs are excluded unless explicitly included. Frozen cohort scores are reused only for the corresponding RNA library; remapped pairs are recomputed. Correlation supports QC and does not establish shared biological material. Changing this filter narrows differential profiles within the released cohort; it does not rebuild the fixed GO reference or restore pairs excluded during atlas construction.

The default “Use all eligible matched biosamples” checkbox avoids rendering hundreds of selected tokens. Manual selection uses server-side search and a scrolling input capped at 160 pixels; biological replicates remain separate.

Data files are external to the git repo: compact atlas, count store, Run mapping, top-50 neighbours and source/output hashes. **265 descriptive study/context contrast requests (149 distinct fits)** have been cached ahead of time. New inferential contrasts are computed on request after replication review. The future/promises dependencies are optional for browsing but required for new background fits; future was installed in this server's R library.

Public fit caches are keyed by actual counts/units, source/reference, full model parameters and engine body. Imported mappings, custom paired objects and access-controlled sessions do not share the public disk fit cache. Changing a cohort invalidates an in-flight result's display. In access-controlled sessions, an atlas must contain only authorized experiments; custom paired-file overrides are disabled until a scoped provider is configured. Permission checks also apply to result downloads.

## Rebuilding releases

Choose a new directory; the builder refuses to overwrite an existing release:

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 Rscript inst/scripts/prepare_translation_atlas.R completed_TE_run new_data_directory
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 Rscript inst/scripts/precompute_translation_contrasts.R new_data_directory
```

For an existing atlas release, create the lightweight display sidecars in a new directory, then regenerate descriptive caches:

```bash
Rscript inst/scripts/prepare_translation_display.R existing_atlas_release new_display_release
Rscript inst/scripts/precompute_translation_contrasts.R new_display_release
```

The display builder splits core and GO data, deduplicates tested-term memberships and creates overview/pair-QC sidecars. It preserves source inputs and existing releases. The builder reads an already validated human TE run. It does not rerun sequencing processing or rewrite source results. Other species can be registered through a `datasets.csv` directory registry (`species,directory`) once compatible atlas/count/mapping bundles exist; only the human bundle has been prepared here. uORF/CAGE/protein integration remains a later extension; this page implements gene-level differential TE, covariation and GO browsing.

Validation includes keyed/technical matching cases, missingness and ambiguity, signed proportionality, fractional-count refusal, confounded designs, synthetic Ribo-specific effects and contrast reversal, real-data pairing/descriptive analysis, exact background-worker agreement, UI rendering and a Chromium smoke test across all three modes. Screenshots and test logs are retained in the matching workspace's versioned implementation evidence directory. An existing access-session call to the unexported `ORFik::organism` was corrected to `BiocGenerics::organism` so access regression tests work with this installed ORFik version.

## Startup benchmark (2026-10-12)

The actual main launcher was tested in Chromium with 1,724 experiments, 69 collections and 13,421 metadata rows. First-session Translation click to gene plots decreased from 11.54 s to 4.05 s (about 65%); the new overview appeared in 2.82 s and opening the gene view then took 0.96 s. Three final sessions showed overview times of 2.82–4.44 s and gene-view times of 0.96–1.07 s. Before/after landing workloads differ because the overview is new; filesystem caches and concurrent load were not controlled.

Direct reads decreased from a 3.01 s monolithic atlas read to 0.002 s for header/overview and 0.039 s for gene matrices. Default differential preview no longer reads and joins the count store on startup. Frozen pair-QC scores were independently reproduced to within 1.12e-15. UI checks confirmed the 160 px selector cap and a change from 715 to 19 eligible profiles when the cutoff was raised to 0.8. Evidence is in `Ribo_RNA_matching/results/benchmarks/20261012_translation_startup/`.

The v4 release includes the37-study count recovery. Another reshift job had47 previously included studies pending at discovery; this release excludes them until authoritative completion returns. All22 TE checks and real-data overview/gene/GO rendering checks passed. Evidence: matching workspace `results/maintenance/20261012_TE_atlas_refresh/`. Existing app sessions need a launcher reload to select v4.

## Stable current release and automatic TE pipeline export (2026-10-15)

`~/livemount/Bio_data/RiboCrypt_translation/current` selects the latest validated app export. The shared `~/livemount/run_RiboCrypt_app.R` launcher uses this path. Future runs of the matching workspace's `scripts/run_human_tec.py` automatically prepare the atlas, display sidecars and descriptive caches, validate keyed counts/QC/neighbours and actual lazy overview/gene/GO rendering, then atomically switch `current`. Export failure leaves the previous app release selected. An unchanged run reuses its exact validated export. The `translation_export` configuration can explicitly disable scratch exports.

Translation sessions resolve the alias once when the module starts, preserving one coherent release throughout that session. Reload the app/browser session to see a later release; the launcher path stays the same. The initial switch from a literal v4 path needs a launcher reload. No package installation is required when using the server fork through `devtools::load_all()`.
