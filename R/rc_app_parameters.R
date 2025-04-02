rc_parameter_setup <- function() {
  with(rlang::caller_env(), {
    time_before <- Sys.time()
    print(paste("- Starting app instance at time:", format(Sys.time(), "%Y-%m-%d %H:%M")))
    stopifnot(is(all_exp, "data.table"))
    stopifnot(!is.null(all_exp$name))
    stopifnot(nrow(all_exp) > 0)

    if (is.null(all_exp_meta)) {
      all_exp_meta <- data.table(name = character(), organism = character())
    } else stopifnot(is(all_exp_meta, "data.table"))
    if (nrow(all_exp_meta) > 0) {
      all_exp <- all_exp[!(name %in% all_exp_meta$name),]
    }
    print(paste("Running with", nrow(all_exp), "experiments"))
    print(paste("Running with", nrow(all_exp_meta), "collections"))

    if (!is.null(metadata)) {
      if (is.character(metadata)) metadata <- fread(metadata)
      stopifnot(is(metadata, "data.table"))
      columns_to_show <- c("study_accession", "Run", "ScientificName", "sample_title", "BioProject",
                           "LIBRARYTYPE", "REPLICATE", "CONDITION", "INHIBITOR",
                           "BATCH", "TIMEPOINT", "TISSUE", "CELL_LINE", "GENE", "FRACTION")
      metadata <- metadata[, colnames(metadata) %in% columns_to_show, with = FALSE]
      print(paste("Running with", nrow(metadata), "metadata rows"))
    }
    # Set environments
    with_readlengths_env <- new.env()
    without_readlengths_env <- new.env()
    #with_cigar_env <- new.env() # Not used for now
    # Add resource directories
    addResourcePath(prefix = "images",
                    directoryPath = system.file("images", package = "RiboCrypt"))
    addResourcePath(prefix = "rmd",
                    directoryPath = system.file("rmd", package = "RiboCrypt"))
    # Setup variables
    if (!isTruthy(browser_options["default_experiment"])) {
      browser_options["default_experiment"] <- all_exp$name[1]
    }

    if (!isTruthy(browser_options["default_experiment_meta"]) &
        nrow(all_exp_meta) > 1) {
      browser_options["default_experiment_meta"] <- all_exp_meta$name[1]
    }
    if (isTruthy(browser_options["default_experiment_meta"]) &
        nrow(all_exp_meta) > 0) {
      stopifnot(browser_options["default_experiment_meta"] %in% all_exp_meta$name)
    }
    if (is.na(browser_options["plot_on_start"])) {
      browser_options["plot_on_start"] <- FALSE
    }
    if (!isTruthy(browser_options["default_view_mode"])) {
      browser_options["default_view_mode"] <- "tx"
    }
    stopifnot(is.character(browser_options["default_view_mode"]) &
              browser_options["default_view_mode"] %in% c("tx", "genomic"))
    if (!isTruthy(browser_options["allow_non_bw"])) {
      browser_options["allow_non_bw"] <- FALSE
    }
    exp_init <- read.experiment(browser_options["default_experiment"],
                                validate = FALSE)
    names_init <- get_gene_name_categories(exp_init)
    if (!isTruthy(browser_options["default_gene"])) {
      browser_options["default_gene"] <- names_init$label[1]
    }
    stopifnot(browser_options["default_gene"] %in% names_init$label)
    if (!isTruthy(browser_options["default_gene_meta"])) {
      browser_options["default_gene_meta"] <- names_init$label[1]
    }
    stopifnot(browser_options["default_gene_meta"] %in% names_init$label)

    gene_isoforms <- names_init[label == browser_options["default_gene"],]
    if (nrow(gene_isoforms) == 0) stop("Selected gene has no isoforms!")
    if (!isTruthy(browser_options["default_isoform"])) {
      browser_options["default_isoform"] <- gene_isoforms$value[1]
    }
    stopifnot(browser_options["default_isoform"] %in% gene_isoforms$value)

    if (!isTruthy(browser_options["default_kmer"])) {
      browser_options["default_kmer"] <- 1
    } else {
      stopifnot(!is.na(as.numeric(browser_options["default_kmer"])))
    }
    if (!isTruthy(browser_options["default_frame_type"])) {
      browser_options["default_frame_type"] <- "lines"
    } else {
      stopifnot(is.character(browser_options["default_frame_type"]))
    }
    libs <- bamVarName(exp_init)
    if (!isTruthy(browser_options["default_libs"])) {
      browser_options["default_libs"] <- libs[1]
    } else {
      default_libs <- unlist(strsplit(browser_options["default_libs"], "\\|"))
      if (!all(default_libs %in% libs))
        stop("You defined default_libs, but some of those are not valid names,",
             " in selected experiment!")
    }
    cat("Done (Parameter setup):"); print(round(Sys.time() - time_before, 2))
  }
  )
}
