rc_parameter_setup <- function() {
  with(rlang::caller_env(), {
    time_before <- Sys.time()
    print(paste("- Starting app instance at time:", format(Sys.time(), "%Y-%m-%d %H:%M")))
    stopifnot(is(all_exp, "data.table"))
    stopifnot(!is.null(all_exp$name))
    stopifnot(nrow(all_exp) > 0)

    if (!isTruthy(browser_options["allow_collections_in_browser"])) {
      browser_options["allow_collections_in_browser"] <- TRUE
    }

    if (is.null(all_exp_meta)) {
      all_exp_meta <- data.table(name = character(), organism = character())
    } else stopifnot(is(all_exp_meta, "data.table"))
    allow_collections_in_browser <- browser_options["allow_collections_in_browser"]
    if (nrow(all_exp_meta) > 0 & !as.logical(allow_collections_in_browser)) {
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

    if (!isTruthy(browser_options["default_experiment_meta"])) {
      browser_options["default_experiment_meta"] <- all_exp_meta$name[1]
    }
    if (isTruthy(browser_options["default_experiment_meta"]) &
        nrow(all_exp_meta) > 0) {
      stopifnot(browser_options["default_experiment_meta"] %in% all_exp_meta$name)
    }
    if (!isTruthy(browser_options["default_experiment_translon"])) {
      default <- "all_merged-Homo_sapiens"
      selected_exp <- ifelse(default %in% all_exp$name, default, "AUTO")
      if (selected_exp == "AUTO") {
        default <- "human_all_merged_l50"
        selected_exp <- ifelse(default %in% all_exp$name, default, "AUTO")
      }
      browser_options["default_experiment_translon"] <- selected_exp
    }


    if (is.na(browser_options["plot_on_start"])) {
      browser_options["plot_on_start"] <- FALSE
    }
    if (!isTruthy(browser_options["default_view_mode"])) {
      browser_options["default_view_mode"] <- "tx"
    }
    if (!isTruthy(browser_options["url_api_call"])) {
      browser_options["url_api_call"] <- ""
    }

    stopifnot(is.character(browser_options["default_view_mode"]) &
              browser_options["default_view_mode"] %in% c("tx", "genomic"))
    if (!isTruthy(browser_options["collapsed_introns_width"])) {
      browser_options["collapsed_introns_width"] <- "30"
    }
    stopifnot(!is.na(as.numeric(browser_options["collapsed_introns_width"])))
    if (!isTruthy(browser_options["full_annotation"])) {
      browser_options["full_annotation"] <- FALSE
    }

    if (!isTruthy(browser_options["translons"])) {
      browser_options["translons"] <- FALSE
    }
    if (!isTruthy(browser_options["translons_transcode"])) {
      browser_options["translons_transcode"] <- FALSE
    }

    if (!isTruthy(browser_options["search_on_init"])) {
      browser_options["search_on_init"] <- ""
    }

    if (!isTruthy(browser_options["allow_non_bw"])) {
      browser_options["allow_non_bw"] <- FALSE
    }
    exps_dir <- ORFik::config()["exp"]
    exp_init <- read.experiment(browser_options["default_experiment"],
                                validate = FALSE, in.dir = exps_dir)
    names_init <- get_gene_name_categories(exp_init)
    if (!isTruthy(browser_options["default_gene"])) {
      browser_options["default_gene"] <- names_init$label[1]
    }
    stopifnot(browser_options["default_gene"] %in% names_init$label)

    tx_init <- loadRegion(exp_init)
    cds_init <- loadRegion(exp_init, "cds")

    exp_init_meta <- names_init_meta <- NULL

    if (nrow(all_exp_meta) > 0) {
      meta_org <- all_exp_meta[name == browser_options["default_experiment_meta"]]$organism[1]
      browser_org <- all_exp[name == browser_options["default_experiment"]]$organism[1]
      exp_init_meta <- read.experiment(browser_options["default_experiment_meta"],
                                       validate = FALSE, in.dir = exps_dir)
      names_init_meta <- if (meta_org == browser_org) {
        copy(names_init)
        } else {
        get_gene_name_categories(exp_init_meta)
      }
      if (!isTruthy(browser_options["default_gene_meta"])) {
        browser_options["default_gene_meta"] <- names_init_meta$label[1]
      }
      stopifnot(browser_options["default_gene_meta"] %in% names_init_meta$label)

      gene_isoforms_meta <- names_init_meta[label == browser_options["default_gene"],]
      if (!isTruthy(browser_options["default_isoform_meta"])) {
        browser_options["default_isoform_meta"] <- gene_isoforms_meta$value[1]
      }
      stopifnot(browser_options["default_isoform_meta"] %in% gene_isoforms_meta$value)
    }

    gene_isoforms <- names_init[label == browser_options["default_gene"],]
    if (nrow(gene_isoforms) == 0) stop("Selected gene has no isoforms!")
    if (!isTruthy(browser_options["default_isoform"])) {
      browser_options["default_isoform"] <- gene_isoforms$value[1]
    }
    stopifnot(browser_options["default_isoform"] %in% gene_isoforms$value)

    if (!isTruthy(browser_options["hide_settings"])) {
      browser_options["hide_settings"] <- TRUE
    }
    if (!isTruthy(browser_options["default_kmer"])) {
      browser_options["default_kmer"] <- 1
    } else {
      stopifnot(!is.na(as.numeric(browser_options["default_kmer"])))
    }
    if (!isTruthy(browser_options["codon_filter_count"])) {
      browser_options["codon_filter_count"] <- 1000
    } else {
      stopifnot(!is.na(as.numeric(browser_options["codon_filter_count"])))
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
      default_libs <- libraries_string_split(browser_options["default_libs"], libs)
    }

    gg_theme <- gg_theme_template()

    cat("Done (Parameter setup):"); print(round(Sys.time() - time_before, 2))
  }
  )
}
