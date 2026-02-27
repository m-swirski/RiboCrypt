#' Create RiboCrypt app
#' @param validate.experiments logical, default TRUE, set to FALSE
#' to allow starting the app with malformed experiments, be careful
#' will crash if you try to load that experiment!
#' @param options list of arguments, default
#'  \code{list("launch.browser" = ifelse(interactive(), TRUE, FALSE))}
#' @param all_exp a data.table, default:
#' \code{list.experiments(validate = validate.experiments)}. Which experiments
#' do you want to allow your app to see, default is all in your system config
#' path.
#' @param browser_options named character vector of browser specific arguments:\cr
#' - default_experiment : Which experiment to select, default: first one\cr
#' - default_gene : Which genes to select, default: first one\cr
#' - default_isoform : Which isoform to select, default: first one\cr
#' - default_libs : Which libraries to select: first one, else a single string,
#' where libs are seperated by ",", like "RFP_WT_r1|RFP_WT_r2". Also support
#'  run ids (SRR... etc)\cr
#' - default_kmer : K-mer windowing size, default: 1\cr
#' - default_frame_type : Ribo-seq line type, default: "lines"\cr
#' - default_view_mode : "tx", alternative "genomic"
#' - default_experiment_meta : Which experiment to select for meta analysis,
#'    default: first one\cr
#' - default_gene_meta : Which genes to select for meta analysis,
#'   default: first one\cr
#' - default_isoform_meta : Which isoform to select for meta analysis,
#'    default: first one\cr
#' - translons : Use translon annotation, default "FALSE"\cr
#' - plot_on_start : Plot when starting, default: "FALSE"\cr
#' - hide_settings : Hide settings bar in browser on start, default "TRUE"\cr
#' @param init_tab_focus character, default "browser". Which tab to open on
#' init.
#' @param metadata a path to csv or a data.table of metadata columns,
#' must contain a "Run" column to merge IDs to ORFik experiments.
#' It is used in the metabrowser tab for grouping of samples.
#' @param all_exp_meta a data.table, default:
#' \code{all_exp[grep("all_samples-", name),]}. Can also be NULL, to ignore
#' the metabrowser completly.
#' It is the subset of all_exp which are collections (the set of
#' all experiments per organism), this will be fed to the metabrowser, while
#' remaining all_exp are used in all other modules.
#' @import shiny bslib ORFik NGLVieweR ggplot2 fst rclipboard data.table
#' @importFrom Biostrings strsplit width
#' @importFrom GenomeInfoDb seqlevelsStyle seqlevelsStyle<-
#' @importFrom shinycssloaders withSpinner
#' @importFrom markdown mark_html
#' @importFrom shinyjqui jqui_resizable jqui_draggable
#' @importFrom shinyjs click useShinyjs hideElement showElement
#' @importFrom shinyWidgets prettySwitch
#' @importFrom knitr knit
#' @importFrom stringr str_sub str_count
#' @importFrom httr GET write_disk
#' @importFrom ComplexHeatmap Heatmap row_order draw
#' @importFrom DT datatable renderDT DTOutput formatStyle styleInterval
#' @importFrom grid grid.grabExpr
#' @importFrom rtracklayer BigWigFile
#' @importFrom writexl write_xlsx
#' @return RiboCrypt shiny app
#' @export
#' @examples
#' run_variable <- 1 # Ignore check test limit
#' ## Default run
#' # RiboCrypt_app()
#' ## Plot on start
#' # RiboCrypt_app(browser_options = c(plot_on_start = "TRUE"))
#' ## Init with an experiment and gene (you must of course have the experiment)
#'
#' # RiboCrypt_app(validate.experiments = FALSE,
#' #       browser_options = c(plot_on_start = "TRUE",
#' #                           default_experiment = "all_merged-Homo_sapiens_2024_8",
#' #                           default_gene = "ATF4-ENSG00000128272"))
#' # RiboCrypt_app(validate.experiments = FALSE, all_exp = all_exp,
#' # browser_options = c(plot_on_start = "TRUE",
#' #                    default_experiment = "human_all_merged_l50",
#' #                    default_gene = "RPL12-ENSG00000197958",
#' #                    default_isoform = "ENST00000361436",
#' #                    default_view_mode = "genomic"))
#' # RiboCrypt_app(validate.experiments = FALSE,
#' #       browser_options = c(plot_on_start = "TRUE",
#' #                           default_experiment = "all_merged-Saccharomyces_cerevisiae",
#' #                           default_gene = "EFM5-YGR001",
#' #                           default_view_mode = "genomic"))
RiboCrypt_app <- function(
  validate.experiments = TRUE,
  options = list("launch.browser" = ifelse(interactive(), TRUE, FALSE)),
  all_exp = list.experiments(validate = validate.experiments),
  browser_options = c(), init_tab_focus = "browser",
  metadata = NULL, all_exp_meta = all_exp[grep("all_samples-", name), ]
) {
  rc_parameter_setup()
  # User interface
  ui <- tagList(
    rc_header_styling(),
    rc_header_image(),
    helper_button_redirect_call(),
    navbarPage(
      id = "navbarID",
      lang = "en",
      windowTitle = "RiboCrypt",
      title = rc_title(),
      theme = rc_theme(),
      selected = init_tab_focus,
      browser_ui("browser", all_exp, browser_options, names_init, libs),
      browser_allsamp_ui(
        "browser_allsamp", all_exp_meta, browser_options, metadata,
        names_init_meta
      ),
      observatory_ui("observatory", all_exp_meta),
      analysis_ui("Analysis", all_exp, browser_options, libs, metadata),
      metadata_ui("metadata", all_exp, all_exp_meta),
      tutorial_ui("tutorial")
    )
  )
  cat("Done (UI setup):")
  print(round(Sys.time() - time_before, 2))

  server <- function(input, output, session) {
    reactive_url()
    org_and_study_changed_checker(input, output, session)
    tutorial_server("tutorial")
    rv <- browser_server(
      "browser", all_exp, without_readlengths_env, df,
      experiments, tx, cds, libs, org, gene_name_list,
      gg_theme, rv, browser_options
    )
    if (nrow(all_exp_meta) > 0) {
      browser_allsamp_server(
        "browser_allsamp", all_exp_meta, df_meta, metadata,
        names_init_meta, browser_options, exp_init_meta,
        exps_dir
      )
    } else {
      print("No MegaBrowser exps given, ignoring MegaBrowser server.")
    }
    rv <- analysis_server(
      "Analysis", all_exp, without_readlengths_env,
      with_readlengths_env, df, df_with, experiments, tx, cds, libs, org,
      gene_name_list, rv, metadata, names_init, browser_options
    )
    metadata_server(
      "metadata",
      all_exp, metadata,
      all_exp_meta,
      browser_options
    )
    observatory_server(
      "observatory",
      meta_experiment_list = all_exp_meta,
      all_libraries_df = metadata
    )

    cat("Server this: ")
    print(round(Sys.time() - this_time_before, 2))
    cat("Server total: ")
    print(round(Sys.time() - time_before, 2))
  }
  shinyApp(ui, server, options = options)
}
