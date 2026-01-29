fetchJS <- function(script_name) {
  # if (script_name == "render_on_zoom.js") {
  #   message("Using local render script!")
  #   script <- "~/Desktop/forks/RiboCrypt/inst/js/render_on_zoom.js"
  # } else
  script <- system.file("js", script_name, package =
                         "RiboCrypt")
  lines <- readLines(script)
  paste(lines, sep = "", collapse = "\n")
}

#' Fetch Javascript sequence
#'
#' @param target_seq the target sequence
#' @param nplots number of plots
#' @param distance numeric, default 250 When to show sequence, when x-range is
#' <= to this number.
#' @param aa_letter_code character, default: "one_letter", alternative: three_letters
#' @param input_id shiny id of the object
#' @return a list of 2 lists, the nt list (per frame, total 3)
#'  and AA list (per frame, total 3)
#' @importFrom Biostrings AMINO_ACID_CODE
fetch_JS_seq <- function(target_seq, nplots, distance = 250,
                         aa_letter_code = "one_letter", input_id, frame_colors = "R") {
  display_dist <- nchar(target_seq)
  fr_colors <- frame_color_themes(frame_colors)
  nt_yaxis <- paste0("y", nplots + 1)
  aa_yaxis <- paste0("y", nplots + 3)
  nts <- lapply(1:3, function(x) seq(x, display_dist, 3))
  rendered_seq <- strsplit(as.character(target_seq),"")[[1]]
  translate_fuzzy_logic <- ifelse(all(rendered_seq %in% DNA_BASES), "error", "X")
  aas <- lapply(1:3, function(x) suppressWarnings(
    strsplit(as.character(translate(target_seq[[1]][x:display_dist],
                                    if.fuzzy.codon = translate_fuzzy_logic)), "")[[1]]))
  if (aa_letter_code == "three_letters") {
    aa_code <- AMINO_ACID_CODE
    aa_code["*"] <- "*"
    aas <- lapply(aas, function(x) aa_code[match(x, names(aa_code))] )
  }
  nts_js_data <- lapply(1:3, function(fr) list(x = nts[[fr]],
                                               text = rendered_seq[nts[[fr]]],
                                               y = rep(0.5, length(nts[[fr]])),
                                               xaxis = "x",
                                               yaxis = nt_yaxis,
                                               distance = distance,
                                               color = fr_colors[fr]))
  aa_js_data <- lapply(1:3, function(fr) list(x = nts[[fr]][seq_along(aas[[fr]])],
                                              text = aas[[fr]],
                                              y = rep(2.4 - fr, length(aas[[fr]])),
                                              xaxis = "x",
                                              yaxis = aa_yaxis,
                                              distance = distance,
                                              color = "grey45"))

  render_on_zoom_data <- c(nts_js_data, aa_js_data)

  return(list(traces = render_on_zoom_data,
              sequence = as.character(target_seq),
              input_id = paste0(input_id, "_copy")))
}

addJSrender <- function(multiomics_plot, target_seq, nplots, seq_render_dist,
                        aa_letter_code, input_id, frame_colors) {
  render_on_zoom_data <- fetch_JS_seq(target_seq = target_seq, nplots = nplots,
                                      distance = seq_render_dist,
                                      aa_letter_code = aa_letter_code, input_id, frame_colors)
  select_region_on_click_data <- list(nplots = nplots, input_id = input_id)
  multiomics_plot <- multiomics_plot %>%
    onRender(fetchJS("render_on_zoom.js"), render_on_zoom_data) %>%
    onRender(fetchJS("select_region_on_click.js"), select_region_on_click_data)
  return(multiomics_plot)
}

helper_button_redirect_call <- function() {
  tabPanel("a", tags$head(tags$script(HTML('
                          var fakeClick = function(tabName, anchorName) {
                            var dropdownList = document.getElementsByTagName("a");
                            for (var i = 0; i < dropdownList.length; i++) {
                              var link = dropdownList[i];
                              if(link.getAttribute("data-value") == tabName) {
                                link.click();
                                console.log("Jumping to: ", anchorName);
                                setTimeout(() => {
                                document.querySelector("iframe").contentDocument.getElementById(anchorName).scrollIntoView({behavior: "smooth"});
                                }, 300);

                              };
                            }
                          };
        '))))
}

remove_y_axis_zero_tick_js <- function(multiomics_plot) {
  multiomics_plot %>% htmlwidgets::onRender("
function(el) {

  function hideZeroTicks() {
    // Select ALL y-axis tick labels for ALL subplots
    var allYTicks = el.querySelectorAll('.yaxislayer-above .ytick text, .yaxislayer-above .y2tick text, .yaxislayer-above .y3tick text, .yaxislayer-above .y4tick text');

    allYTicks.forEach(function(label) {
      if (label.textContent.trim() === '0') {
        label.style.display = 'none';   // hide zero tick
      } else {
        label.style.display = '';       // show all others
      }
    });
  }

  // Hide zero tick on initial draw
  hideZeroTicks();

  // Hide zero tick after any zoom/pan/resize/autoscale
  el.on('plotly_relayout', hideZeroTicks);
  el.on('plotly_afterplot', hideZeroTicks);
}
")
}


