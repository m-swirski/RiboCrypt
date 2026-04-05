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

  rendered_seq <- strsplit(as.character(target_seq),"")[[1]]
  translate_fuzzy_logic <- ifelse(all(rendered_seq %in% DNA_BASES), "error", "X")

  ir <- IRanges(seq(3), display_dist) # 3 AA frames
  aas <- suppressWarnings(translate(extractAt(target_seq[[1]], ir), if.fuzzy.codon = translate_fuzzy_logic))
  aas <- strsplit(as.character(aas), "")
  if (aa_letter_code == "three_letters") {
    aa_code <- AMINO_ACID_CODE
    aa_code["*"] <- "*"
    aas <- lapply(aas, function(x) aa_code[match(x, names(aa_code))] )
  }

  nts <- lapply(seq(3), function(x) seq(x, display_dist, 3))
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

addColumnsZoomSwitch <- function(multiomics_plot, threshold = columns_zoom_switch_threshold()) {
  multiomics_plot %>% htmlwidgets::onRender(
    "
function(el, x, data) {
  function getRange(ev) {
    function normalizeRange(lo, hi) {
      lo = Number(lo);
      hi = Number(hi);
      if (!isFinite(lo) || !isFinite(hi)) return null;
      return lo <= hi ? [lo, hi] : [hi, lo];
    }

    if (ev) {
      var keys = Object.keys(ev);
      for (var i = 0; i < keys.length; i++) {
        var key = keys[i];
        if (/^xaxis[0-9]*\\.range$/.test(key)) {
          var direct = ev[key];
          if (Array.isArray(direct) && direct.length === 2) {
            var normalized = normalizeRange(direct[0], direct[1]);
            if (normalized !== null) return normalized;
          }
        }
      }
      for (var j = 0; j < keys.length; j++) {
        var key0 = keys[j];
        var m = key0.match(/^(xaxis[0-9]*)\\.range\\[0\\]$/);
        if (m) {
          var base = m[1];
          var normalized2 = normalizeRange(
            ev[base + '.range[0]'],
            ev[base + '.range[1]']
          );
          if (normalized2 !== null) return normalized2;
        }
      }
    }

    var layout = el._fullLayout || {};
    var axisNames = Object.keys(layout).filter(function(k) {
      return /^xaxis[0-9]*$/.test(k);
    });
    for (var a = 0; a < axisNames.length; a++) {
      var axis = layout[axisNames[a]];
      if (axis && Array.isArray(axis.range) && axis.range.length === 2) {
        var normalized3 = normalizeRange(axis.range[0], axis.range[1]);
        if (normalized3 !== null) return normalized3;
      }
    }
    return null;
  }

  function buildSegments(trace, lo, hi) {
    var xs = Array.isArray(trace.x) ? trace.x : [];
    var ys = Array.isArray(trace.y) ? trace.y : [];
    var frame = (trace.meta && trace.meta.rc_columns_frame) || trace.name || '';
    var xOut = [];
    var yOut = [];
    var textOut = [];

    for (var i = 0; i < xs.length && i < ys.length; i++) {
      var xVal = Number(xs[i]);
      var yVal = Number(ys[i]);
      if (!isFinite(xVal) || !isFinite(yVal)) continue;
      if (xVal < lo || xVal > hi) continue;
      xOut.push(xVal, xVal, null);
      yOut.push(0, yVal, null);
      var label = 'position: ' + xVal + '<br>count: ' + yVal + '<br>frame: ' + frame;
      textOut.push(label, label, null);
    }

    if (!xOut.length) {
      xOut = [null];
      yOut = [null];
      textOut = [null];
    }

    return {x: xOut, y: yOut, text: textOut};
  }

  function getPanelWidthPx() {
    var layout = el._fullLayout || {};
    var candidates = ['xaxis', 'xaxis2', 'xaxis3', 'xaxis4'];
    for (var i = 0; i < candidates.length; i++) {
      var axis = layout[candidates[i]];
      if (axis && axis._length && isFinite(axis._length)) return axis._length;
    }
    if (el.clientWidth && isFinite(el.clientWidth)) return el.clientWidth;
    return 800;
  }

  function getGlLineWidth(span) {
    if (!isFinite(span) || span <= 0) return 6;
    var panelWidth = getPanelWidthPx();
    var pxPerNt = panelWidth / span;
    return Math.max(1, Math.min(12, Math.round(pxPerNt * 0.8 * 10) / 10));
  }

  function updateColumns(ev) {
    var dataTraces = (el.data || []);
    var groups = {};
    for (var i = 0; i < dataTraces.length; i++) {
      var meta = dataTraces[i].meta;
      if (!meta || !meta.rc_columns_switch) continue;
      var key = (meta.rc_columns_group || '') + '::' + (meta.rc_columns_frame || '');
      if (!groups[key]) groups[key] = {};
      if (meta.rc_columns_switch === 'gl_subset') groups[key].gl = i;
      if (meta.rc_columns_switch === 'line') groups[key].line = i;
    }
    var groupKeys = Object.keys(groups).filter(function(key) {
      return groups[key].gl !== undefined && groups[key].line !== undefined;
    });
    if (!groupKeys.length) return;

    var range = getRange(ev);
    if (range === null) return;
    var lo = range[0];
    var hi = range[1];
    var span = Math.abs(hi - lo);
    var showGl = span < data.threshold;
    var rangeKey = lo + ':' + hi;
    var glLineWidth = getGlLineWidth(span);
    if (el.__rcColumnsMode === (showGl ? 'gl' : 'line') &&
        el.__rcColumnsRangeKey === rangeKey &&
        (!showGl || el.__rcColumnsWidth === glLineWidth)) {
      return;
    }

    if (showGl) {
      var glIdx = [];
      var glX = [];
      var glY = [];
      var glText = [];
      var glVisible = [];
      var glWidth = [];
      var lineIdx = [];
      var lineVisible = [];
      for (var k = 0; k < groupKeys.length; k++) {
        var pair = groups[groupKeys[k]];
        var glTraceIndex = pair.gl;
        var lineTraceIndex = pair.line;
        var seg = buildSegments(dataTraces[lineTraceIndex], lo, hi);
        glIdx.push(glTraceIndex);
        glX.push(seg.x);
        glY.push(seg.y);
        glText.push(seg.text);
        glVisible.push(true);
        glWidth.push(glLineWidth);
        lineIdx.push(lineTraceIndex);
        lineVisible.push(false);
      }
      Plotly.restyle(el, {
        x: glX,
        y: glY,
        text: glText,
        visible: glVisible,
        hoverinfo: glVisible.map(function() { return 'text'; }),
        hovertemplate: glVisible.map(function() { return '%{text}<extra></extra>'; }),
        'line.width': glWidth
      }, glIdx);
      Plotly.restyle(el, {visible: lineVisible}, lineIdx);
    } else {
      var glIdx2 = [];
      var glVisible2 = [];
      var lineIdx2 = [];
      var lineVisible2 = [];
      for (var k2 = 0; k2 < groupKeys.length; k2++) {
        var pair2 = groups[groupKeys[k2]];
        glIdx2.push(pair2.gl);
        glVisible2.push('legendonly');
        lineIdx2.push(pair2.line);
        lineVisible2.push(true);
      }
      Plotly.restyle(el, {visible: glVisible2}, glIdx2);
      Plotly.restyle(el, {visible: lineVisible2}, lineIdx2);
    }
    el.__rcColumnsMode = showGl ? 'gl' : 'line';
    el.__rcColumnsRangeKey = rangeKey;
    el.__rcColumnsWidth = glLineWidth;
  }

  el.on('plotly_relayout', updateColumns);
  el.on('plotly_afterplot', updateColumns);
  updateColumns();
}
",
    list(threshold = threshold)
  )
}

addBrowserXRangeClamp <- function(multiomics_plot, min_x = 1, max_x) {
  multiomics_plot %>% htmlwidgets::onRender(
    "
function(el, x, data) {
  function normalizeRange(lo, hi) {
    lo = Number(lo);
    hi = Number(hi);
    if (!isFinite(lo) || !isFinite(hi)) return null;
    return lo <= hi ? [lo, hi] : [hi, lo];
  }

  function clampRange(range) {
    var parsed = normalizeRange(range[0], range[1]);
    if (!parsed) return null;
    var lo = Math.max(data.min_x, parsed[0]);
    var hi = Math.min(data.max_x, parsed[1]);
    if (hi < lo) {
      hi = lo;
    }
    return [lo, hi];
  }

  function sameRange(a, b) {
    return Array.isArray(a) && Array.isArray(b) &&
      a.length === 2 && b.length === 2 &&
      Math.abs(Number(a[0]) - Number(b[0])) < 1e-9 &&
      Math.abs(Number(a[1]) - Number(b[1])) < 1e-9;
  }

  function collectRanges(ev) {
    var updates = {};
    var layout = el._fullLayout || {};
    var axisNames = Object.keys(layout).filter(function(k) {
      return /^xaxis[0-9]*$/.test(k);
    });

    function maybeClamp(axisName, range) {
      if (!Array.isArray(range) || range.length !== 2) return;
      var clamped = clampRange(range);
      if (!clamped || sameRange(range, clamped)) return;
      updates[axisName + '.range'] = clamped;
      updates[axisName + '.autorange'] = false;
    }

    if (ev) {
      Object.keys(ev).forEach(function(key) {
        if (/^xaxis[0-9]*\\.range$/.test(key) && Array.isArray(ev[key])) {
          var axisName = key.replace(/\\.range$/, '');
          maybeClamp(axisName, ev[key]);
        }
      });

      axisNames.forEach(function(axisName) {
        var loKey = axisName + '.range[0]';
        var hiKey = axisName + '.range[1]';
        if (loKey in ev || hiKey in ev) {
          var range = normalizeRange(ev[loKey], ev[hiKey]);
          if (range) maybeClamp(axisName, range);
        } else if (ev[axisName + '.autorange']) {
          var axis = layout[axisName];
          if (axis && Array.isArray(axis.range)) maybeClamp(axisName, axis.range);
        }
      });
    } else {
      axisNames.forEach(function(axisName) {
        var axis = layout[axisName];
        if (axis && Array.isArray(axis.range)) maybeClamp(axisName, axis.range);
      });
    }

    return updates;
  }

  function enforceClamp(ev) {
    if (el.__rcClampBusy) return;
    var updates = collectRanges(ev);
    if (!updates || !Object.keys(updates).length) return;
    el.__rcClampBusy = true;
    Plotly.relayout(el, updates).then(function() {
      el.__rcClampBusy = false;
    }).catch(function() {
      el.__rcClampBusy = false;
    });
  }

  el.on('plotly_relayout', enforceClamp);
  el.on('plotly_afterplot', function() { enforceClamp(null); });
  enforceClamp(null);
}
",
    list(min_x = min_x, max_x = max_x)
  )
}

addMegabrowserDoubleClickReset <- function(plot_object, reset_range, peer_ids = character(),
                                           reset_layout = NULL,
                                           peer_reset_layout = NULL) {
  if (is.null(reset_layout)) {
    reset_layout <- list(
      "xaxis.range" = reset_range,
      "xaxis.autorange" = FALSE
    )
  }
  if (is.null(peer_reset_layout)) {
    peer_reset_layout <- list(
      "xaxis.range" = reset_range,
      "xaxis.autorange" = FALSE
    )
  }
  plot_object %>% htmlwidgets::onRender(
    "
function(el, x, data) {
  function cloneRange(range) {
    return Array.isArray(range) ? [range[0], range[1]] : null;
  }

  function suppressEvent(evt) {
    if (!evt) return;
    if (typeof evt.preventDefault === 'function') evt.preventDefault();
    if (typeof evt.stopPropagation === 'function') evt.stopPropagation();
    if (typeof evt.stopImmediatePropagation === 'function') evt.stopImmediatePropagation();
  }

  function captureInitialResetLayout(target) {
    if (!target || target.__rcMegabrowserInitialResetLayout) return;
    var layout = target._fullLayout || {};
    var resetLayout = {};

    if (layout.xaxis && Array.isArray(layout.xaxis.range)) {
      resetLayout['xaxis.range'] = cloneRange(layout.xaxis.range);
      resetLayout['xaxis.autorange'] = false;
    }
    if (layout.yaxis && Array.isArray(layout.yaxis.range)) {
      resetLayout['yaxis.range'] = cloneRange(layout.yaxis.range);
      resetLayout['yaxis.autorange'] = false;
    }

    if (Object.keys(resetLayout).length) {
      target.__rcMegabrowserInitialResetLayout = resetLayout;
    }
  }

  function buildResetLayout(layout) {
    var cloned = {};
    var source = layout || {};

    Object.keys(source).forEach(function(key) {
      var value = source[key];
      cloned[key] = Array.isArray(value) ? value.slice() : value;
    });

    return cloned;
  }

  function applyReset(target, layout) {
    if (!target || typeof Plotly === 'undefined') return;
    captureInitialResetLayout(target);
    Plotly.relayout(target, buildResetLayout(layout));
  }

  function triggerReset() {
    if (el.__rcMegabrowserResetBusy) return false;
    el.__rcMegabrowserResetBusy = true;
    setTimeout(function() {
      applyReset(el, data.reset_layout);
      (data.peer_ids || []).forEach(function(id) {
        var peer = document.getElementById(id);
        if (peer) applyReset(peer, data.peer_reset_layout);
      });
      setTimeout(function() {
        el.__rcMegabrowserResetBusy = false;
      }, 0);
    }, 0);
    return false;
  }

  el.on('plotly_doubleclick', function(evt) {
    suppressEvent(evt);
    if (evt && evt.event) suppressEvent(evt.event);
    return triggerReset();
  });

  el.addEventListener('dblclick', function(evt) {
    suppressEvent(evt);
    return triggerReset();
  }, true);

  var attachInnerDblclick = function() {
    var inner = el.querySelectorAll('.gl-container, .nsewdrag, .plotly .user-select-none, canvas');
    Array.prototype.forEach.call(inner, function(node) {
      if (!node || node.__rcMegabrowserDblclickBound) return;
      node.__rcMegabrowserDblclickBound = true;
      node.addEventListener('dblclick', function(evt) {
        suppressEvent(evt);
        return triggerReset();
      }, true);
    });
  };

  attachInnerDblclick();
  captureInitialResetLayout(el);
  el.on('plotly_afterplot', function() { captureInitialResetLayout(el); });
  el.on('plotly_afterplot', attachInnerDblclick);
}
",
    list(
      reset_range = reset_range,
      reset_layout = reset_layout,
      peer_reset_layout = peer_reset_layout,
      peer_ids = as.list(peer_ids)
    )
  )
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
