# Observatory UMAP below one second (in-session)

Follow-up to `2026-11-10-observatory-umap-render.md` (claude-megafst, 2026-11-10). Same server, headless Chromium 156,
public human 2D embedding (4,522 samples, 177 tissue/cell-line groups), `?go=FALSE` so the Browser's default gene plot
does not compete. Six interleaved paired trials per mode after one warm-up session per app; "shown" is the definition
of the earlier report (Plotly data, full layout and plot container present).

| Measurement (median, range) | d1c2da7 | this branch |
| --- | ---: | ---: |
| Open Observatory from an existing session: plot shown | 1.25 s (1.18–1.33) | **0.70 s (0.65–0.78)** |
| Same, legend complete | 1.57 s | 1.59 s |
| Load a new page directly at Observatory: plot shown | 3.57 s | 2.51 s |

Below one second is reached for opening Observatory within a session. A cold page load is not below one second:
about 1 s of it is downloading and parsing 1.5 MB of JavaScript (plotly.js 1.1 MB) before Shiny connects.

## Where the time went and what changed

Measured with Chrome traces, CPU profiles and an R profile of the server during the click.

- **Server, ~250 ms per first visit: creating the Observatory module.** `observatory_server()` and `collection_pages()`
  were created on the first click of the tab; Shiny's per-observer setup (labels, deparsing, OpenTelemetry hooks) over
  the many reactives took ~200–250 ms before any output could be computed. They are now created in the background
  1 s after the first page flush (`on_first_tab(prewarm=)`; `options(ribocrypt.observatory_prewarm = FALSE)` disables
  it, a number changes the delay). Hidden outputs stay suspended until the tab is shown, so the background cost is the
  module set-up only (~0.2–0.3 s server CPU per session). The click then reaches the plot message in ~230 ms instead of
  ~500–650 ms.
- **Browser, ~215 ms: the 177-entry legend.** `drawMarginPushers` (legend) was the largest part of `Plotly.newPlot`.
  The plot is now sent with `showlegend = FALSE` and `layout.meta.rcDeferredLegend`; `umap_plot_extension.js` switches
  the native legend on in the task after the first paint. Legend toggling, hover and selection are unchanged.
  The first render task fell from ~730 ms to ~370 ms (it also reduces the work of shinyjs's DOM mutation observer).
- **Payload 556 KB → 357 KB.** Coordinates are rounded to 3 decimals (Plotly serializes doubles with 50 digits;
  0.001 on a UMAP spanning tens of units is far below a pixel); the tissue | cell-line hover field is read from the
  trace name (`%{fullData.name}`) when it equals the group, as it always does for the default colouring (0/4,522
  mismatches), instead of per point; empty per-trace `error_x`/`error_y`/`line` from `plotly_build()` are dropped.
  The cache schema is 7 (`..._2d_v7.rds`): rebuild caches with `inst/scripts/precompute_observatory_umaps.R`.
- **Page HTML for app URLs.** The static-UI response cache (`cache_static_app_ui()`) bypassed every request with a
  query string, so links such as `?go=FALSE` or `?dff=..&gene=..` rendered the 268 KB page in R each time
  (0.6–0.85 s to first byte). The HTML is byte-identical for such queries, so they now share the cached page; Shiny's
  own restore/showcase parameters (`_state_id_`, `_inputs_`, `_values_`, `showcase`) and bookmarking still render per
  request. Apps with the Home page enabled (as ribocrypt.org) still do not cache the shell (see the comment in
  `RiboCrypt_app()`); their direct loads keep that cost.

WebGL note: this server has no GPU, so headless Chromium renders WebGL in software (an empty `scattergl` plot takes
~570 ms to paint here, an empty SVG plot ~60 ms). Absolute render times on users' machines with a GPU are lower; the
trace type was left as WebGL.

Checks: browser interactions (plot→Shiny selection, table→plot selection, reset, legend group toggle, lazy Browse tab),
hover text, real box drag (3,831 runs), R tests (cache, static UI, startup) and JS tests pass; the one failing JS test
(`Shiny updates tolerate document targets…`) fails identically on d1c2da7.
Scripts and raw numbers: `~/shared_workspace/hakon_playground/observatory_umap_20261109/sub1s_20261110/`
(`paired_sub1s.csv`, `timeline*.py`, `trace.py`, `cpuprof*.py`, `bench*.py`, `prof_nav.py`).
