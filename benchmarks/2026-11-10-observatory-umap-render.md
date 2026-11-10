# Observatory default UMAP: end-to-end benchmark

Tested on the local server with an already-started R app, Chromium headless 156, a fresh browser context per trial, one BLAS thread, the public human 2D embedding (4,522 samples, 177 tissue/cell-line groups), and a local connection. `paired_release.csv` records four interleaved before/after trials in each mode. Plot shown means the real Plotly widget has data, full layout and a plot container; this is not app-process startup time.

| Measurement | Before median | After median |
| --- | ---: | ---: |
| Open Observatory from an existing Browser session | 2.79 s | 1.29 s |
| Load a new session directly at Observatory | 3.58 s | 2.83 s |

A separate six-trial in-browser measurement from click to the registered Plotly selection hook was 1.03–1.13 s after optimization. The requested <1 second **total** was not met: first-load JavaScript parsing, Shiny initialization, and rendering 177 real legend groups remain. Numbers vary with CPU, browser, network, app configuration and cache warmth; they are not guarantees for ribocrypt.org.

## What cost time and what changed

- The original cache read was already about 10 ms. However Plotly rebuilt the cached R widget when rendering; a direct R microbenchmark measured 0.30–0.40 s for that rebuild. The new cache stores a completely built Plotly widget and a narrowly scoped S3 method skips the redundant rebuild. Plotly's normal per-session event registration still runs. `plotly_version`, source fingerprints and Run order still invalidate stale cache entries; access-filtered sessions do not use the shared cache.
- The plot had 177 native traces. WebGL reduces plotting work while preserving the real legend, hover and selection. Compact Run IDs and structured hover metadata reduced the plot WebSocket JSON from approximately 849 KB to 574 KB. All 4,522 coordinates, category labels, colors and hover fields were checked against the previous cache.
- The hidden Observatory browser now initializes only when its Browse tab opens. The sample table is built after the plot reports readiness; if the requested UMAP file is missing, the table is built immediately, and a three-second timer ensures a table remains available if plotting fails. The table-selection ready handler still restores selections.
- The Plotly library is loaded eagerly during initial app startup when collections exist. This moves its parse/execute cost earlier: median time to initial Browser Shiny connection increased from 0.70 to 1.39 s in this test. To favor the Browser's initial startup instead, set `options(ribocrypt.observatory_plotly_preload = "prefetch")` before constructing the app; in a separate paired run its Observatory-transition median was 1.95 s and initial connection 0.66 s. Use `"none"` to disable priming. Browser WebGL compatibility can be forced to the tested SVG path with `options(ribocrypt.umap_webgl = FALSE)`.

The cache builder regenerated 47 eligible species caches. Two historical collections lacking source provenance safely use the table-rendering fallback; 38 have no UMAP. Browser tests covered a real box drag (3,831 selected runs), plot/table selection synchronization, resets, native legend toggles, actual hover text, lazy Browse initialization, all 2D–6D choices and the SVG fallback. A missing 3D UMAP still allowed its sample table to load. Unit tests and an isolated package installation passed. One existing unrelated Observatory URL test fails identically on pre-change master (see prior `independent_dimensions/tests_master_baseline.log`).

To get reliably below one second **including a cold page load**, a larger change is needed: benchmark an alternative lightweight canvas/WebGL widget with a virtualized, searchable category legend and the same selection, hover, zoom and accessibility behavior, or a carefully hydrated pre-rendered view. A screenshot-only preview is insufficient because it would not yet be interactive. Compare cold and warm sessions, network conditions, memory and all existing interactions before replacing Plotly.

Raw scripts, browser timelines, paired measurements and test logs: `~/shared_workspace/hakon_playground/observatory_umap_20261109/performance_20261110/`.
