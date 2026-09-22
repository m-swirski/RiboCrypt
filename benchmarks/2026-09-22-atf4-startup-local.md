# RiboCrypt ATF4 startup profile

Measured 2026-09-22 on clean local checkout
`fdea24ce4972382928ee8f8530c236e783e07d1e` (RiboCrypt 1.15.2).
No application code was changed for profiling.
Local Chrome headless, 1440 x 1000, localhost port 7819, one R process,
four sequential independent browser contexts (empty browser caches).
RiboCrypt and local ORFik loaded with devtools::load_all().
Data and metadata setup from ~/Desktop/run_ribocrypt.R; Browser selected,
plot_on_start TRUE, human_all_merged_l50, ATF4-ENSG00000128272,
ENST00000674920, RFP, translons enabled. Collection modules remained enabled.
All four sessions displayed the same 21-trace ATF4 plot without JavaScript errors.
The locally captured screenshot `browser-3.png` confirmed that the coverage
and sequence/gene tracks rendered.

## User-visible timings

| Session | Plot visible | Plot visible and Shiny idle |
|---|---:|---:|
| First user | 6.119 s | 6.256 s |
| Second user | 4.502 s | 4.800 s |
| Third user | 3.428 s | 3.715 s |
| Fourth user | 3.549 s | 3.841 s |

Warm-session median: 3.549 s to plot, 3.841 s to idle.
The immediately following second user was slower than the later warm users;
do not describe all subsequent loads as 3.8 s.

Completion checks: browser-c has Plotly fullData and main SVG, followed by
Shiny busy class clearing. Browser events, WebSocket messages, navigation
timings and screenshots were captured. There was a two-second observation
period after completion, excluded from reported load times.

## Warm critical path (fourth user)

| Sequential interval | Time |
|---|---:|
| Navigation to completed HTML response | 0.593 s |
| HTML response to Shiny connected (assets/client setup) | 0.587 s |
| Shiny connected to browser-c value received | 1.229 s |
| Plot value received to visible plot | 1.140 s |
| Visible plot to observed Shiny idle | 0.292 s |
| Total | 3.841 s |

The rendering interval includes widget processing and the browser completion
check, rather than an isolated measurement of Plotly.newPlot CPU time.

## R function timings

Fourth-session wall times, nested within the critical path above:

| Function/stage | Time |
|---|---:|
| Entire server registration | 0.633 s |
| browser_server | 0.103 s |
| browser_allsamp_server (hidden MegaBrowser) | 0.110 s |
| observatory_server (hidden) | 0.221 s |
| analysis_server (hidden) | 0.120 s |
| metadata_server (hidden) | 0.033 s |
| get_exp, four calls combined | 0.016 s |

Hidden module registration totals 0.484 s. This excludes their subsequent
observer execution. Traces show Heatmap, Codon and quality isoform updates
after the initial Browser output despite those tabs not being opened.
The second user's entire server registration took 1.484 s; the third and
fourth both took 0.633 s. Sampling contains compiler activity in session 2,
but the full difference has not been causally isolated.

First-session plot calculation costs:

- click_plot_browser_main_controller: 0.416 s.
- bottom_panel_shiny: 0.442 s.
- browser_track_panel_shiny: 0.203 s.

None of those three functions executed in sessions 2-4: app-level caching
already avoids this approximately 1.06-second calculation on identical loads.

## Before accepting users

- Launcher's experiment/metadata preparation: 1.946 s.
- RiboCrypt app construction: 3.223 s total, including parameter setup
  2.143 s and subsequent UI construction approximately 1.080 s.
- These exclude loading R/packages and are separate from page-navigation times.

## Payload

Each clean browser fetched 85 recorded subresources, about 1.78 MB transferred.
Plotly JS alone was approximately 1.14 MB transferred (3.58 MB decoded).
Initial HTML was 18.8 kB transferred (189.6 kB decoded).
WebSocket payloads totaled about 176 kB, including a 162 kB plot message.
On localhost, assets/client setup and plot rendering are appreciable costs;
actual remote users also have network latency and different device speeds.

## Interpretation

The measurements support investigating hidden module initialization, per-request
HTML handling, and client widget rendering. Moving experiment reads before
server would save little in this measured warm path. Plot construction is
already cached effectively. No optimization or repository edit was made.

This is a local instrumented baseline, not a production-site benchmark.
Function tracing, Rprof sampling and browser polling add some overhead.
Three warm observations characterize this run, not a production percentile.
Raw logs, browser JSON, Rprof files, screenshots and the temporary profiling
scripts were saved locally in `/tmp/ribocrypt-load-profile/`. They are not
included in this repository and may disappear when temporary files are cleaned.
Use explicit wall-time trace records for function timings; sampled Rprof
summaries include Shiny stack labels with embedded newlines and should not
be treated as precise additive stage totals.

## Hardware and software baseline

Environment information was collected on the same machine after the measurement,
using the same R executable and source checkouts. Memory availability and thread
settings below are follow-up observations, not telemetry recorded during the run.

| Component | Local baseline |
|---|---|
| CPU | Intel Core i7-8850H, x86_64 |
| CPU topology | 1 socket, 6 physical cores, 12 logical CPUs |
| CPU frequency | 2.60 GHz nominal; reported maximum 4.30 GHz |
| Memory | 32,532,192 KiB OS-visible RAM (31.03 GiB) |
| Swap | 2 GiB; no swap used when inspected after profiling |
| System/repository storage | Samsung PM981 NVMe 512 GB; root filesystem ext4 |
| BigWig/metadata storage | ST1000LX015-1U7172, 1 TB, rotational; NTFS mounted at /media/roler/S |
| Graphics hardware | Intel UHD Graphics 630 and NVIDIA Quadro P1000 Mobile |
| OS | Ubuntu 24.04.2 LTS |
| Kernel | Linux 7.0.0-31-generic |
| R | R Under development (unstable), 2024-03-18, revision r86148 |
| R platform | x86_64-pc-linux-gnu |
| R executable | /home/roler/R/x86_64-pc-linux-gnu-library/R-4.4-devel/bin/exec/R |
| BLAS | /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 |
| LAPACK | /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0 |
| Chrome | Google Chrome 146.0.7680.164, headless |
| Automation | Playwright 1.51.1; Node.js 18.19.1 |
| Browser viewport | 1440 x 1000 |
| Network | Loopback HTTP; browser and R process on the same machine |
| data.table threads | 6 in the follow-up environment check |
| Time zone | Europe/Oslo |
| Locale | LC_ALL unset; LC_CTYPE C.UTF-8; LC_COLLATE en_US.UTF-8; LC_NUMERIC C |

GPU hardware is listed for context; the renderer used by headless Chrome was
not recorded. CPU utilization, power mode, thermal throttling, peak RSS and disk
cache state were not measured or controlled. "First user" means first session
in this R process, not a cold operating-system file cache.

| R package | Version |
|---|---|
| RiboCrypt | 1.15.2 |
| ORFik | 1.29.13 |
| shiny | 1.11.1 |
| plotly | 4.12.0 |
| DT | 0.32 |
| data.table | 1.18.99 |
| htmlwidgets | 1.6.4 |
| httpuv | 1.6.16 |
| bslib | 0.9.0 |
| cachem | 1.1.0 |
| devtools | 2.4.5 |

Plotly's browser asset was version 2.25.2, distinct from the R package version.
ORFik was loaded from the clean local source checkout at
`e43cf261bb85a5698a71e8349a1974bd3fbc3d44`.

## Comparing with the server

Use the same RiboCrypt/ORFik commits, data and settings before attributing any
difference to hardware. The local configuration exposed 27 experiments,
3 collections and 13,845 metadata rows. Metadata included the extended QC CSV
and dominant-state annotations from the local launch script. Exact input-file
checksums were not captured, so experiment names alone do not guarantee identical
data on a second machine.

1. Launch one R process with Browser selected and automatic ATF4 plotting enabled.
   Load the source package with `devtools::load_all(".")`, using
   `env -u LC_ALL R --vanilla -q -e 'devtools::load_all("."); ...'`.
2. Match the measured settings: `default_experiment = "human_all_merged_l50"`,
   `default_experiment_meta = "all_samples-Homo_sapiens"`,
   `default_gene = "ATF4-ENSG00000128272"`,
   `default_gene_meta = "AMD1-ENSG00000123505"`,
   `plot_on_start = "TRUE"`, `default_view_mode = "tx"`,
   `translons = "TRUE"`, `collapsed_introns_width = "30"`,
   `full_annotation = "FALSE"`, `columns_plotly_backend = "gl"`.
   Keep collection and analysis modules enabled.
3. Open the first page, verify transcript ENST00000674920 and the RFP plot,
   then close that browser context. Keep the R process running.
4. Open three fresh browser contexts sequentially at the same viewport.
   Measure navigation to visible plot and to Shiny idle separately. Closing
   each context between trials matches the original run; concurrent users
   are a different workload.
5. Record the HTML response, Shiny connection, plot-value receipt, visible plot
   and idle timestamps, plus R function wall times. Keep startup time separate.
6. Record server hardware, R/package versions, storage, browser machine and
   network path. A remote browser comparison includes network costs that this
   localhost baseline does not. Keep local-browser and remote-browser results
   separately labeled.

The original profiler used function entry/exit traces and `Rprof(interval = 0.01)`
per session. Browser readiness checked for `browser-c._fullData` and a
`.main-svg`, followed by clearing of the document's `shiny-busy` class.
No application code changes are needed to reproduce those observations.
