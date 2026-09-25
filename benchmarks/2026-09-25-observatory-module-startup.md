# Defer the unused collection module

Date: 2026-09-25. Branch: `perf/observatory-startup`.
Control: master `7166831`. Candidate: independent first-use initialization of
MegaBrowser and Observatory, retaining their shared collection state.

## Result

Keep this small change, but it does not solve the whole startup problem.
Warm median time to both populated widgets improved from **3.034 to 2.802 s**
(232 ms / 7.6%). Time to the same readiness check plus Shiny idle improved from
**3.034 to 2.887 s** (147 ms / 4.9%). One pair was slower for the candidate;
the remaining five favored it. More work is needed for a substantially faster page.

The measured hidden MegaBrowser data requests fell from **seven to zero** on
every candidate startup. First-use cost moves to opening MegaBrowser; it is not
eliminated globally. Revisiting either tab preserves initialized modules and state.

## Why this candidate

The launcher now has `plot_on_start=TRUE` and `init_tab_focus="Observatory"`.
Therefore earlier Plotly loading, already in master, is active. The previous
[3.64-second live-server profile](2026-09-24-observatory-startup.md) used a launcher
with FALSE and initial MegaBrowser focus. It is not a valid unchanged control
for today's module change. Do not credit the full historical-to-current difference
to this patch.

The existing outer first-tab callback still prepares collection data, URL-derived
defaults and shared reactive state once per session. Inside that closure, each
module now has its own `on_first_tab()` callback instead of both being started
immediately. Both callbacks retain the same shared bindings. No data preparation,
URL parser, scientific settings, UMAP renderer or DT selection logic was changed.
Observatory Browse initialization remains unchanged in this iteration.

## Measurement

Two source-loaded processes on 7820 (control) and 7821 (candidate); the control
was fully started before editing the R function. Both use the same local launcher
metadata, human merged Browser default, human collection and default gene settings,
automatic plotting enabled and initial Observatory focus. Both RiboCrypt and local
ORFik loaded with `devtools::load_all()` using `env -u LC_ALL R --vanilla -q -e`.
Hardware/software baseline is documented in the earlier
[local startup report](2026-09-22-atf4-startup-local.md#hardware-and-software-baseline).
The user's old port 4723 was not listening; no user server was modified or stopped.

Each app received three primer visits. Then one headless Chrome process alternated
seven control/candidate rounds, using a fresh empty-cache 1440 x 1000 context
for every page. Round 0 is excluded as browser warmup. No R tests or other browser
jobs ran during the comparison. Completion criteria are unchanged from the
Observatory profile; each visit verified 3,857 total libraries, 15 first-page rows,
populated visible UMAP/DT and no JavaScript or Shiny errors.

Milliseconds from navigation:

| Round | Control both ready | Candidate both ready | Control idle | Candidate idle |
|---|---:|---:|---:|---:|
| 0, excluded | 3090.7 | 3004.8 | 3201.3 | 3083.7 |
| 1 | 2857.4 | 2669.5 | 2933.3 | 2745.6 |
| 2 | 3021.0 | 2703.6 | 3021.0 | 2794.8 |
| 3 | 3206.9 | 2805.9 | 3206.9 | 2889.9 |
| 4 | 3047.8 | 2858.6 | 3047.8 | 2931.3 |
| 5 | 2846.2 | 2798.7 | 2923.1 | 2884.5 |
| 6 | 3847.9 | 3883.2 | 3848.0 | 3929.1 |
| Median, 1-6 | 3034.4 | 2802.3 | 3034.4 | 2887.2 |

The last round was slower for both apps and is retained. These are small local
samples with fixed order, not production percentiles. Readiness polling cannot
run during main-thread blocking and is not an exact browser paint timestamp.

## Validation

The full R suite passed 2,805 assertions with no failures, warnings or skips.
Added tests exercise independent nested initialization in both entry orders,
shared reactive updates and state retention without repeated module registration.
The real-browser URL smoke passed restored groups/settings, highlighting,
Get URL round trip, page/selected-row subsets, Remove subset and `go=FALSE`.
Screenshot inspection confirmed the populated UMAP and first 15 library rows.
The owned benchmark servers were stopped afterward.

## Reproduction

The Observatory harness now accepts comma-separated URLs for alternating trials:

```sh
RESULTS=/tmp/obs-isolated-comparison.json \
  NODE_PATH=/tmp/observatory-url-browser/node_modules \
  node benchmarks/observatory-startup.cjs \
  'http://127.0.0.1:7820/#Observatory,http://127.0.0.1:7821/#Observatory' 7
NODE_PATH=/tmp/observatory-url-browser/node_modules \
  node benchmarks/collection-startup-smoke.cjs http://127.0.0.1:7821/
```

The new smoke test must fail against the old control: it sees seven hidden-module
requests where zero are expected. It passes against the candidate in both tab-entry
orders, including first-use initialization, subset preservation across tab switches,
Remove subset and retained MegaBrowser gene/transcript choices. It does not generate
a full MegaBrowser coverage heatmap.
