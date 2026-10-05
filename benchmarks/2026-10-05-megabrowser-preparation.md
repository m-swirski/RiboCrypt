# MegaBrowser preparation and group inspection

Date: 2026-10-05. Working tree on master, based on `ba01e05` and the preceding
MegaBrowser iterations. Not committed or pushed.

## Changes

- Plotly matrix preparation now applies the saved order once instead of reversing
  it twice. Full precision, row order, renderer selection and plot ranges are unchanged.
- Display metadata no longer eagerly computes enrichment. The result-filter reset
  observer watches analysis keys, not computed statistics. Enrichment remains cached
  and computed from individual libraries when an enrichment/statistics/result output
  needs it. Existing non-server callers retain eager statistics by default.
- Group summary tab reports per-group library counts, categorical modes and agreement
  percentages, missing counts, and numeric means/minima/maxima. CSV exports provide
  summaries and individual library memberships for currently visible groups.
- Summaries are independent of collapse. Enrichment and the original result table
  remain full-library analyses. Summary downloads retain numeric precision.

## Benchmark

Same i7-8850H machine (6 cores / 12 threads, 31 GiB RAM) as preceding reports.
R devel 2024-03-18 r86148; RiboCrypt source-loaded with `devtools::load_all()`.
Deterministic synthetic matrices with five saved groups; no data loading or clustering
inside timing. First call warms each function; each trial runs ten preparations,
with garbage collection before, not inside, the timed trial. Five trials per case.
The baseline reproduces the preceding preparation function. Exact complete returned
lists are compared with `identical()` before timing. These are **backend timings**,
not app refresh times or measurements of real AMD1 coverage.

| Positions x libraries | Before (ms/call) | After (ms/call) | Median change |
|---|---|---|---|
| 483 x 1,754 (AMD1-sized) | 13.9, 14.1, 14.3, 14.4, 14.0 | 12.0, 11.8, 11.6, 11.7, 11.9 | 14.1 to 11.8, 16% faster |
| 2,000 x 3,000 | 118.1, 119.5, 119.6, 120.7, 123.4 | 86.1, 85.7, 131.3, 86.8, 89.0 | 119.6 to 86.8, 27% faster |

One updated large-matrix trial is slower than baseline; memory allocation/GC introduces
variance. This removes one full matrix copy but does not reduce the serialized heatmap
payload. The previous report's ~8.55 MB full payload and ~300 ms drawing cost remain
larger opportunities. No unmeasured end-to-end speedup is claimed for lazy statistics.

Reproduce from repo root:

```bash
env -u LC_ALL R --vanilla -q -e 'devtools::load_all("."); source("benchmarks/megabrowser-preparation.R")'
```

Machine-readable output: `/tmp/megabrowser-preparation.json`.

## Verification and limits

Tests cover exact matrix/data.table row order, deferred/eager statistics equivalence,
full-library analysis invariants, visible-group membership, categorical agreement with
missing values, numeric summaries, CSV precision and immutable original metadata.
Tutorial Rmd and served HTML are updated.

Final full suite: **3,091 assertions**, zero failures, warnings or skips.
Focused MegaBrowser suite: **357 assertions**, all passing. `git diff --check` clean.

The full real-data app cannot be launched under the current filesystem restrictions:
ORFik's database lock lies outside writable roots. Live warm-refresh comparisons against
the preceding app are therefore pending; synthetic preparation timings must not be
substituted for browser measurements.

An isolated Shiny fixture and Chrome check were also attempted, but this environment
blocks TCP server creation and Chrome socket operations. No browser verification is
claimed for this iteration.
