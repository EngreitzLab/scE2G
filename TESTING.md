# Testing

This document describes what's tested in scE2G, how, and what's known to be untested. See
`CLAUDE.md` for project overview and `TODO.md` for the backlog (including test-coverage gaps
tracked as TODO items, cross-referenced below).

## Overview: three test tiers, always run

CircleCI (`.circleci/config.yml`, a single static job) always runs, in order:

1. **Build rule conda envs** — `snakemake --use-conda --conda-create-envs-only`, isolates env
   solve/build failures into their own labeled step. Cached (see below).
2. **Run unit tests** (tier 1) — fast, no full pipeline run.
3. **Run tests** (tiers 2 and 3) — the full end-to-end pipeline on chr22 test data, which also
   exercises the tier-2 stage checkpoints as part of the same run.

Every push, on every branch, runs all three steps unconditionally — there is no path-based
skipping. This is deliberate: CircleCI has no reliable, documented way to detect "this push is
part of an open PR" from `config.yml`, and diff-based skipping has a real failure mode (it diffs
each push against its own previous commit, not against the PR's target branch, so a PR's *final*
push could go green without ever validating the PR's full cumulative diff). Given the full E2E
test only takes a few minutes on the chr22 test dataset, always running it is cheaper than the
risk of a silent gap.

`ABC/` and `ENCODE_rE2G/` are git submodules with their **own** CircleCI configs — MACS2 peak
calling, ABC contact prediction, and genomewide feature assembly are tested independently
upstream. scE2G's own tests focus on the code it owns: `workflow/rules/*.smk` and
`workflow/scripts/{feature_computation,model_application,prediction_qc}/*`.

## Env-build caching

The "Build rule conda envs" step caches `.snakemake/conda` (`restore_cache`/`save_cache`), keyed
on a checksum of the three env files that feed it: `workflow/envs/sc_e2g.yml`,
`ENCODE_rE2G/workflow/envs/encode_re2g.yml`, `ENCODE_rE2G/ABC/workflow/envs/abcenv.yml`. Unchanged
env files restore instead of re-solving from scratch; changing any of them naturally busts the
cache key, so a real env change always gets a fresh, fully-tested solve.

## Nightly scheduled builds

A `nightly` pipeline parameter skips the cache restore entirely, forcing a from-scratch env solve
even when the env files haven't changed:

- Set via a CircleCI **Scheduled Pipeline** trigger (Project Settings > Triggers) passing custom
  parameter `nightly: true` — configure separately for `main` and `dev`.
- This is the only thing nightly changes — unit tests and the full E2E test already run on every
  push regardless, so nightly's whole purpose is to periodically re-verify the env solve itself
  (catching upstream conda/bioconda/pip drift, e.g. a dependency publishing a breaking version
  bump) without needing to wait for someone to touch an env file.
- `save_cache` still runs afterward on a nightly build, refreshing the cache with a
  freshly-verified solve for the next regular push.

## Tier 1: unit tests

`workflow/rules/unit_tests.smk` adds `test_python_unit`/`test_r_unit` Snakemake rules that run
`pytest tests/unit/python` / `testthat::test_dir("tests/unit/r")` through Snakemake's own
`--use-conda` env resolution for `sc_e2g.yml` (no second conda env to build). Real tests exist for
`training_functions.py`, `run_e2g_cv.py` (Python), and `get_fill_values.R` (R).

Gotchas:
- Snakemake's bash strict mode treats pytest's exit code 5 ("no tests collected") and testthat's
  "No test files found" error as rule failures — both rules guard against exactly those two
  conditions (and only those) so real test failures still fail the build.
- `testthat::test_dir()` runs each test file with its own directory as the working directory, not
  the repo root — R test files `source()` the script under test via a relative `../../../workflow/...`
  path, not a repo-root-relative one.

**Follow-up refactor pattern (not yet applied)**: `get_fill_values.R` needed no changes since it
was already a plain sourceable function file. The more monolithic rule scripts
(`compute_kendall.R`, `compute_arc_e2g.R`, `merge_features_with_crispr_data_apply.R`) interleave
pure logic with `snakemake@input`-style script-level I/O and aren't unit-testable as-is. When
tackling those, wrap the script-level I/O in `if (exists("snakemake")) { ... }` so the file stays
sourceable (for tests) while still working as a Snakemake `script:` when `snakemake` the object is
present.

## Tiers 2 and 3: stage-level checkpoints and the full E2E test

`tests/test_sce2g_apply.py` checks, in addition to the final thresholded predictions file:

- `Kendall/Pairs.Kendall.tsv.gz` — small enough (~500KB) to compare as a full fixture.
- `ARC/EnhancerPredictionsAllPutative_ARC.tsv.gz` and pre-threshold
  `multiome_powerlaw_v3/scE2G_predictions.tsv.gz` — both still several MB even after chr22-scoping
  (one row per candidate pair genome-wide), so these are checked via `hash_large_file()` in
  `tests/utils.py`: row count + MD5 of identity columns and score columns rounded to 6 decimals,
  immune to cross-node BLAS floating-point noise while still catching real regressions.
- Distal (non-promoter) element scores on the pre-threshold file — the final thresholded output on
  this chr22 test dataset is, and always has been, 100% "promoter"-class rows (distal elements get
  real, varying, non-degenerate scores, just genuinely below the 0.177 threshold at this data
  scale). `check_distal_elements_have_nonzero_scores()` guards against distal scoring silently
  going degenerate (all-zero/NaN), which the thresholded-only check can never catch since no
  distal row is ever in it.

The test's `genes`/`gene_TSS500` references are scoped to chr22 only (via the `alt_genes`/`alt_TSS`
per-biosample config columns in `tests/config/test_config_cell_clusters.tsv`) rather than
genome-wide — gene-enhancer pairing is strictly same-chromosome, so a genome-wide gene/TSS
reference against chr22-only test data was inflating intermediate files (pre-threshold
predictions, ARC, genomewide features, candidate regions) with rows that could never produce
output pairs anyway. `alt_genes` restricts pairing; `alt_TSS` separately restricts the
TSS/promoter reference merged into candidate-region generation
(`ABC/workflow/scripts/peaks.py`'s `--regions_includelist`) — the two are driven by different
config keys and both need scoping.

To regenerate fixtures after an intentional pipeline change: run the full E2E test once to produce
real output under `tests/test_output/`, then `bash tests/replace_expected_output.sh
<results_dir>`.

## Known test-coverage gaps (not prioritized/scheduled)

- The `alt_TSS`/`alt_genes` fallback branch (no override → falls back to
  `config['ref']['genome_tss']`/`config['ref']['genes']`) has no test coverage — the test dataset
  always uses the override branch. Low priority: it's a trivial 2-line conditional in ABC's own
  `_configure_tss_and_gene_files()` (`ABC/workflow/rules/utils.smk`), not scE2G's code.
- `crispr_benchmarking: True` (CV model application, `benchmark_performance.py`) — never exercised
  (test config sets `benchmark_performance: False`).
- A populated `HiC_file` (HiC-based ABC scoring instead of the powerlaw model) — never exercised.
- `fragments_preprocessed: True` (pre-tabixed fragment input path) — never exercised.
- A hard rank-cutoff tie-break in `ABC/workflow/scripts/peaks.py`'s MACS2 peak ranking
  (`sort -nr -k 4 | head -n {nStrongestPeaks}`) is nondeterministic when many peaks tie at the
  cutoff boundary — invisible on chr22-scoped test data since candidate region counts here never
  approach the default `nStrongestPeaks=175,000`. This is a scale-dependent bug class (memory
  pressure, integer overflow, sort tie-breaking, performance regressions) that's structurally
  invisible to fast chr22-subset tests — see `TODO.md` item 1a. Blocked on patching the upstream
  `ABC` submodule (external org, not EngreitzLab-owned).
