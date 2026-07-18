# scE2G confidence interval workflow

**Status: implemented and tested** (February 2026)

## Overview

A Snakemake workflow (`workflow/Snakefile_ci`) that runs on top of a completed scE2G run. It generates confidence intervals on E2G scores by:
1. Subsampling cells from each cluster N times
2. Re-running the scE2G pipeline per subsample (two modes)
3. Annotating the full run's predictions with CI statistics from the subsamples

---

## Design decisions

### Mode A vs mode B

Two pipeline depth modes, both supported (run separately with different configs):
- **Mode A (`full`)**: Each subsample runs the complete pipeline: fragment → tagAlign → MACS2 peaks → ABC → features → model. Captures all sources of variability.
- **Mode B (`abc`)**: Each subsample creates a subsampled tagAlign, then reuses the full run's MACS2 peaks/candidate regions, re-runs ABC neighborhoods + predictions + all downstream features + model. ~50% faster per subsample; tests whether peak calling variability is a major contributor to CI width.

### Architecture: nested Snakemake

The CI Snakefile launches each subsample pipeline as an independent `snakemake` shell call (each with its own `--directory` and `--configfile`). This allows the existing `workflow/Snakefile` to run unchanged, avoids duplicating ABC rules, and cleanly separates CI orchestration from pipeline execution. `--nolock` is passed to nested calls since all subsample runs share the same working directory.

### Mode B: handling peak reuse

Before launching the subsample pipeline in mode B:
1. Create the subsample's `{results_dir}/{cluster}/Peaks/` directory
2. Copy full run's peak files there and `touch` them (update mtime to now)
3. Launch pipeline with `fragments_preprocessed: True`
4. Snakemake sees peaks exist and are newer than tagAlign → skips MACS2 → re-runs from neighborhoods onward

### Cell-level subsampling

Randomly select a fraction of cell barcodes from the fragment file. Filter the fragment file to those barcodes. Set `RNA_matrix_filtered: False` in subsample config so the pipeline auto-subsets the RNA matrix to matching barcodes.

### Adaptive subsample fraction

Default target: 70% of cells per subsample (creates meaningful variation while staying well above minimums).

Adaptive logic:
- Compute `n_cells_total` from the fragment file
- Compute `n_fragments_total` from the fragment file
- Compute `actual_fraction = min(1.0, max(target_fraction, min_cells/n_cells, min_fragments/n_fragments))`
- Warn if `actual_fraction >= 0.95` (subsamples will be near-identical, CIs will be artificially tight)

### E-G pair matching

Use `bedtools intersect -f 0.5 -r` (50% reciprocal overlap) on element coordinates, filtered by exact `TargetGene` match. If a full-run pair has no overlapping pair in a subsample, its score is 0.

---

## Configuration

See `config/config_ci_template.yaml` for a fully annotated template. Key parameters:

```yaml
full_run_results_dir: "/path/to/full/run/results"

scratch_dir: "/scratch/path"     # per-subsample predictions (large, temporary)
results_dir: "/path/to/ci/out"   # final annotated predictions and reports

clusters:
  - "K562_cluster1"
model: "multiome_powerlaw_v3"

mode: "full"          # "full" (mode A) or "abc" (mode B)

n_subsamples: 10
subsample_seed: 42
target_fraction: 0.70
min_cells: 100
min_fragments: 2000000

ci_level: 0.95
overlap_fraction: 0.50

snakemake_jobs: 50
```

---

## Files

```
workflow/
├── Snakefile_ci                             # Main CI workflow
└── rules/
│   ├── ci_subsample.smk                     # Subsample generation rules
│   ├── ci_run.smk                           # Launch subsample pipeline rules
│   └── ci_analysis.smk                      # CI computation rules
└── scripts/
    └── ci/
        ├── determine_subsample_params.py    # Compute adaptive fraction
        ├── subsample_cells.py               # Generate subsampled fragment file
        ├── generate_subsample_config.py     # Write per-subsample config.yaml
        ├── compute_ci_annotations.R         # CI statistics and annotation
        └── plot_ci.R                        # CI QC report (HTML)
config/
└── config_ci_template.yaml                  # Template CI config
tests/
└── config/
    └── test_config_ci.yaml                  # Test config (K562_cluster1_chr22p, 2 subsamples)
```

---

## Pipeline steps

### Step 1: Load full run configs (`Snakefile_ci` init)

Read from `{full_run_results_dir}/config/`:
- `scE2G_config.yml` → get `model_dir`, `gene_TSS500`, `chr_sizes`, `gene_annotations`, etc.
- `expanded_biosample_config.tsv` → get `atac_frag_file`, `rna_matrix_file`, `model_threshold`, `model_dir_base` per cluster
- Validate that requested `clusters` and `model` exist in the saved config
- Build `CI_BIOSAMPLE_DF` with per-cluster input paths and thresholds

### Step 2: Determine subsample parameters (per cluster)

**Rule `determine_subsample_params`** → `{scratch}/{cluster}/subsample_params.json`

Counts unique barcodes and total fragments; computes `actual_fraction`; writes JSON with `{n_cells, n_fragments, actual_fraction, n_cells_per_subsample}`.

### Step 3: Generate subsampled fragment files (per cluster × subsample)

**Rule `subsample_cells`** → `{scratch}/{cluster}/subsample_{i}/atac_fragments.tsv.gz` + `.tbi`

Randomly selects `round(n_cells × actual_fraction)` barcodes (seed = `subsample_seed + i`), filters fragment file, bgzips, and tabix-indexes.

### Step 4: Generate per-subsample configs

**Rule `generate_subsample_config`** → `{scratch}/{cluster}/subsample_{i}/config.yaml`

Copies full run's `scE2G_config.yml` and overrides:
- `results_dir` → `{scratch}/{cluster}/subsample_{i}/pipeline_output/`
- `cell_clusters` → new TSV pointing to subsampled fragment file
- `fragments_preprocessed: True` (correct key name per `workflow/Snakefile:28`)
- `RNA_matrix_filtered: False`
- `make_IGV_tracks: False`, `benchmark_performance: False`

### Step 5a (mode B only): Seed peaks from full run

**Rule `seed_peaks_mode_b`** → `{scratch}/{cluster}/subsample_{i}/.peaks_seeded`

Copies `{full_run_results}/{cluster}/Peaks/` into the subsample's pipeline output dir and touches all files to update timestamps, causing Snakemake to skip MACS2.

### Step 6: Launch subsample pipeline

**Rule `run_subsample_pipeline`** → `{scratch}/{cluster}/subsample_{i}/pipeline_output/{cluster}/{model}/encode_e2g_predictions.tsv.gz`

Launches the existing `workflow/Snakefile` as a nested `snakemake` call with `--nolock --rerun-incomplete`.

### Step 7: Compute CI annotations (per cluster)

**Rule `compute_ci_annotations`** → `{results}/{cluster}/{model}/encode_e2g_predictions_ci_annotated.tsv.gz` + `ci_summary_stats.tsv`

Annotates the **full (unthresholded)** predictions. For each subsample, uses `bedtools intersect -f {overlap_fraction} -r` to match elements, then joins by `TargetGene`. Builds a score matrix (rows = E-G pairs, cols = subsamples) and computes per-pair statistics:
- `E2G.Score.subsample_mean`
- `E2G.Score.subsample_sd`
- `E2G.Score.CI_low` / `E2G.Score.CI_high` (quantiles at `(1 ± ci_level) / 2`)
- `E2G.Score.pct_above_threshold`
- `n_subsamples_with_match`
- `E2G.Score.subsample_scores` (comma-separated individual scores, rounded to 6 decimals)

`ci_summary_stats.tsv` columns: `cluster`, `model`, `mode`, `n_subsamples`, `actual_fraction`, `n_cells`, `n_fragments`, `n_cells_per_subsample`, `n_pairs_full`, `n_pairs_above_threshold`, `median_CI_width`, `pct_pairs_CI_width_lt_0.1`, `median_pct_above_threshold`, `median_n_subsamples_matched`.

### Step 8: Generate CI QC report (per cluster)

**Rule `plot_ci`** → `{results}/{cluster}/ci_qc_report.html`

Four ggplot panels: CI width distribution, score stability scatter (full vs subsample mean), prediction stability by score bin, element coverage bar chart. Rendered to HTML via rmarkdown.

---

## Output directory structure

```
{scratch_dir}/
└── {cluster}/
    ├── subsample_params.json
    ├── subsample_0/
    │   ├── atac_fragments.tsv.gz
    │   ├── atac_fragments.tsv.gz.tbi
    │   ├── config.yaml
    │   ├── config_cell_clusters.tsv
    │   ├── .peaks_seeded               # mode B only
    │   └── pipeline_output/
    │       └── {cluster}/
    │           └── {model}/
    │               └── encode_e2g_predictions.tsv.gz
    └── subsample_1/ ...

{results_dir}/
└── {cluster}/
    ├── {model}/
    │   ├── encode_e2g_predictions_ci_annotated.tsv.gz
    │   └── ci_summary_stats.tsv
    └── ci_qc_report.html
```

---

## Known considerations

- **Nested Snakemake locking**: `--nolock` is required on nested calls because all subsample runs share the same working directory (`.snakemake/` lock would otherwise conflict). Each run gets its own `--directory` pointing to `{scratch}/{cluster}/subsample_{i}/pipeline_output/`.
- **Conda environments**: The conda environments are already shared across the pipeline; no `--conda-prefix` override is needed.
- **Small clusters**: If `actual_fraction` is forced to 1.0 because the cluster is too small, subsamples will be near-identical and CIs will be artificially tight. This is flagged with a warning in `determine_subsample_params.py` when `actual_fraction >= 0.95`.
- **Mode B correctness vs mode A**: Running both on the same clusters enables direct comparison of CI width to assess the contribution of peak calling variability.
- **Config key**: The correct config key for pre-processed fragments is `fragments_preprocessed` (not `preprocessed_fragments` — the latter may appear in stale `.snakemake/scripts/` temp files from older runs).

---

## Testing

Tested with `tests/config/test_config_ci.yaml` on `K562_cluster1_chr22p` (chr22p only; 480 cells, ~4.4k fragments) with `n_subsamples: 2`. All 10 Snakemake jobs completed successfully. Output verified:
- Both subsample pipelines ran end-to-end (fragments → predictions)
- Annotated predictions file contains all 6 CI columns appended to the standard 27-column schema
- Summary stats TSV written correctly
- HTML QC report generated (1.6 MB)

To run the test:
```bash
conda activate run_snakemake
snakemake -s workflow/Snakefile_ci \
  --configfile tests/config/test_config_ci.yaml \
  --use-conda -j 4
```
