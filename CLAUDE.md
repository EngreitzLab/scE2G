# scE2G - Single-cell Enhancer-to-Gene

A Snakemake pipeline for predicting genome-wide enhancer-gene regulatory links from single-cell ATAC-seq and/or multiome data. Developed by Engreitz and Andersson Labs at Stanford University.

**Publication:** bioRxiv 2024.11.23.624931v1

## Session handover

**Keep `HANDOVER.md` up to date as you work.** It exists so a new Claude Code session can pick
up seamlessly if the current one is killed or restarted. Whenever you complete a meaningful
step (finish an investigation, commit something, kick off a long-running job, discover a
non-obvious gotcha), update `HANDOVER.md` with: current branch/commit, what's done, what's
in progress, uncommitted files and why, any running/pending jobs, and next steps. Read
`HANDOVER.md` (along with this file and `TODO.md`) at the start of a new session before doing
anything else.

## Branch overview

- **`main`** — stable/production line. Submodule architecture: `ENCODE_rE2G` (submodule) nests its own `ABC` submodule. Has CircleCI testing (unit tests + chr22-scoped E2E, see `TESTING.md`).
- **`dev`** (this branch) — integration branch for staging larger changes before they land on `main`. As of 2026-08-01, its testing setup was reset to match `main`'s (see `TESTING.md`) after an earlier, more elaborate dev-only CI design (path-based stage skipping, env-build caching, nightly rebuilds) was set aside in favor of the simpler single-static-job approach that actually landed on `main`. This branch otherwise carries project docs (`CLAUDE.md`/`TODO.md`/`HANDOVER.md`) and other non-testing work not yet on `main`.
- **`flatten-encode-re2g-submodule`** — flattens `ENCODE_rE2G` directly into scE2G, leaving `ABC` as the only (direct) submodule; also consolidates fragment-processing rules and removes Snakemake checkpoint indirection in favor of upfront config evaluation. Design notes in `docs/plan_flatten_encode_re2g.md` (on that branch only). Not yet merged into `dev`.
- **`confidence-intervals-cell-type-comparison`** — branched from `flatten-encode-re2g-submodule` (its `Snakefile_ci` calls `ABC.enable_retry(...)` directly, so it depends on the flattened namespace). Adds the CI subsampling workflow and cross-cell-type comparison scripts described later in this doc — those sections describe code that lives on **this branch**, not on `dev`/`main`. Not yet merged.

## Quick start

```bash
# Setup
conda env create -f workflow/envs/run_snakemake.yml
conda activate run_snakemake

# Run pipeline
snakemake -j1 --use-conda --configfile config/config.yaml

# Run tests
pytest -s tests/test_sce2g_apply.py
```

## Project structure

```
config/                 # Configuration files (config.yaml, config_training.yaml)
workflow/
├── Snakefile           # Main workflow (imports ENCODE_rE2G as a Snakemake module)
├── Snakefile_training  # Model training workflow
├── rules/              # Modular Snakemake rules
├── scripts/            # Analysis scripts (R, Python)
└── envs/               # Conda environments
models/                 # Pre-trained logistic regression models
ENCODE_rE2G/            # Git submodule; nests its own ABC submodule inside it
resources/              # Reference data, annotations, example data
tests/                  # Unit tests + end-to-end test with expected outputs
```

(On `flatten-encode-re2g-submodule` and later, `ENCODE_rE2G/` is gone and `ABC/` is a direct top-level submodule instead.)

## Technologies

- **Workflow:** Snakemake 7+
- **Languages:** Python 3.6+, R (Tidyverse/Seurat/Signac), Bash
- **ML:** scikit-learn 1.2.1 (logistic regression)
- **Genomics:** bedtools, tabix, samtools, MACS2

## Pipeline steps

1. Fragment processing → tagAlign conversion
2. MACS2 peak calling
3. ABC (Accessibility-Chromatin-Based) contact prediction via ABC submodule
4. Kendall correlation between ATAC and RNA (multiome only)
5. ARC-E2G scoring (multiome only)
6. Feature generation and integration
7. Logistic regression model application
8. Quantile normalization against reference distribution

## Key inputs

- Pseudobulk ATAC fragments (bgzip + tabix indexed, 5-column TSV)
- RNA count matrices (CSV, H5AD, H5, or sparse matrix)
- Configuration: `config/config.yaml` and `config/config_cell_clusters.tsv`

## Key outputs

- `scE2G_predictions.tsv.gz` - All enhancer-gene predictions with scores
- `scE2G_predictions_threshold{X}.tsv.gz` - Filtered predictions
- `predictions_qc_report.html` - QC report

## Configuration

Main config in `config/config.yaml`:
- `data_type`: "multiome" or "scATAC"
- `model_dir`: path to model (default: `models/multiome_powerlaw_v3/`)
- Cell clusters defined in `config/config_cell_clusters.tsv`

## Testing

See `TESTING.md` for the full picture (test tiers, CircleCI caching/nightly builds, known
coverage gaps). Quick local run:

```bash
pytest -s tests/test_sce2g_apply.py
```

Tests use chr22-scoped example data with expected outputs in `tests/expected_output/`.

## Code conventions

- Snakemake rules in separate `.smk` files under `workflow/rules/`
- R scripts use data.table, dplyr, tidyverse style
- Python scripts use click for CLI arguments
- Log transformations use epsilon=0.01
- NaN/Inf values replaced with 0 during feature processing

## CI workflow (`workflow/Snakefile_ci`) — on `confidence-intervals-cell-type-comparison`, not `dev`/`main`

Separate Snakemake workflow that quantifies score uncertainty via cell subsampling. Run after a full scE2G run.

```bash
snakemake -j1 --use-conda --snakefile workflow/Snakefile_ci --configfile config/config_ci.yaml
```

**Key config fields:** `full_run_results_dir`, `scratch_dir`, `results_dir`, `clusters`, `model`, `n_subsamples`, `mode` ("full" or "abc").

**Mode A ("full"):** reruns entire pipeline per subsample from fragments.
**Mode B ("abc"):** seeds MACS2 peaks from the full run; only reruns from neighborhoods onward (faster).

**Outputs per cluster:** `encode_e2g_predictions_ci_annotated.tsv.gz` (with CI columns: `E2G.Score.CI_low/high`, `E2G.Score.pct_above_threshold`) + `ci_qc_report.html`.

**Rule files:** `workflow/rules/ci_subsample.smk`, `ci_run.smk`, `ci_analysis.smk`.
**Scripts:** `workflow/scripts/ci/` — `determine_subsample_params.py`, `subsample_cells.py`, `compute_ci_annotations.R`, `plot_ci.R`.

## Cross-cell-type comparison (`workflow/scripts/visualization/`) — on `confidence-intervals-cell-type-comparison`, not `dev`/`main`

Standalone scripts (not yet in a Snakemake rule) for comparing predictions across cell types. Uses a **pixi** environment — run with `pixi run --manifest-path workflow/scripts/visualization/pixi.toml Rscript ...`.

**Step 1 — precompute** (run once, slow):
```bash
Rscript workflow/scripts/visualization/precompute_comparisons.R \
  --ci_dir <ci_results_dir> --model <model_name> --threshold 0.177 --output_dir <comparison_dir>
```
Outputs: `consensus_elements.tsv`, `consensus_scores.tsv.gz`, `link_specificity.tsv`, `gene_summary.tsv`, `pairwise_metrics.tsv`.

**Consensus elements:** All peak elements across all cell types are pooled and merged with `GenomicRanges::reduce()` (default gap = 0, so only overlapping peaks merge). Each original element maps to exactly one consensus element. Links are then compared across cell types using `(consensus_id, TargetGene)` as the shared key — slightly more permissive than strict reciprocal overlap but transitive and consistent. Use `--merge_gap 500` (or larger) for a neighborhood-level catalog.

**Step 2 — per-gene bubble plot + locus plot:**
```bash
Rscript workflow/scripts/visualization/plot_gene_comparison.R \
  --gene <GENE> --comparison_dir <comparison_dir> --threshold 0.177 --max_elements 75
```
Produces two PDFs: `{gene}.pdf` (bubble matrix) and `{gene}_locus.pdf` (locus plot).

- **Bubble matrix:** rows = consensus elements (genomic order), columns = cell types, dot size = ATAC signal, dot color = E2G score (purple gradient), dot border = black if above threshold. Includes distance-from-TSS panel on the left. Top strip shows RNA TPM and promoter ATAC per cell type.
- **Locus plot:** continuous genomic x-axis, one facet per cell type. ATAC signal shown as downward bars (height = normalized ATAC, color = element status). Above-threshold elements outlined in black. Arcs drawn above the axis connecting element midpoints to TSS, colored and weighted by E2G score. Optional `--bw_dir` adds bigwig coverage track above the bars.

**Step 3 — rank genes by enhancer landscape variability:**
```bash
Rscript workflow/scripts/visualization/rank_variable_genes.R \
  --comparison_dir <comparison_dir> --output <comparison_dir>/variable_genes.tsv
```
Ranks all genes by how much their enhancer landscape varies across cell types, using two metrics averaged into a combined rank:
- **Score variance (option 2):** for each (gene, element), compute variance of E2G scores across cell types (absent = 0); aggregate per gene as mean variance across elements.
- **n_links CV (option 3):** coefficient of variation of `n_links_above_threshold` across cell types per gene.
- **CI-weighted score variance (option 4):** same as option 2 but scores are multiplied by `E2G.Score.pct_above_threshold` before computing variance, downweighting elements with unstable CI. Skipped gracefully if CI columns are absent from `consensus_scores.tsv.gz`.

Combined rank is the equal-weight average of all available option ranks.

**Other visualization scripts:** `plot_global_comparison.R` (pairwise metrics), `plot_expression_vs_enhancers.R`.

## Test data

WTC11 differentiation timecourse (iPSC → EC, 5 cell types): `results/2026_0210_WTC11_EC/`
- Cell types: `d0_ipsc`, `d1_ps`, `d2_meso`, `d3_ec`, `d4_ec`
- Precomputed comparison outputs: `results/2026_0210_WTC11_EC/comparison/`
- Bigwigs: `results/2026_0210_WTC11_EC/{cell_type}/ATAC_norm.bw`
- Example gene for locus plots: KDR (highly expressed in d4_ec)

## Important notes

- `ENCODE_rE2G` is a git submodule that itself nests `ABC` as a submodule - run `git submodule update --init --recursive` if missing
- Conda environments auto-install on first Snakemake run with `--use-conda`
- Model training requires CRISPR validation data (K562 included in resources)

## Troubleshooting

### SLURM cpus-per-task conflict error

**Error message:**
```
srun: fatal: cpus-per-task set by two different environment variables SLURM_CPUS_PER_TASK=4 != SLURM_TRES_PER_TASK=cpu=1
```

**Cause:** When running Snakemake with a SLURM profile that uses `--export=ALL`, the parent shell's `SLURM_CPUS_PER_TASK` variable (e.g., from an interactive job) is exported to submitted jobs. This conflicts with `--cpus-per-task {threads}` which sets `SLURM_TRES_PER_TASK` on newer SLURM versions.

**Solution:** Wrap the sbatch command in your Snakemake SLURM profile (`~/.config/snakemake/slurm/config.yaml`) to unset the conflicting variables before submission:

```yaml
# Before (broken):
cluster: "sbatch --parsable --export=ALL --ntasks 1 --cpus-per-task {threads} ..."

# After (fixed):
cluster: "bash -c 'unset SLURM_CPUS_PER_TASK SLURM_TRES_PER_TASK; sbatch --parsable --export=ALL --ntasks 1 --cpus-per-task {threads} ...'"
```

**What doesn't work:**
- `unset SLURM_CPUS_PER_TASK && snakemake ...` — doesn't reliably fix it
- Adding `unset` inside a rule's `shell:` block — error occurs at job submission, before the shell runs
