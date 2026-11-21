# Plan: Split Kendall correlation by chromosome

## Overview

Split the `compute_kendall` rule by chromosome when cell count exceeds a configurable threshold (default 10,000). This reduces wall time from 12-36 hours to ~1-2 hours and per-job memory from 63-250+ GB to ~16-30 GB for large datasets. Small datasets (below threshold) use the existing single-job behavior unchanged.

## Architecture

A new Snakemake **checkpoint** (`kendall_precompute`) loads the matrices, computes RNA metrics, saves pre-processed matrices to disk, and writes a routing flag ("single" or "split"). A routing input function reads the flag and directs downstream rules to either the existing single-job `compute_kendall` output or merged per-chromosome outputs.

```
                    kendall_precompute (checkpoint)
                    ├── umi_count.txt, cell_count.txt, gene_expression_metrics.tsv.gz
                    ├── rna_matrix_filtered.rds, atac_matrix_binary.rds
                    └── kendall_routing.txt ("single" or "split")
                              │
                ┌─────────────┴─────────────┐
                │ routing="single"          │ routing="split"
                │                           │
        compute_kendall (existing)    compute_kendall_per_chr × N chromosomes
        → single/Pairs.Kendall.tsv.gz       → split/Pairs.Kendall.{chr}.tsv.gz
                │                           │
                │                     merge_kendall_chunks
                │                     → split/Pairs.Kendall.merged.tsv.gz
                │                           │
                └─────────┬─────────────────┘
                  route_kendall_output
                  → Pairs.Kendall.tsv.gz (canonical path, consumed by arc_e2g etc.)
```

Metadata files (`umi_count.txt`, `cell_count.txt`, `gene_expression_metrics.tsv.gz`) are always produced by the precompute step and copied to canonical paths by `route_kendall_metadata`. No downstream rules need changes.

## Files to modify

### 1. `config/config.yaml`
Add: `kendall_split_threshold: 10000`

### 2. `tests/config/test_config.yml`
Add: `kendall_split_threshold: 10000`

### 3. `workflow/Snakefile` (~line 32)
Add default:
```python
if "kendall_split_threshold" not in config:
    config["kendall_split_threshold"] = 10000
```

### 4. `workflow/rules/compute_kendall.smk` (full rewrite)
Replace the single `compute_kendall` rule with 5 rules + 2 input functions:

**a. `checkpoint kendall_precompute`**
- Input: `Pairs.tsv.gz`, `atac_matrix.rds`, `rna_matrix`
- Output (in `precompute/` subdir): `umi_count.txt`, `cell_count.txt`, `gene_expression_metrics.tsv.gz`, `chromosomes.txt`, `kendall_routing.txt`, `rna_matrix_filtered.rds`, `atac_matrix_binary.rds`
- Params: `gene_gtf`, `abc_genes`, `cell_threshold`
- Resources: `min_gb=63`, runtime 2h
- Script: `kendall_precompute.R` (new)

**b. `rule compute_kendall`** (existing script, new output paths)
- Same inputs/params as current rule
- Output paths moved to `single/` subdir: `single/Pairs.Kendall.tsv.gz`, `single/umi_count.txt`, `single/cell_count.txt`, `single/gene_expression_metrics.tsv.gz`
- Script: existing `compute_kendall.R` **unchanged** (reads paths from `snakemake@output`)
- Resources: same as current (63 GB min, 12h)

**c. `rule compute_kendall_per_chr`**
- Input: `Pairs.tsv.gz`, `precompute/rna_matrix_filtered.rds`, `precompute/atac_matrix_binary.rds`
- Params: `chromosome = "{chr}"`
- Output: `split/Pairs.Kendall.{chr}.tsv.gz`
- Script: `compute_kendall_per_chr.R` (new)
- Resources: `min_gb=16`, runtime 6h

**d. `rule merge_kendall_chunks`**
- Input: all per-chr Kendall files (via aggregate function reading `chromosomes.txt`), `precompute/gene_expression_metrics.tsv.gz`
- Output: `split/Pairs.Kendall.merged.tsv.gz`
- Script: `merge_kendall_chunks.py` (new)
- Resources: 32 GB

**e. `rule route_kendall_output`**
- Input: `get_kendall_output()` function reads checkpoint routing flag, returns either `single/` or `split/` path
- Output: canonical `Pairs.Kendall.tsv.gz`
- Shell: `cp {input} {output}`

**f. `rule route_kendall_metadata`**
- Input: `precompute/umi_count.txt`, `precompute/cell_count.txt`, `precompute/gene_expression_metrics.tsv.gz`
- Output: canonical `umi_count.txt`, `cell_count.txt`, `gene_expression_metrics.tsv.gz`
- Shell: `cp` each file

## Files to create

### 5. `workflow/scripts/feature_computation/kendall_precompute.R`
Extracted from `compute_kendall.R` lines 1-231 (everything except the Kendall computation). Additional logic:
- Extract unique chromosomes from `Pairs.tsv.gz` → write `chromosomes.txt`
- Compare `ncol(rna_matrix)` to threshold → write `kendall_routing.txt` ("single" or "split")
- Save binarized ATAC matrix and filtered/normalized RNA matrix as RDS
- If routing="single", save placeholder (empty) matrix files to avoid writing multi-GB files that won't be used
- Includes `extract_attributes()` and `map_gene_names()` functions (copied from `compute_kendall.R`)

### 6. `workflow/scripts/feature_computation/compute_kendall_per_chr.R`
Streamlined Kendall computation for a single chromosome:
- Loads pre-saved RDS matrices from precompute step (no re-normalization needed)
- Filters `Pairs.tsv.gz` to specified chromosome
- Handles empty chromosome case (writes header-only file)
- Includes `count_diff` Rcpp function, `kendall_one_gene`, `kendall_multiple_genes` (copied from `compute_kendall.R`)
- Output columns: `chr, start, end, TargetGene, PeakName, PairName, Kendall` (no expression metrics)

### 7. `workflow/scripts/feature_computation/merge_kendall_chunks.py`
Python script that:
- Concatenates per-chromosome Kendall TSV files
- Reads `gene_expression_metrics.tsv.gz`, joins expression columns by TargetGene
- Maps column names: `RNA_meanLogNorm` → `mean_log_normalized_rna`, `RNA_pseudobulkTPM` → `RnaPseudobulkTPM`, `RNA_percentCellsDetected` → `RnaDetectedPercent`
- Writes final output matching the exact column order of the original `Pairs.Kendall.tsv.gz`

## Key design decisions

1. **Chromosomes derived from `Pairs.tsv.gz`**, not hardcoded. Only chromosomes with actual pairs get jobs (no wasted jobs for chrM/chrY if no pairs exist).

2. **Pre-saved matrices** in the precompute step avoid 22× redundant loading/normalization in per-chr jobs. The precompute step saves the gene-mapped, normalized RNA matrix and binarized ATAC matrix as RDS. For single-path routing, empty placeholder files are saved instead.

3. **`compute_kendall.R` is unchanged.** The single-job path uses the exact same script — only the Snakemake output paths change (moved to `single/` subdir). Since the R script reads paths from `snakemake@output`, this is transparent.

4. **Downstream rules unaffected.** All canonical output paths remain the same. `arc_e2g`, `get_gex_file()`, `get_count_file()`, and `features_to_generate()` continue to reference the same paths.

5. **Single-path overhead is small.** The precompute step adds ~5-10 minutes of matrix loading for small datasets. The `compute_kendall` rule then re-loads matrices (since it uses the unmodified script). This duplication is the tradeoff for keeping `compute_kendall.R` unchanged.

## Testing and verification

1. **Run existing test**: `pytest -s tests/test_sce2g_apply.py` — test data has ~500 cells, below 10k threshold, takes single path. Expected outputs should match exactly.

2. **Test split path**: Create a test config with `kendall_split_threshold: 100` (or temporarily modify the test config), run the pipeline, and verify the merged output matches the single-path output (identical up to row ordering).

3. **Verify canonical paths**: Confirm `umi_count.txt`, `cell_count.txt`, `gene_expression_metrics.tsv.gz`, `Pairs.Kendall.tsv.gz` all appear at the expected canonical locations.

## Implementation order

1. Add config parameter to `config/config.yaml`, `tests/config/test_config.yml`, and `workflow/Snakefile`
2. Create `kendall_precompute.R`
3. Create `compute_kendall_per_chr.R`
4. Create `merge_kendall_chunks.py`
5. Rewrite `compute_kendall.smk`
6. Run existing test (single path)
7. Test split path with low threshold
