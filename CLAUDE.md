# scE2G permutation testing fork

**Project location**: `/oak/stanford/groups/engreitz/Users/sheth/scE2G_permute_cell_types/scE2G`

This is a modified fork of scE2G for testing whether enhancer-gene predictions are cell type-specific by comparing real vs permuted cell samples across a hierarchical cell type structure.

For general scE2G documentation (setup, inputs, outputs, model training), see [README.md](README.md).

---

## Overview

**Goal**: Test if E2G predictions from specific cell types (L3) differ from predictions made using cells randomly sampled from broader groupings (L2 or L1).

**Approach**:
1. Define cell type hierarchy: L3 (finest) → L2 (parent) → L1 (merged all)
2. For N iterations, create three sample types per L3 cell type:
   - `real.L3` - true cells from true cell type
   - `permute.L3_from_L2` - cells sampled from parent L2
   - `permute.L3_from_L1` - cells sampled from merged L1
3. Run E2G predictions on all samples
4. Compare unique E-G pairs across sample types

---

## Modified/new files

### Config
- `config/config.yaml` - Added sampling parameters
- `config/config_cell_types.tsv` - Cell type hierarchy with file paths

### Workflow
- `workflow/Snakefile` - Iteration logic, new rule includes
- `workflow/rules/utils.smk` - Added `get_cell_cluster_config()`
- `workflow/rules/add_external_features.smk` - Global feature table path
- `workflow/rules/sc_predictions.smk` - Changed QC output paths

### New rules
- `workflow/rules/generate_sample_inputs.smk` - Cell sampling
- `workflow/rules/compare_permutations.smk` - Intersection analysis

---

## Configuration

### config.yaml additions

```yaml
cell_type_config: "config/config_cell_types.tsv"
n_iterations: 10
seed_multiplier: 13
n_sample: 30000          # Total cells per iteration
n_L1_sample: 10000       # Subset for L1 analysis
model_dir_use: models/multiome_powerlaw_v3
scratch_dir: /scratch/users/shethm/scE2G_permute_and_compare
results_dir: results/2025_0220_iter10_seed13
```

### config_cell_types.tsv format

```
L3              L2       L3_fragment_file           L2_fragment_file         L1_fragment_file          L3_rna_matrix              L2_rna_matrix            L1_rna_matrix
B1_B            B        /path/to/B1_B.frag.gz      /path/to/B.frag.gz       /path/to/all.frag.gz      /path/to/B1_B.rna.csv.gz   /path/to/B.rna.csv.gz    /path/to/all.rna.csv.gz
Naive_CD20_pos  B        /path/to/Naive.frag.gz     ...                      ...                       ...
CD4_pos_T       T        /path/to/CD4.frag.gz       /path/to/T.frag.gz       ...                       ...
```

- L3 = finest cell type (B1_B, CD4_pos_T_naive, etc.)
- L2 = parent grouping (B, T, Myeloid, etc.)
- L1 = implicit merged level (all cells)
- Only first row per L2 needs L2 files; only one row needs L1 files

---

## Sampling workflow

### Sample naming convention

`iter{N}.{type}.{level}.{identity}`

Examples:
- `iter1.real.L3.B1_B` - Iteration 1, real B1_B cells
- `iter1.permute.L3_from_L2.B1_B` - Cells sampled from B (parent L2), sized to match B1_B
- `iter1.permute.L3_from_L1.B1_B` - Cells sampled from all (L1), sized to match B1_B

### Rule: save_reference_barcodes

Extract unique cell barcodes from each L3 fragment file:
```bash
zcat {L3_fragments} | cut -f 4 | sort | uniq | gzip > {l3_id}_barcodes.txt.gz
```

### Rule: calculate_cell_type_counts

1. Count total cells per L3 from reference barcodes
2. Calculate proportional sampling fraction: `n_sample / total_cells`
3. Output counts per L3, L2, and L1 as pickle files

Proportional sampling ensures cell type composition is preserved (if A has 10% of cells, it contributes 10% of n_sample).

### Rule: get_barcodes_for_iteration

For iteration i with `seed = iter * seed_multiplier`:

```python
# 1. REAL L3: sample from true cell type
sample_l3 = {l3_id: random.sample(all_barcodes[l3_id], l3_counts[l3_id])
             for l3_id in l3_counts}

# 2. Pool L3 samples into L2 groups
sample_l2 = {l2_id: chain(sample_l3[l3_id] for l3_id in l2_dict[l2_id])
             for l2_id in l2_dict}

# 3. PERMUTE L3_from_L2: sample from parent L2 pool (breaks L3 identity)
for l2_id, l2_barcodes in sample_l2.items():
    l2_pool = set(l2_barcodes)
    for l3_id in l2_dict[l2_id]:
        permuted_l3_from_l2[l3_id] = random.sample(l2_pool, l3_counts[l3_id])
        l2_pool -= set(permuted_l3_from_l2[l3_id])  # No duplicates

# 4. PERMUTE L3_from_L1: sample from entire merged pool (breaks L2+L3 identity)
all_pool = set(chain(*sample_l3.values()))
for l3_id in l3_counts:
    permuted_l3_from_l1[l3_id] = random.sample(all_pool, l3_counts[l3_id])
    all_pool -= set(permuted_l3_from_l1[l3_id])
```

### Rule: create_sample_fragment_file

Filter parent fragment file to sample barcodes:
```bash
zgrep -F -f <(zcat sample_barcodes.txt.gz) {input.fragments} | \
    sort -k1,1 -k2,2n --parallel 8 | bgzip > atac_fragments.tsv.gz
tabix -p bed atac_fragments.tsv.gz
```

Key logic in `get_input_file()`:
- `real.L3.X` → use `L3_fragment_file` for X
- `permute.L3_from_L2.X` → use `L2_fragment_file` for parent L2 of X
- `permute.L3_from_L1.X` → use `L1_fragment_file`

### Rule: create_sample_rna_matrix

Filter RNA matrix columns to sample barcodes (preserves gene row order).

---

## Comparison analysis

### Rule: make_pred_intersectable

Extract unique E-G pairs from thresholded predictions:
```bash
zcat predictions.tsv.gz | csvtk cut -t -f TargetGene,start,end | \
    sed 1d | sort | uniq | gzip > gene_start_end_thresholded.tsv.gz
```

### Rule: compute_L1_totals

For each iteration, count unique E-G pairs by taking union across all cell types:
```bash
# Iteratively add non-overlapping pairs using bedtools
for pred in predictions; do
    bedtools intersect -a $pred -b union.txt -f 0.5 -F 0.5 -wa -v >> union.txt
done
```

Outputs per iteration: `L3_real_total`, `L3_from_L2_total`, `L3_from_L1_total`

### Rule: combine_L1_totals

Aggregate across iterations into `L1_totals.tsv`.

---

## Efficiency issues and recommendations

### Summary of issues

| Rule | Severity | Issue |
|------|----------|-------|
| `get_barcodes_for_iteration` | High | Reloads all barcodes from disk for every iteration |
| `get_barcodes_for_iteration` | Medium | `sorted()` called repeatedly inside loops |
| `create_sample_fragment_file` | High | `zgrep` scans entire fragment file for each sample; L1 file read 15× per iteration |
| `create_sample_rna_matrix` | Medium | Matrix decompressed twice; full matrix read for column subsetting |

### Detailed analysis

#### `get_barcodes_for_iteration` (lines 71-136)

**Memory - redundant loading:**
```python
# Line 84: Loads ALL barcodes for ALL cell types into memory
all_barcodes = {l3_id: gzip.open(barcode_file, "rt").read().splitlines()
                for l3_id, barcode_file in zip(ct_config["L3"], input.ref_barcodes)}
```
Runs for every iteration (10×), re-reading the same files. With 15 cell types × 30k cells = ~450k strings loaded repeatedly.

**Compute - unnecessary sorting:**
```python
# Lines 105, 115: sorted() called inside loop
permuted_l3_from_l1[l3_id] = random.sample(sorted(sample_merge), l3_counts[l3_id])
```
`sorted()` is O(n log n) for each L3 cell type. With 15 cell types, sorting 30k items 15 times per permutation type.

**Dead code:** Lines 106-107, 116-117 create unused `this_set` variable.

#### `create_sample_fragment_file` (lines 154-173)

**Major I/O bottleneck:**
```bash
zgrep -F -f <(zcat {input.sample_barcodes}) {input.fragments}
```
O(fragments × barcodes) - for every line in fragment file, grep checks all barcodes. Fragment files have millions of lines.

**Redundant reads of L1 file:**
- 15 L3 cell types × `permute.L3_from_L1` all read from same L1 fragment file
- Per iteration: L1 file scanned 15 times with different barcode subsets
- Total: 10 iterations × 15 = 150 scans of the same large file

#### `create_sample_rna_matrix` (lines 176-190)

**Double decompression:**
```bash
cols=1,$(grep -Fx -n -f <(zcat barcodes) <(zcat {input.matrix} | head -1 | tr ',' '\n') | ...)
zcat {input.matrix} | cut -d',' -f"$cols" | gzip > {output}
```
Matrix decompressed twice: once for header column matching, once to cut columns.

### Recommendations

#### Quick wins (minimal code change)

1. **Sort once, not in loop** in `get_barcodes_for_iteration`:
```python
# Before the loop
sample_merge_sorted = sorted(sample_merge)
for l3_id in l3_counts.keys():
    permuted_l3_from_l1[l3_id] = random.sample(sample_merge_sorted, l3_counts[l3_id])
    # Update the sorted list by removing sampled items
```

2. **Use awk instead of zgrep** for fragment filtering:
```bash
zcat {input.sample_barcodes} | awk 'NR==FNR {bc[$1]; next} $4 in bc' - <(zcat {input.fragments}) | ...
```
O(fragments) with O(1) lookup per line instead of O(fragments × barcodes).

#### Larger refactors (for confidence interval pipeline)

3. **Cache reference barcodes**: Load once at pipeline start, save as single pickle, reuse across iterations.

4. **Batch filter L1 fragments**: Filter L1 file once with all sampled barcodes, then split output by sample assignment:
```python
# Single pass: output to multiple files based on barcode→sample mapping
```

5. **Use Python for RNA matrix subsetting**: pandas/polars can select columns efficiently:
```python
# Read only needed columns with usecols parameter
```

### Estimated impact

| Change | Time savings |
|--------|--------------|
| Sort once instead of in loop | ~15× faster for permutation step |
| awk instead of zgrep | ~10-100× faster for fragment filtering |
| Batch L1 filtering | ~15× fewer reads of L1 file per iteration |
| Cache barcode loading | ~10× fewer file reads across iterations |

---

## Adapting for confidence intervals

For bootstrap-based confidence intervals on E2G predictions:

| Permutation (this fork) | Confidence intervals |
|------------------------|----------------------|
| Sample from parent level | Sample **with replacement** from same cell type |
| Compare real vs permuted | Aggregate scores across bootstrap samples |
| Count unique pairs | Compute mean + CI per E-G pair |

### Reusable components

- `save_reference_barcodes` - Extract barcodes per cell type
- `calculate_cell_type_counts` - Proportional sampling counts
- `create_sample_fragment_file` - Filter fragments to barcodes
- `create_sample_rna_matrix` - Filter RNA matrix to barcodes
- Scratch/results directory separation
- Sample naming: `iter{N}.bootstrap.{cell_type}`

### Key change for bootstrap

```python
# With replacement (bootstrap) instead of without replacement (permutation)
bootstrap_sample = random.choices(all_barcodes[l3_id], k=n_cells)
```

### Aggregation for CI

```python
scores = [get_score(eg_pair, iter) for iter in range(n_bootstrap)]
mean_score = np.mean(scores)
ci_low, ci_high = np.percentile(scores, [2.5, 97.5])
```

---

## Running

```bash
# Activate environment
conda activate run_snakemake

# Dry run
snakemake -n

# Run with SLURM
snakemake --use-conda --profile slurm -j 100
```
