# scE2G workflow review

*Review date: 2026-02-09*

This document provides a comprehensive review of the scE2G Snakemake workflow, particularly following the integration of ENCODE_rE2G rules and the refactoring to use ABC as a submodule. The review evaluates correctness, robustness, efficiency, and simplicity.

---

## Executive summary

The scE2G workflow is well-structured with clear separation of concerns across rule files. The recent flattening of ENCODE_rE2G into the main repository has improved maintainability.

**Completed optimizations:**
- ✓ Replaced checkpoint-based conditional execution with upfront configuration evaluation
- ✓ Consolidated fragment processing rules (counting merged into tagAlign)
- ✓ Combined bigWig generation rules (single `bedtools genomecov` call)
- ✓ Added explicit `.tbi` index outputs to prevent silent failures

**Remaining issues:**
- Some join operations have inconsistent keys that could lead to subtle bugs
- Job granularity could be further improved for some operations

**Remaining recommended actions:**
1. Standardize resource specification pattern (simplicity)
2. Batch feature generation into single rule (efficiency)
3. Add integration tests for feature tables (robustness)

---

## 1. Architecture overview

### 1.1 Module structure

```
Snakefile (main)
├── ABC module (submodule at ABC/)
│   ├── rules/macs2.smk
│   ├── rules/candidate_regions.smk
│   ├── rules/neighborhoods.smk
│   └── rules/predictions.smk
└── scE2G rules (workflow/rules/)
    ├── utils.smk              # Helper functions, config handling
    ├── frag_to_tagAlign.smk   # Fragment processing (4 rules)
    ├── e2g_genomewide_features.smk  # Feature generation (7 rules)
    ├── e2g_predictions.smk    # Prediction rules (4 rules)
    ├── add_external_features.smk    # External feature integration (1 rule)
    ├── make_kendall_pairs.smk       # Kendall prep (1 rule)
    ├── generate_atac_matrix.smk     # Matrix generation (2 rules)
    ├── compute_kendall.smk          # Kendall computation (1 rule)
    ├── arc_e2g.smk                  # ARC-E2G integration (1 rule)
    ├── sc_predictions.smk           # Model application + QC (8 rules)
    └── save_configs.smk             # Config export (1 rule)
```

### 1.2 Data flow

```
Fragment files + RNA matrices
         │
         ▼
┌─────────────────────────────────────┐
│  Fragment processing                │
│  - frag_to_tagAlign                 │
│  - process_fragment_file (optional) │
│  - fragment_count                   │
└─────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────┐
│  ABC pipeline (submodule)           │
│  - MACS2 peak calling               │
│  - Candidate regions                │
│  - Neighborhoods                    │
│  - EnhancerPredictionsAllPutative   │
└─────────────────────────────────────┘
         │
         ├──────────────────────────────┐
         ▼                              ▼
┌──────────────────────┐    ┌──────────────────────┐
│  Basic features      │    │  Kendall/ARC         │
│  (config-gated)      │    │  (config-gated)      │
│  - numCandidateEnhG  │    │  - ATAC matrix       │
│  - numTSSEnhGene     │    │  - Kendall corr      │
│  - numNearbyEnh      │    │  - ARC-E2G score     │
└──────────────────────┘    └──────────────────────┘
         │                              │
         └──────────────┬───────────────┘
                        ▼
┌─────────────────────────────────────┐
│  Feature integration                │
│  - activity_only_features           │
│  - add_external_features            │
│  - gen_final_features               │
└─────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────┐
│  Model application                  │
│  - run_e2g_qnorm                    │
│  - filter_e2g_predictions           │
│  - write_predictions_bedpe          │
└─────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────┐
│  QC and reporting                   │
│  - element_and_gene_summaries       │
│  - get_stats_per_model_per_cluster  │
│  - plot_stats                       │
│  - hover_plots → HTML report        │
└─────────────────────────────────────┘
```

---

## 2. Correctness issues

### 2.1 Inconsistent join keys in activity_only_features.R

**Location:** `workflow/scripts/feature_tables/activity_only_features.R:62-82`

**Issue:** Different feature tables are joined with different keys:
- `NumCandidateEnhGene` and `NumTSSEnhGene`: joined by `(name, TargetGene)`
- `NumEnhancersEG5kb` and `SumEnhancersEG5kb`: joined by `name` only

**Impact:** An enhancer with multiple target genes will share the same `numNearbyEnhancers` and `sumNearbyEnhancers` values across all its E-G pairs. This is likely intentional (these features describe the enhancer's neighborhood, not the E-G relationship), but should be documented.

**Recommendation:** Add comment clarifying this is intentional behavior, or verify this is the desired logic.

### 2.2 Checkpoint feature detection is string-based

**Status:** ✓ RESOLVED

**Previous issue:** Feature requirement detection used string substring matching in checkpoints.

**Resolution:** Checkpoints replaced with upfront configuration evaluation in `workflow/rules/utils.smk`. Now uses exact feature name matching via `compute_biosample_feature_requirements()`.

### 2.3 ARC score detection logic

**Status:** ✓ RESOLVED

**Previous issue:** Dead code path where checkpoint never wrote "Kendall" but `features_to_generate` checked for it.

**Resolution:** Checkpoints removed entirely. Feature requirements computed upfront using `BIOSAMPLE_NEEDS_ARC` dictionary lookup.

### 2.4 Fill value edge case

**Location:** `workflow/scripts/feature_tables/gen_final_features.R` (via `get_fill_values.R`)

**Issue:** When `fill_value = "mean"`, the mean is computed excluding Inf/-Inf values. If all values are Inf/-Inf, this returns NaN, which then fails to replace NAs properly.

**Recommendation:** Add defensive handling for edge cases where all values are non-finite.

---

## 3. Robustness issues

### 3.1 Checkpoint output used as input without validation

**Status:** ✓ PARTIALLY RESOLVED

**Previous issue:** Functions returned `RESULTS_DIR` (a directory) when features not required.

**Resolution:** Checkpoints removed. The pattern of returning `RESULTS_DIR` as a placeholder is still used, but is now determined at DAG construction time rather than runtime. The R scripts still check `file_test("-f", ...)` before reading.

**Remaining consideration:** Consider using Snakemake's `touch()` or proper optional inputs for cleaner semantics, but current approach works reliably.

### 3.2 Missing error handling in external feature merging

**Location:** `workflow/scripts/feature_tables/merge_external_features.R:115-130`

**Issue:** The script doesn't validate that source files exist before attempting to read them:
```r
this_source = fread(unique_sources[i])
```

**Impact:** Cryptic error messages if a source file is missing or malformed.

**Recommendation:** Add file existence check and informative error message.

### 3.3 Memory allocation may exceed cluster limits

**Location:** `workflow/rules/utils.smk:23-30`

**Issue:** Memory calculation assumes 8x compression ratio for gzipped files, which may underestimate for highly compressed data or overestimate for already-decompressed intermediate files.

```python
if ".gz" in str(input):
    input_size_mb *= 8  # assume gz compressed the file <= 8x
```

**Impact:** Jobs may OOM with insufficient memory or request excessive resources.

**Recommendation:** Consider making the compression ratio configurable or using actual decompressed size from file headers where available.

### 3.4 Hardcoded 5kb distance threshold

**Location:** `workflow/rules/e2g_genomewide_features.smk:124-125`

**Issue:** The `generate_num_sum_enhancers` rule has a `{kb}` wildcard, but the inputs to `activity_only_features` are hardcoded to `5kb`:
```python
NumEnhancersEG5kb = get_numNearbyEnhancers_file,  # Always 5kb
```

**Impact:** If multiple distance thresholds were intended, this would silently use only 5kb.

**Recommendation:** Either make this configurable or remove the wildcard if only 5kb is supported.

---

## 4. Efficiency issues

### 4.1 Job granularity analysis

The following table summarizes rules where job granularity could be improved:

| Rule group | Current jobs (N clusters, M models) | Recommended | Reduction |
|------------|-------------------------------------|-------------|-----------|
| Fragment counting | N | Merge into tagAlign | N → 0 |
| BigWig generation | 2N | Combine norm/unnorm | 50% |
| Kendall/ARC chain | 5N | Single rule | 80% |
| Feature generation | 3N (conditional) | Batch per cluster | 67% |
| QC aggregation | 2 | Single rule | 50% |
| Stats per cluster | N×M | Batch stats | 80-90% |

**Estimated total reduction:** For 10 clusters × 2 models: ~150 jobs → ~50 jobs

### 4.2 Fragment processing consolidation

**Status:** ✓ RESOLVED

**Previous state:** Three separate operations for fragments with conditional rule definitions.

**Resolution:** Fragment processing consolidated from 5 rules to 4 rules:
1. `process_fragment_file`: Filter chromosomes (only runs when `fragments_preprocessed=False`)
2. `ensure_fragment_index`: Create index on-the-fly (only runs when `fragments_preprocessed=True`)
3. `frag_to_tagAlign`: Convert to tagAlign + count fragments (always runs)
4. `frag_to_bigwig`: Generate both normalized and unnormalized bigWig files (when `make_IGV_tracks=True`)

All rules now explicitly output `.tbi` index files to prevent silent failures.

### 4.3 BigWig generation duplication

**Status:** ✓ RESOLVED

**Previous state:** Two rules (`frag_to_bigWig`, `frag_to_norm_bigWig`) each ran `bedtools genomecov` separately.

**Resolution:** Consolidated into single `frag_to_bigwig` rule that:
1. Generates sorted bedgraph once
2. Creates unnormalized bigWig from bedgraph
3. Scales bedgraph values and creates normalized bigWig
4. Cleans up temporary bedgraph file

### 4.4 Checkpoint serialization barrier

**Status:** ✓ RESOLVED

**Previous state:** Two checkpoints (`basic_features_required`, `features_required`) created sequential barriers.

**Resolution:** Checkpoints replaced with upfront configuration evaluation:
- `compute_biosample_feature_requirements()` in `workflow/rules/utils.smk` reads model feature tables at initialization
- Creates dictionaries: `BIOSAMPLE_NEEDS_ARC`, `BIOSAMPLE_NEEDS_NUMCANDIDATEENHGENE`, etc.
- Getter functions use simple dict lookups instead of checkpoint outputs
- Full DAG is now constructed at startup, enabling better parallelization planning

### 4.5 QC stats aggregation chain

**Current state:**
1. `get_stats_per_model_per_cluster`: Per (cluster, model) stats
2. `plot_stats`: Aggregate all stats, generate plots
3. `hover_plots`: Generate HTML report from aggregated stats

**Issue:** Two sequential aggregation steps.

**Recommendation:** Combine `plot_stats` and `hover_plots` into a single rule:
```python
rule generate_qc_report:
    input:
        stats_files = [...]
    output:
        all_stats = "all_qc_stats.tsv",
        html_report = "predictions_qc_report.html",
        pdf_plots = [...]
    script:
        "scripts/prediction_qc/generate_full_qc_report.R"
```

---

## 5. Simplicity issues

### 5.1 Redundant biosample expansion

**Location:** `workflow/Snakefile:95-96` and `workflow/rules/sc_predictions.smk:124`

**Issue:** `biosample_model_threshold` tuple is created twice with identical logic.

**Recommendation:** Define once in Snakefile and use consistently.

### 5.2 Inconsistent resource specification

**Location:** Various rule files

**Issue:** Resources are specified differently across rules:
- `mem_mb=32*1000` (hardcoded)
- `mem_mb=ABC.determine_mem_mb` (dynamic)
- `mem_mb=partial(determine_mem_mb, min_gb=16)` (partial)
- `mem_mb=determine_mem_mb` (local)

**Recommendation:** Standardize on one approach. Prefer `partial(determine_mem_mb, min_gb=X)` for consistency.

### 5.3 Mixed use of lambda and direct function calls

**Location:** Throughout rule files

**Issue:** Some rules use lambda for params:
```python
params:
    threshold = lambda wildcards: get_model_threshold(...)
```

While others use direct function calls in the Snakefile body.

**Recommendation:** Document the convention. Generally:
- Lambda: When wildcards are needed
- Direct call: When config values suffice

### 5.4 Overly complex path construction

**Location:** `workflow/Snakefile:93-108`

**Issue:** Output file lists are constructed with list comprehensions that duplicate path patterns:
```python
encode_e2g_predictions = [os.path.join(RESULTS_DIR, biosample, model_name, "encode_e2g_predictions.tsv.gz") for biosample,model_name in biosample_model]
prediction_stats = [os.path.join(RESULTS_DIR, biosample, model_name, f"encode_e2g_predictions_threshold{threshold}_stats.tsv") for biosample,model_name,threshold in biosample_model_threshold]
```

**Recommendation:** Use Snakemake's `expand()` for consistency:
```python
encode_e2g_predictions = expand(
    os.path.join(RESULTS_DIR, "{biosample}", "{model_name}", "encode_e2g_predictions.tsv.gz"),
    zip, biosample=BIOSAMPLE_DF["biosample"], model_name=BIOSAMPLE_DF["model_dir_base"]
)
```

---

## 6. Proposed refactoring

### 6.1 Immediate changes (low effort, high impact)

#### 6.1.1 Consolidate fragment counting

**Files:** `workflow/rules/frag_to_tagAlign.smk`

Merge `get_fragment_count` into `frag_to_tagAlign` by adding fragment_count as an additional output:
```python
rule frag_to_tagAlign:
    input:
        frag_file = get_processed_fragment_file
    output:
        tagAlign = "...",
        fragment_count = "..."  # New output
    shell:
        """
        # Existing tagAlign logic...

        # Add fragment counting
        zcat {input.frag_file} | wc -l > {output.fragment_count}
        """
```

**Impact:** Eliminates N extra jobs per run.

#### 6.1.2 Combine QC report generation

**Files:** `workflow/rules/sc_predictions.smk`

Merge `plot_stats` and `hover_plots` into single rule.

**Impact:** Eliminates 1 job and reduces complexity.

#### 6.1.3 Fix dead code path

**Files:** `workflow/rules/add_external_features.smk`

Remove the `val == "Kendall"` branch or fix the checkpoint logic.

### 6.2 Medium-term changes (medium effort, high impact)

#### 6.2.1 Consolidate bigWig generation

**Files:** `workflow/rules/frag_to_tagAlign.smk`

Create single rule that generates both normalized and unnormalized bigWig files from one bedgraph computation.

**Impact:** Reduces I/O by 50% for bigWig generation.

#### 6.2.2 Batch feature generation

**Files:** `workflow/rules/e2g_genomewide_features.smk`

Create a combined feature generation script that:
1. Reads feature table once
2. Generates all required features in a single pass
3. Outputs individual feature files

```python
rule generate_all_derived_features:
    input:
        feature_table = "...",
        abc_predictions = "...",
        gene_tss = "..."
    output:
        numCandidateEnhGene = "...",
        numTSSEnhGene = "...",
        numNearbyEnhancers = "...",
        sumNearbyEnhancers = "..."
    script:
        "scripts/feature_tables/generate_all_derived_features.py"
```

**Impact:** Reduces 3-4 jobs per cluster to 1.

### 6.3 Long-term changes (high effort, high impact)

#### 6.3.1 Replace checkpoints with upfront configuration

**Files:** `workflow/Snakefile`, `workflow/rules/add_external_features.smk`, `workflow/rules/e2g_genomewide_features.smk`

Move feature requirement detection from checkpoints to Snakefile initialization:

```python
# In Snakefile, after loading model config
REQUIRED_FEATURES = {}
for biosample in BIOSAMPLE_DF["biosample"]:
    model_dir = get_model_dir(biosample)
    feature_table = pd.read_csv(os.path.join(model_dir, "feature_table.tsv"), sep="\t")
    REQUIRED_FEATURES[biosample] = {
        "numCandidateEnhGene": "numCandidateEnhGene" in feature_table["input_col"].values,
        "kendall": "Kendall" in feature_table["input_col"].values,
        # etc.
    }
```

Then use conditional rule inclusion or input functions that reference this dict.

**Impact:** Removes serialization barriers, enables full DAG optimization.

#### 6.3.2 Implement scatter-gather for Kendall computation

**Files:** `workflow/rules/compute_kendall.smk`

If Kendall computation is a bottleneck, split by chromosome:

```python
rule compute_kendall_per_chr:
    input:
        pairs = "Pairs_{chr}.tsv.gz",
        # ...
    output:
        kendall = "Pairs.Kendall_{chr}.tsv.gz"

rule merge_kendall:
    input:
        expand("Pairs.Kendall_{chr}.tsv.gz", chr=CHROMOSOMES)
    output:
        "Pairs.Kendall.tsv.gz"
```

**Impact:** Enables parallel Kendall computation across chromosomes.

---

## 7. Testing recommendations

### 7.1 Add integration tests for feature table construction

The feature table construction is complex with multiple join operations. Add tests that verify:
- All expected features are present in final output
- No unexpected NaN values in non-optional features
- Row counts match expectations

### 7.2 Add validation for checkpoint outputs

Add validation rules that check checkpoint output files contain expected values:
```python
rule validate_checkpoint_outputs:
    input:
        "to_generate.txt"
    run:
        with open(input[0]) as f:
            val = f.read().strip()
            assert val in ["Kendall", "ARC", "Neither"], f"Invalid checkpoint value: {val}"
```

### 7.3 Add dry-run CI check

Add a CI job that runs `snakemake -n` to verify the DAG can be constructed without errors.

---

## 8. Summary of recommendations

### Priority 1 (implement first)
1. **DONE:** Fix dead code path in `features_to_generate` (correctness) - Replaced checkpoints with upfront config evaluation
2. **DONE:** Consolidate fragment counting into tagAlign rule (efficiency) - Fragment count now produced as output of `frag_to_tagAlign`
3. Combine QC report rules (simplicity)

### Priority 2 (implement soon)
4. **DONE:** Consolidate bigWig generation (efficiency) - Single `frag_to_bigwig` rule now produces both normalized and unnormalized outputs
5. Standardize resource specification pattern (simplicity)
6. Add comments for intentional join key differences (correctness)

### Priority 3 (implement when feasible)
7. Batch feature generation into single rule (efficiency)
8. **DONE:** Replace checkpoints with upfront config evaluation (efficiency) - Feature requirements computed at Snakefile initialization
9. Add integration tests for feature tables (robustness)

---

## Appendix A: Rule inventory

| File | Rule | Type | Inputs | Outputs | Memory |
|------|------|------|--------|---------|--------|
| frag_to_tagAlign.smk | frag_to_tagAlign | shell | frag_file, frag_file_index | tagAlign.sort.gz, tagAlign.sort.gz.tbi, fragment_count.txt | dynamic |
| frag_to_tagAlign.smk | process_fragment_file | shell | frag_file | fragments_filtered.tsv.gz, fragments_filtered.tsv.gz.tbi | dynamic |
| frag_to_tagAlign.smk | ensure_fragment_index | shell | frag_file | input_fragments.tsv.gz.tbi | dynamic |
| frag_to_tagAlign.smk | frag_to_bigwig | shell | frag_file, fragment_count | ATAC.bw, ATAC_norm.bw | dynamic |
| e2g_predictions.smk | make_biosample_feature_table | script | ABC_BIOSAMPLES | feature_table.tsv | 4GB |
| e2g_predictions.smk | filter_e2g_predictions | shell | predictions | thresholded.tsv.gz | dynamic |
| e2g_predictions.smk | write_predictions_bedpe | shell | thresholded | .bedpe | dynamic |
| e2g_genomewide_features.smk | generate_num_candidate_enh_gene | shell | abc_predictions | NumCandidateEnhGene.tsv | 16GB min |
| e2g_genomewide_features.smk | generate_num_tss_enh_gene | shell | abc_predictions, gene_tss | NumTSSEnhGene.tsv | 32GB min |
| e2g_genomewide_features.smk | generate_num_sum_enhancers | shell | abc_predictions, enhancer_list | NumEnhancersEG.txt, SumEnhancersEG.txt | 16GB min |
| e2g_genomewide_features.smk | activity_only_features | script | feature_table, abc, derived features | ActivityOnly_features.tsv.gz | dynamic |
| e2g_genomewide_features.smk | add_external_features | script | predictions_extended, feature_table, external_config | plus_external_features.tsv.gz | 8-32GB |
| e2g_genomewide_features.smk | gen_final_features | script | plus_external_features, feature_table | genomewide_features.tsv.gz | dynamic |
| add_external_features.smk | make_external_features_config | script | feature_inputs | external_features_config.tsv | 8GB |
| compute_kendall.smk | compute_kendall | script | pairs, atac_matrix, rna_matrix | Pairs.Kendall.tsv.gz | 63GB min |
| arc_e2g.smk | arc_e2g | script | abc_predictions, kendall_predictions | ARC predictions | 64GB |
| sc_predictions.smk | overlap_features_crispr_apply | script | predictions, crispr | benchmark features | 32GB |
| sc_predictions.smk | crispr_benchmarking | shell | crispr_features | performance summary | 64GB |
| sc_predictions.smk | run_e2g_qnorm | shell | final_features | predictions.tsv.gz | dynamic |
| sc_predictions.smk | element_and_gene_summaries | script | predictions | gene/element lists | dynamic |
| sc_predictions.smk | get_stats_per_model_per_cluster | script | predictions | stats.tsv | 32GB |
| sc_predictions.smk | plot_stats | script | stats_files | QC plots | dynamic |
| sc_predictions.smk | hover_plots | script | qc_stats | HTML report | dynamic |
| save_configs.smk | save_reference_configs | run | - | config files | - |

---

## Appendix B: Configuration dependencies

```
config.yaml
├── cell_clusters → CELL_CLUSTER_DF
├── model_dir → MODEL_DIR (default model directory)
├── abc_dir → ABC module location
├── results_dir → RESULTS_DIR
├── gene_TSS500 → Reference TSS annotations
├── chr_sizes → Chromosome sizes
├── gene_annotations → GTF file
├── gene_classes → Ubiquitously expressed genes
├── crispr_dataset → Validation data
├── qc_reference → Reference clusters for QC
├── data_type → "multiome" or "scATAC"
├── final_score_col → Score column for filtering
├── include_self_promoter → Boolean
├── make_IGV_tracks → Boolean
├── benchmark_performance → Boolean
├── fragments_preprocessed → Boolean
└── RNA_matrix_filtered → Boolean

config_cell_clusters.tsv
├── cluster → Cell cluster name
├── atac_frag_file → Path to fragment file
├── rna_matrix_file → Path to RNA matrix
├── HiC_file → Optional HiC data
└── HiC_type → "avg", "hic", or null
```
