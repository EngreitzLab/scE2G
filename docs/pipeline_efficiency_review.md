# Pipeline efficiency and structure review

## 1. Submodule architecture

**Current structure:** scE2G → ENCODE_rE2G (submodule) → ABC (nested submodule)

Each layer uses Snakemake's `module` directive to import the layer below, prefixing rules (`abc_*`) and passing config down through transformation functions.

### What works

- **Logical separation of concerns.** ABC handles contact-based enhancer prediction, ENCODE_rE2G adds feature engineering and logistic regression, scE2G adds single-cell specific steps (Kendall, ARC-E2G, fragment processing). This maps cleanly to the scientific method layers.
- **ENCODE_rE2G and ABC are independently useful.** They serve the bulk/ENCODE use case without scE2G. Keeping them as separate repos is justified.
- **Config persistence** (`save_configs.smk`) dumps all three resolved configs to the results directory, which is good for reproducibility.

### What doesn't work

**Config plumbing is the main pain point.** Three transformation functions chain configs across layers:

| Function | Location | Maps |
|---|---|---|
| `make_biosample_config()` | `workflow/rules/utils.smk:20` | cell clusters → ABC biosamples |
| `get_e2g_config()` | `workflow/rules/utils.smk:52` | scE2G config → ENCODE_rE2G config |
| `get_abc_config()` | `ENCODE_rE2G/workflow/rules/utils.smk` | ENCODE_rE2G config → ABC config |

Problems with this chain:

1. **Naming inconsistencies.** scE2G calls it `encode_re2g_dir`, ENCODE_rE2G calls it `E2G_DIR_PATH`. scE2G passes `ABC_BIOSAMPLES`, which becomes `biosamplesTable` at the ABC layer. `gene_TSS500` maps to `ref.genome_tss` in ABC's nested config structure.

2. **Flat vs nested config.** ABC nests reference files under a `ref:` key; the other two layers use flat keys. `get_abc_config()` manually maps between them — every new reference file requires changes in multiple places.

3. **Duplicated utility functions.** `make_paths_absolute()` is copy-pasted identically in all three `utils.smk` files. `determine_mem_mb()` exists in both ENCODE_rE2G and ABC.

4. **Fragile path resolution.** Each layer calls `make_paths_absolute()` with a different base path (`os.getcwd()`, `E2G_DIR_PATH`, `ABC_DIR_PATH`). This works because scE2G resolves paths first and subsequent calls are no-ops on absolute paths, but it's not obvious and would break if the call order changed.

5. **No branch pinning in `.gitmodules`.** Neither submodule pins a branch or tag. If ENCODE_rE2G or ABC's default branch changes (e.g., new config keys, renamed functions), scE2G breaks silently.

6. **Tight cross-layer coupling.** scE2G directly calls functions through module namespaces:
   ```python
   encode_e2g.ABC.determine_mem_mb
   encode_e2g.get_feature_table_file(...)
   encode_e2g._get_model_dir_from_wildcards(...)
   ```
   This means scE2G depends on internal implementation details of both submodules.

### Recommendation

The three-repo structure is reasonable given that ENCODE_rE2G and ABC have independent users. But the integration layer needs cleanup:

- **Pin submodule commits or tags** in `.gitmodules` (or document which commits are compatible).
- **Standardize config key naming** across layers, or create a single schema definition.
- **Extract shared utilities** (`make_paths_absolute`, `determine_mem_mb`) into a small shared module rather than copy-pasting.
- **Reduce cross-layer function calls** — scE2G reaching into `encode_e2g.ABC.*` is fragile. Wrap needed functionality in scE2G-level functions.

An alternative worth considering: **flatten ENCODE_rE2G into scE2G** and keep only ABC as a submodule. scE2G already overrides 5 ENCODE_rE2G rules and adds significant custom logic. The ENCODE_rE2G layer mostly provides feature engineering and model application, which could be `include`d directly. This would eliminate one config transformation layer entirely. The tradeoff is that standalone ENCODE_rE2G users would need to use that repo independently.

---

## 2. Computational bottlenecks

### Primary bottleneck: Kendall correlation

- **Rule:** `compute_kendall` (`workflow/rules/compute_kendall.smk`)
- **Resources:** 63 GB minimum memory, 12-36 hour runtime
- **Problem:** Loads full ATAC and RNA matrices into memory, computes Kendall tau-b for all enhancer-gene pairs genome-wide in a single job. The R script uses Rcpp for the inner loop but processes all chromosomes together.
- **Opportunity:** Split by chromosome. 22 parallel jobs at ~3-6 GB each would reduce wall time from 12-36 hours to the runtime of the slowest chromosome (~1-2 hours), and reduce per-job memory by ~10x.

### Secondary bottleneck: ABC predictions

- **Rule:** `abc_create_predictions` (ABC's `predict.py`)
- **Problem:** Loops over chromosomes sequentially within a single Snakemake job. No rule-level parallelism.
- **Opportunity:** Expose per-chromosome ABC prediction as separate Snakemake rules. This is an ABC-level change.

### Tertiary bottleneck: fragment processing

- **Rules:** `frag_to_tagAlign`, `frag_to_bigWig`, `frag_to_norm_bigWig`
- **Resources:** 8-16 threads, 250 GB memory cap, 24-48 hour runtime
- **Problem:** Fragment files are read and sorted multiple times independently. The tagAlign conversion and both bigwig rules each sort the same underlying data.
- **Opportunity:** Share a single sorted intermediate across these rules. The bigwig rules already use `temp()` for bedgraph files, which is good.

### Feature generation

- **Rules in** `genomewide_features.smk`: `activity_only_features`, `add_external_features`, `gen_final_features`
- **Problem:** Process genome-wide prediction files (potentially 20-100M rows) as single monolithic jobs with 16-32 GB memory.
- **Opportunity:** Per-chromosome splitting would reduce memory and enable parallelism, but this requires changes in ENCODE_rE2G.

---

## 3. Memory usage

The global `MAX_MEM_MB = 250000` (250 GB) is used as a cap by `determine_mem_mb()`, which scales memory as 4x input file size (8x for gzipped). This is reasonable as a safety cap, but:

- The Kendall rule's 63 GB minimum is the only rule that genuinely needs high memory in its current form.
- Most other rules would function fine with 8-32 GB.
- The dynamic allocation function is sensible, but it's duplicated across layers.

---

## 4. I/O patterns

**Repeated file reads:**
- RNA count matrix is read independently by `generate_atac_matrix.R`, `compute_kendall.R`, and expression metric scripts. Each performs its own parsing of potentially large H5AD/CSV files.
- ATAC fragments are read 4+ times: fragment counting, cell barcode extraction, tagAlign conversion, bigwig generation, normalized bigwig generation.
- ABC's `EnhancerPredictionsAllPutative.tsv.gz` is read by Kendall pair creation, ARC-E2G, feature generation, and model application.

**Intermediate file sizes:**
- `atac_matrix.rds` (Seurat object): 0.5-5 GB per cluster, persists after pipeline completes
- `genomewide_features.tsv.gz`: can be very large for whole-genome runs
- Temporary bedgraph files correctly use `temp()` markers

**Suggestion:** Mark large intermediates like `atac_matrix.rds` and `kendall_pairs.tsv.gz` as `temp()` if they're only consumed downstream, or document which outputs users should expect to keep.

---

## 5. Rule organization

**scE2G layer** (`workflow/rules/`): 8 `.smk` files, well-organized by pipeline stage. Each file is focused and small (30-160 lines).

**ENCODE_rE2G layer** (`ENCODE_rE2G/workflow/rules/`): 8 `.smk` files covering features, predictions, QC, training, and model comparison. Clean separation.

**ABC layer** (`ENCODE_rE2G/ABC/workflow/rules/`): 6 `.smk` files for peaks, candidates, neighborhoods, predictions, QC. Straightforward.

The rule organization itself is clean. The main structural issue is the 5 rules scE2G excludes from ENCODE_rE2G and reimplements:
- `generate_e2g_predictions` → replaced by `run_e2g_qnorm` (adds quantile normalization)
- `format_external_features_config` → replaced by scE2G version (adds Kendall/ARC features)
- `get_stats`, `generate_plots` → replaced with scE2G-specific QC
- `all` → replaced with scE2G's own target rule

This override pattern is the right approach, but the excluded rules and their replacements should be documented explicitly (e.g., a comment at the `use rule` line explaining why each is excluded).

---

## 6. Version pinning

| Dependency | scE2G env | ENCODE_rE2G env | ABC env |
|---|---|---|---|
| bedtools | 2.29.2 | 2.29.2 | **2.26.0** |
| scikit-learn | 1.2.1 | 1.2.1 | — |
| python | >=3.6 | >=3.6 | 3.10.14 |
| snakemake | >=7.0 | >=7.0 | 7.32.4 |

The bedtools version mismatch (2.29.2 vs 2.26.0) is a concern. Since each layer uses `--use-conda` with its own environment, rules execute in the correct environment for their layer. But this assumes Snakemake correctly resolves which conda env to use for module-imported rules — worth verifying.

---

## 7. Summary of recommendations (prioritized)

### High impact
1. **Split Kendall correlation by chromosome** — largest single performance win
2. **Pin submodule commits** — prevent silent breakage
3. **Standardize config key naming** across layers

### Medium impact
4. **Share sorted fragment intermediate** across tagAlign/bigwig rules
5. **Document rule exclusion/override rationale** at the `use rule` line
6. **Extract shared utility functions** to avoid triple-duplication
7. **Mark large intermediates as `temp()`** where appropriate

### Lower priority
8. Consider flattening ENCODE_rE2G into scE2G (major refactor, weigh against standalone ENCODE_rE2G usage)
9. Per-chromosome ABC prediction splitting (requires ABC-level changes)
10. Verify bedtools version isolation across module-imported rules
