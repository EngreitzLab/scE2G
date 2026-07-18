# Plan: Flatten ENCODE_rE2G submodule into scE2G

## Goal

Eliminate the nested submodule architecture by merging ENCODE_rE2G directly into scE2G, keeping only ABC as a direct submodule.

**Current:** `scE2G/ENCODE_rE2G/ABC/` (nested)
**Target:** `scE2G/ABC/` (direct)

## Rationale

- Eliminates one layer of config transformation (`get_e2g_config()` → `get_abc_config()` becomes just `get_abc_config()`)
- Removes fragile cross-layer function calls (`encode_e2g.ABC.*`)
- Simplifies maintenance (one less repo to coordinate)
- ENCODE_rE2G rules are already heavily customized by scE2G (5 rules excluded/replaced)

## Files to modify

### 1. `workflow/Snakefile`

**Changes:**
- Replace `module encode_e2g` import with direct `module ABC` import
- Replace `use rule * from encode_e2g exclude ...` with `use rule * from ABC exclude all as abc_*`
- Add `include:` statements for merged ENCODE_rE2G rule files
- Update function calls:
  - `encode_e2g.ABC.determine_mem_mb` → `ABC.determine_mem_mb`
  - `encode_e2g.ABC.enable_retry` → `ABC.enable_retry`
  - `encode_e2g.expand_biosample_df` → `expand_biosample_df` (local)
- Remove `encode_re2g_dir` references, add `abc_dir`

### 2. `workflow/rules/utils.smk`

**Add functions from ENCODE_rE2G:**
- `expand_biosample_df()` (lines 98-120)
- `get_feature_table_file()` (lines 186-188)
- `get_trained_model()` (lines 190-192)
- `get_model_threshold()` (lines 194-199)
- `get_tpm_threshold()` (lines 201-209)
- `_get_model_dir_from_wildcards()` (lines 151-157)
- `_get_biosample_model_dir_from_row()` (lines 159-184)
- `_validate_model_dir()` (lines 140-149)
- `make_accessibility_file_df()` (lines 70-92)

**Replace:**
- `get_e2g_config()` → new `get_abc_config()` that maps scE2G config directly to ABC config

### 3. `workflow/rules/sc_predictions.smk`

**Update 12 function references:**
- Lines 48-51: `encode_e2g.get_feature_table_file`, `encode_e2g.get_trained_model`, `encode_e2g._get_model_dir_from_wildcards`, `encode_e2g.get_tpm_threshold`
- Lines 57, 91, 133, 157: `encode_e2g.ABC.determine_mem_mb`

### 4. `workflow/rules/save_configs.smk`

**Update:**
- Line 9: `encode_e2g.get_abc_config(encode_e2g_config)` → `get_abc_config(config)`

### 5. `config/config.yaml`

**Remove:**
```yaml
encode_re2g_dir: "ENCODE_rE2G"
```

**Add:**
```yaml
abc_dir: "ABC"
```

## Files to create (copy from ENCODE_rE2G)

### Rule files → `workflow/rules/`

| Source | Destination | Notes |
|--------|-------------|-------|
| `ENCODE_rE2G/workflow/rules/genomewide_features.smk` | `workflow/rules/e2g_genomewide_features.smk` | Feature generation rules |
| `ENCODE_rE2G/workflow/rules/predictions.smk` | `workflow/rules/e2g_predictions.smk` | Remove `generate_e2g_predictions`, `format_external_features_config` (scE2G overrides) |

### Scripts → `workflow/scripts/`

```
ENCODE_rE2G/workflow/scripts/feature_tables/* → workflow/scripts/feature_tables/
ENCODE_rE2G/workflow/scripts/model_application/* → workflow/scripts/model_application/
```

### Environment

```
ENCODE_rE2G/workflow/envs/encode_re2g.yml → workflow/envs/encode_re2g.yml
```

### Resources

```
ENCODE_rE2G/resources/external_features/gene_promoter_class*.tsv → resources/external_features/
```

## Git operations

### Step 1: Copy files before removing submodule

```bash
# Copy rule files
cp ENCODE_rE2G/workflow/rules/genomewide_features.smk workflow/rules/e2g_genomewide_features.smk
cp ENCODE_rE2G/workflow/rules/predictions.smk workflow/rules/e2g_predictions.smk

# Copy scripts
cp -r ENCODE_rE2G/workflow/scripts/feature_tables/* workflow/scripts/feature_tables/
cp -r ENCODE_rE2G/workflow/scripts/model_application/* workflow/scripts/model_application/

# Copy env
cp ENCODE_rE2G/workflow/envs/encode_re2g.yml workflow/envs/

# Copy resources
cp ENCODE_rE2G/resources/external_features/gene_promoter_class*.tsv resources/external_features/
```

### Step 2: Note ABC commit hash

```bash
cd ENCODE_rE2G/ABC && git rev-parse HEAD
# Save this hash for step 4
```

### Step 3: Remove ENCODE_rE2G submodule

```bash
git submodule deinit ENCODE_rE2G
git rm ENCODE_rE2G
rm -rf .git/modules/ENCODE_rE2G
```

### Step 4: Add ABC as direct submodule

```bash
git submodule add https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction.git ABC
cd ABC && git checkout <commit-hash-from-step-2>
cd .. && git add ABC
```

## Key config transformation change

**Current chain:**
```
scE2G config
    ↓ get_e2g_config()
ENCODE_rE2G config
    ↓ get_abc_config()
ABC config
```

**Target:**
```
scE2G config
    ↓ get_abc_config()
ABC config
```

The new `get_abc_config()` must handle the full mapping:
- `ABC_BIOSAMPLES` → `biosamplesTable`
- `gene_TSS500` → `ref.genome_tss`
- `genes` → `ref.genes`
- `chr_sizes` → `ref.chrom_sizes`
- `regions_blocklist` → `ref.regions_blocklist`
- `macs2_genomesize` → `params_macs.genome_size`

## Rules to handle

### Rules scE2G keeps from ENCODE_rE2G (include in copied files)

- `generate_num_candidate_enh_gene`
- `generate_num_tss_enh_gene`
- `generate_num_sum_enhancers`
- `activity_only_features`
- `add_external_features`
- `gen_final_features`
- `make_biosample_feature_table`
- `filter_e2g_predictions`
- `write_predictions_bedpe`

### Rules scE2G overrides (remove from copied files)

- `generate_e2g_predictions` → `run_e2g_qnorm`
- `format_external_features_config` → `make_external_features_config`
- `get_stats` → `get_stats_per_model_per_cluster`
- `generate_plots` → `plot_stats` + `hover_plots`
- `all` → scE2G's own target

## Verification

### 1. Dry run

```bash
snakemake -n --configfile config/config.yaml
```

### 2. Run tests

```bash
pytest -s tests/test_sce2g_apply.py
```

### 3. Full pipeline test

```bash
snakemake -j4 --use-conda --configfile config/config.yaml
```

### 4. Compare outputs

Compare `scE2G_predictions.tsv.gz` before and after - should be identical.

## Risks and mitigations

| Risk | Mitigation |
|------|------------|
| Function references to globals like `BIOSAMPLE_DF` | Ensure globals defined before `include:` statements |
| Conda env path issues | Verify paths in copied env file |
| Rule name conflicts | Delete overridden rules from copied files |
| Missing script dependencies | Run dry-run to catch import errors |

## Implementation order

1. Copy all files from ENCODE_rE2G (before removal)
2. Update `workflow/rules/utils.smk` with consolidated functions
3. Edit copied rule files to remove overridden rules
4. Update `workflow/Snakefile` (module import, includes, function refs)
5. Update `workflow/rules/sc_predictions.smk` function references
6. Update `config/config.yaml`
7. Git operations (remove ENCODE_rE2G, add ABC)
8. Test
