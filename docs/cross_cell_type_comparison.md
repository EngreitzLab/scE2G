# Cross-cell-type comparison workflow

Tools for comparing E2G predictions across multiple cell types. All scripts live in
`workflow/scripts/visualization/` and use the pixi environment defined in
`workflow/scripts/visualization/pixi.toml`.

## Setup

```bash
module load pixi
pixi install --manifest-path workflow/scripts/visualization/pixi.toml
```

Run any script via:
```bash
pixi run --manifest-path workflow/scripts/visualization/pixi.toml Rscript workflow/scripts/visualization/<script>.R [args]
```

---

## Step 1: Precompute comparison tables

Run once per set of CI-annotated results. All downstream visualization scripts read
from its outputs.

```bash
pixi run --manifest-path workflow/scripts/visualization/pixi.toml \
  Rscript workflow/scripts/visualization/precompute_comparisons.R \
    --ci_dir results/2026_0210_WTC11_EC \
    --model multiome_powerlaw_v3 \
    --threshold 0.177 \
    --output_dir results/2026_0210_WTC11_EC/comparison
```

**Key parameters:**
- `--ci_dir`: directory containing per-cell-type subdirs (each with `{model}/encode_e2g_predictions_ci_annotated.tsv.gz` or `encode_e2g_predictions.tsv.gz`)
- `--model`: model subdirectory name (default: `multiome_powerlaw_v3`)
- `--threshold`: score threshold for above/below calls (default: 0.177)
- `--cell_types`: comma-separated list to restrict analysis; auto-discovers from `ci_dir` if omitted
- `--merge_gap`: gap in bp for consensus element merging (default: 0 = overlapping only; try 500 for neighborhood-level analysis)
- `--output_dir`: where to write all output tables

**Outputs:**
| File | Description |
|------|-------------|
| `consensus_elements.tsv` | Consensus element coordinates + `n_cell_types_with_element` |
| `consensus_scores.tsv.gz` | Long-format score matrix: one row per (consensus_id, TargetGene, cell_type) |
| `link_specificity.tsv` | Per-link: `n_cell_types_above`, `specificity_score` (1 − pct), `cell_types_above` |
| `gene_summary.tsv` | Per (gene, cell_type): `n_links_above_threshold`, `max_score`, `sum_score`, `RNA_pseudobulkTPM`, `normalizedATAC_prom` |
| `pairwise_metrics.tsv` | Per cell-type pair: `jaccard`, `Pearson`, `Pearson_log1p`, `Spearman` |

**Consensus element definition:** all elements across all cell types are pooled and
merged with `GenomicRanges::reduce()`. Shared E-G links are defined as matching on
(consensus_id, TargetGene) — slightly more permissive than 50% reciprocal overlap but
consistent and transitive. Self-promoters are excluded from all comparison metrics but
retained in `consensus_scores.tsv.gz`.

---

## Step 2a: Global pairwise similarity heatmap

Double-triangle heatmap of pairwise cell-type similarity. Lower triangle = Pearson
correlation of scores (blue); upper triangle = Jaccard overlap of above-threshold links
(red). Cell types are clustered by average-linkage on Jaccard × Pearson product.

```bash
pixi run --manifest-path workflow/scripts/visualization/pixi.toml \
  Rscript workflow/scripts/visualization/plot_global_comparison.R \
    --pairwise_metrics results/2026_0210_WTC11_EC/comparison/pairwise_metrics.tsv \
    --output results/2026_0210_WTC11_EC/comparison/global_comparison.pdf
```

**Key parameters:**
- `--lower_metric`: metric for lower triangle (default: `Pearson`; options: `Pearson`, `Pearson_log1p`, `Spearman`)
- `--upper_metric`: metric for upper triangle (default: `jaccard`)
- `--lower_limits` / `--upper_limits`: color scale limits as `min,max` (default: `0,0.6`)
- `--sample_key`: optional TSV with columns `biosample`, `display_name`, `dataset` — adds display names and a dataset color bar alongside the heatmap
- `--width` / `--height`: output dimensions in inches (defaults: 7 × 6)

---

## Step 2b: Single-gene enhancer landscape comparison

Bubble matrix showing the enhancer landscape for one gene across all cell types.
Rows = consensus elements (ordered by genomic position); columns = cell types.
Dot size = ATAC accessibility; dot color = E2G score (purple gradient); grey = below
threshold or not called.

```bash
pixi run --manifest-path workflow/scripts/visualization/pixi.toml \
  Rscript workflow/scripts/visualization/plot_gene_comparison.R \
    --gene GATA2 \
    --comparison_dir results/2026_0210_WTC11_EC/comparison \
    --threshold 0.177 \
    --cell_types d0_ipsc,d1_ps,d2_meso,d3_ec,d4_ec \
    --output results/2026_0210_WTC11_EC/comparison/gene_plots/GATA2.pdf
```

**Key parameters:**
- `--gene`: target gene name (required)
- `--comparison_dir`: directory containing `consensus_scores.tsv.gz` from step 1
- `--cell_types`: comma-separated ordered cell types (default: alphabetical)
- `--max_elements`: maximum elements to show; keeps highest-scoring (default: 75)
- `--output`: output PDF path (default: `{comparison_dir}/gene_plots/{gene}.pdf`)
- `--width` / `--height`: output dimensions in inches (auto-sized if omitted)

A top strip shows RNA TPM and promoter ATAC signal per cell type.

---

## Step 2c: Expression vs. enhancer landscape changes

Scatter plots relating gene expression changes to enhancer landscape changes between
two cell types. One point per gene; colored by expression direction (up/stable/down).

Three panels:
1. log2FC(TPM) vs delta max E2G score
2. log2FC(TPM) vs delta links above threshold
3. log2FC(TPM) vs delta promoter ATAC signal

```bash
pixi run --manifest-path workflow/scripts/visualization/pixi.toml \
  Rscript workflow/scripts/visualization/plot_expression_vs_enhancers.R \
    --gene_summary results/2026_0210_WTC11_EC/comparison/gene_summary.tsv \
    --cell_type_A d0_ipsc \
    --cell_type_B d4_ec \
    --output results/2026_0210_WTC11_EC/comparison/expression_vs_enhancers_d0_vs_d4.pdf
```

**Key parameters:**
- `--cell_type_A` / `--cell_type_B`: reference and comparison cell type
- `--min_tpm`: minimum TPM in at least one cell type to include gene (default: 1.0)
- `--log2fc_threshold`: threshold to call a gene up or down (default: 1.0)
- `--tpm_pseudocount`: pseudocount before log2FC (default: 0.1)
- `--highlight_genes`: comma-separated gene names to label on plots
- `--width` / `--height`: output dimensions in inches (default: 10 × 4)

Pearson correlations between log2FC(TPM) and each enhancer metric are printed to stdout.

---

## Full example (WTC11 → EC differentiation)

```bash
CMP=results/2026_0210_WTC11_EC/comparison
PIXI="pixi run --manifest-path workflow/scripts/visualization/pixi.toml"
CT=d0_ipsc,d1_ps,d2_meso,d3_ec,d4_ec

# Step 1: precompute (run once)
$PIXI Rscript workflow/scripts/visualization/precompute_comparisons.R \
  --ci_dir results/2026_0210_WTC11_EC --model multiome_powerlaw_v3 \
  --threshold 0.177 --output_dir $CMP

# Step 2a: global heatmap
$PIXI Rscript workflow/scripts/visualization/plot_global_comparison.R \
  --pairwise_metrics $CMP/pairwise_metrics.tsv \
  --output $CMP/global_comparison.pdf

# Step 2b: gene bubble matrix (repeat for any gene of interest)
$PIXI Rscript workflow/scripts/visualization/plot_gene_comparison.R \
  --gene GATA2 --comparison_dir $CMP --threshold 0.177 --cell_types $CT

# Step 2c: expression vs enhancers (repeat for any pair of cell types)
$PIXI Rscript workflow/scripts/visualization/plot_expression_vs_enhancers.R \
  --gene_summary $CMP/gene_summary.tsv \
  --cell_type_A d0_ipsc --cell_type_B d4_ec \
  --output $CMP/expression_vs_enhancers_d0_vs_d4.pdf
```
