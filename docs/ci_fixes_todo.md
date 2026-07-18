# CI pipeline fixes

Tracking updates and fixes for the confidence interval workflow after first real dataset run.

## To do

10. [ ] QC report feedback from rerun (awaiting user input)

### Downstream analysis scripts

12. [ ] Write script/rule to identify cell-type-specific E-G pairs using CI information
    - Purpose: Find E-G pairs specific to a cell type or group of cell types
    - Approach options:
      a. Non-overlapping CIs between cell types
      b. Score in target cell type significantly higher than others
      c. Above threshold in target, below in others (with CI consideration)
    - Output: Table of specific pairs with statistical support

### Deferred

13. [ ] Fix race condition in main Snakefile (defer until all other changes complete)

## Completed

1. [x] Confirm subsample scores come from E2G.Score.qnorm column
   - Verified at `compute_ci_annotations.R:92`
2. [x] Add comma-separated list of individual subsample scores column to CI annotations output
   - Added `E2G.Score.subsample_scores` column in `compute_ci_annotations.R:127`
   - Format: comma-separated scores rounded to 6 decimal places
3. [x] Consolidate report into single document covering all cell types
   - Changed `Snakefile_ci` rule all to request single `ci_qc_report.html`
   - Updated `ci_analysis.smk` plot_ci rule to aggregate all clusters
   - Rewrote `plot_ci.R` to handle multi-cluster input
4. [x] Section 1: Basic parameters table
   - Labels: "Cell types:", "Model name:", "Number of subsamples:", "Fraction of cells sampled:"
5. [x] Section 2: Subsample sizes per cell type
   - 3 bar plots (cells, fragments, UMIs) with separate y-axes via patchwork
   - Darker bar = full cluster, lighter bar = average subsample
   - Dashed red lines at y=100 (cells), y=1M (fragments), y=2M (UMIs)
   - Log-scale y-axis
   - Gracefully hides UMI plot if model contains "scATAC"
   - Added depth summary table
6. [x] Section 3: CI width distribution
   - Density plots instead of histograms (fixes y-axis scaling)
   - Two-column facet: all pairs (left), above threshold (right)
   - One plot per cell type with n= in subtitle
   - Pagination at ~6 cell types per page
7. [x] Section 4: Score stability
   - 2D histogram (geom_bin2d) with log-scaled color, 3 square plots per row
   - Quantile ribbon plots with adaptive binning (50 bins in 0.1 range around threshold, 100 bins elsewhere)
   - Ribbons show 50%, 80%, 95% intervals; median as solid line
   - Identity line and threshold lines overlaid
   - Pagination at 3 plots per row
8. [x] Section 5: Prediction stability
   - Overlaid lines (one per cell type, colored) showing % subsamples above threshold vs score
   - Side-by-side: full range (left) + zoomed around threshold ±0.1 (right)
   - Point size indicates bin sample size
   - Shared legend at bottom
9. [x] Section 6: Element and gene coverage
   - Element coverage: unique elements, % y-axis, grouped bars per cell type
   - Gene coverage now bidirectional (multiome only):
     a. Genes in full run: how many subsamples have them above threshold
     b. Genes only in subsamples: genes above threshold in subsamples but not full run
   - Uses "above threshold" as proxy for expressed (captures all expression-based filtering)
   - Added gene_coverage_stats.tsv output from compute_ci_annotations.R
10. [x] Write script to plot predictions with CIs for a given locus/gene across cell types
   - Created `workflow/scripts/visualization/plot_locus_ci.R`
   - CLI args: --locus, --genes, --ci_dir, --output_prefix, --cell_types, --model, --threshold, --order_by
   - Outputs: `*_ci_bars.pdf` (grouped bar + error bars), `*_ci_dots.pdf` (dot plot), `*_arcs.pdf` (arc plot), `*_summary.tsv`
   - Auto-discovers cell types from CI results directory
   - Nature color palette integration
