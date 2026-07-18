#!/usr/bin/env Rscript
#
# rank_variable_genes.R
#
# Ranks genes by variability of enhancer landscape across cell types using three metrics:
#
#   Option 2 (score variance): for each (gene, element), compute variance of E2G scores
#     across cell types (absent elements treated as score=0). Aggregate per gene as
#     mean variance across elements.
#
#   Option 3 (n_links CV): coefficient of variation of n_links_above_threshold
#     across cell types per gene.
#
#   Option 4 (CI-weighted score variance): same as option 2 but each score is weighted
#     by E2G.Score.pct_above_threshold before computing variance, so elements with
#     unstable CI scores contribute less. Requires CI columns in consensus_scores.tsv.gz;
#     skipped gracefully if absent.
#
#   Combined rank: equal-weight average of available option ranks.
#
# Usage:
#   Rscript rank_variable_genes.R \
#     --comparison_dir results/.../comparison \
#     --output results/.../comparison/variable_genes.tsv

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
})

option_list <- list(
    make_option("--comparison_dir", type="character", help="Directory with precompute outputs"),
    make_option("--min_cell_types", type="integer", default=2L,
                help="Min cell types a gene must appear in [default: %default]"),
    make_option("--output", type="character", default=NULL,
                help="Output TSV [default: {comparison_dir}/variable_genes.tsv]")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$comparison_dir)) stop("--comparison_dir is required")
if (is.null(opt$output)) opt$output <- file.path(opt$comparison_dir, "variable_genes.tsv")

## ============================================================================
## Option 3: n_links variance from gene_summary.tsv
## ============================================================================

message("Option 3: loading gene_summary.tsv...")

gene_summary <- fread(file.path(opt$comparison_dir, "gene_summary.tsv"))
all_cell_types <- unique(gene_summary$cell_type)
n_ct_total <- length(all_cell_types)

n_links_stats <- gene_summary[, .(
    n_cell_types   = .N,
    n_links_mean   = mean(n_links_above_threshold, na.rm=TRUE),
    n_links_sd     = sd(n_links_above_threshold,   na.rm=TRUE),
    n_links_cv     = {
        m <- mean(n_links_above_threshold, na.rm=TRUE)
        s <- sd(n_links_above_threshold,   na.rm=TRUE)
        if (!is.na(m) && m > 0) s / m else NA_real_
    },
    n_links_range  = max(n_links_above_threshold, na.rm=TRUE) - min(n_links_above_threshold, na.rm=TRUE),
    n_links_max    = max(n_links_above_threshold, na.rm=TRUE),
    mean_TPM       = mean(RNA_pseudobulkTPM, na.rm=TRUE)
), by=TargetGene]

message("  ", nrow(n_links_stats), " genes")

## ============================================================================
## Option 2: score variance from consensus_scores.tsv.gz
## ============================================================================

message("Option 2 / Option 4: loading consensus_scores.tsv.gz (may take a moment)...")

scores_path <- file.path(opt$comparison_dir, "consensus_scores.tsv.gz")
ci_col      <- "E2G.Score.pct_above_threshold"
header      <- names(fread(scores_path, nrows=0))
has_ci      <- ci_col %in% header
load_cols   <- c("consensus_id", "TargetGene", "cell_type", "score", "isSelfPromoter",
                 if (has_ci) ci_col)

scores <- fread(scores_path, select=load_cols)

if ("isSelfPromoter" %in% names(scores)) {
    scores <- scores[isSelfPromoter != TRUE]
    scores[, isSelfPromoter := NULL]
}

if (has_ci) {
    message("  CI column found: will compute CI-weighted score variance (option 4)")
} else {
    message("  No CI column found; skipping CI-weighted variance")
}

## For each (gene, element): fill absent cell types with score=0, then compute variance
all_keys <- unique(scores[, .(consensus_id, TargetGene)])
full_grid <- all_keys[, .(cell_type = all_cell_types), by=.(consensus_id, TargetGene)]
full_grid  <- merge(full_grid, scores, by=c("consensus_id", "TargetGene", "cell_type"), all.x=TRUE)
full_grid[is.na(score), score := 0]
if (has_ci) full_grid[is.na(get(ci_col)), (ci_col) := 0]

## Option 2: raw score variance
element_var <- full_grid[, .(element_score_var = var(score, na.rm=TRUE)),
                          by=.(TargetGene, consensus_id)]

score_var_stats <- element_var[, .(
    score_var_mean = mean(element_score_var, na.rm=TRUE),
    score_var_max  = max(element_score_var,  na.rm=TRUE),
    n_elements     = .N
), by=TargetGene]

message("  ", nrow(score_var_stats), " genes")

## Option 4: CI-weighted score variance (score * pct_above_threshold; absent = 0 for both)
if (has_ci) {
    element_var_ci <- full_grid[, .(
        element_score_var_ci = var(score * get(ci_col), na.rm=TRUE)
    ), by=.(TargetGene, consensus_id)]

    score_var_ci_stats <- element_var_ci[, .(
        score_var_ci_mean = mean(element_score_var_ci, na.rm=TRUE),
        score_var_ci_max  = max(element_score_var_ci,  na.rm=TRUE)
    ), by=TargetGene]
}

## ============================================================================
## Combine and rank
## ============================================================================

message("Combining and ranking...")

combined <- merge(n_links_stats, score_var_stats, by="TargetGene", all=TRUE)
if (has_ci) combined <- merge(combined, score_var_ci_stats, by="TargetGene", all=TRUE)
combined  <- combined[n_cell_types >= opt$min_cell_types]

combined[, rank_n_links_cv  := frankv(-n_links_cv,    ties.method="average", na.last=TRUE)]
combined[, rank_score_var   := frankv(-score_var_mean, ties.method="average", na.last=TRUE)]
if (has_ci) {
    combined[, rank_score_var_ci := frankv(-score_var_ci_mean, ties.method="average", na.last=TRUE)]
    combined[, combined_rank     := (rank_n_links_cv + rank_score_var + rank_score_var_ci) / 3]
} else {
    combined[, combined_rank := (rank_n_links_cv + rank_score_var) / 2]
}
setorder(combined, combined_rank)

fwrite(combined, opt$output, sep="\t")
message("Wrote ", nrow(combined), " genes to: ", opt$output)

message("\nTop 10 most variable genes:")
show_cols <- c("TargetGene", "n_links_cv", "n_links_range", "score_var_mean",
               if (has_ci) "score_var_ci_mean", "combined_rank")
print(combined[seq_len(min(10L, .N)), ..show_cols])
