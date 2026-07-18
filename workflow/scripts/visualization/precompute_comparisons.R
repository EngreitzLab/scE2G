#!/usr/bin/env Rscript
#
# precompute_comparisons.R
#
# Builds a consensus element catalog from all cell types' predictions and derives
# summary tables used by all cross-cell-type visualization scripts. Run once;
# visualization scripts read from the outputs.
#
# Usage:
#   Rscript precompute_comparisons.R \
#     --ci_dir /path/to/ci_results \
#     --model  multiome_powerlaw_v3 \
#     --threshold 0.177 \
#     --output_dir /path/to/output
#
# Outputs (all in --output_dir):
#   consensus_elements.tsv    consensus element coordinates + n_cell_types_with_element
#   consensus_scores.tsv.gz   long-format score matrix: (consensus_id, TargetGene, cell_type)
#   link_specificity.tsv      per-link: how many cell types call it above threshold
#   gene_summary.tsv          per (gene, cell type): n_links, max_score, TPM
#   pairwise_metrics.tsv      per cell type pair: Jaccard, Pearson, Spearman
#
# Consensus element definition:
#   All elements across all cell types are pooled and merged with GenomicRanges::reduce()
#   using --merge_gap (default 0 = overlapping-only merge). Each original element maps
#   to exactly one consensus element. Shared links are defined as matching on
#   (consensus_id, TargetGene) — slightly more permissive than strict 50% reciprocal
#   overlap but consistent and transitive.

suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(optparse)
})

## ============================================================================
## CLI
## ============================================================================

option_list <- list(
    make_option("--ci_dir",     type="character", help="CI results directory"),
    make_option("--model",      type="character", default="multiome_powerlaw_v3",
                help="Model name [default: %default]"),
    make_option("--threshold",  type="numeric",   default=0.177,
                help="Score threshold for above/below calls [default: %default]"),
    make_option("--cell_types", type="character", default=NULL,
                help="Comma-separated cell types (default: auto-discover from ci_dir)"),
    make_option("--score_col",  type="character", default="E2G.Score.qnorm",
                help="Score column in prediction files [default: %default]"),
    make_option("--merge_gap",  type="integer",   default=0L,
                help="Gap in bp for consensus element merging; 0 = overlapping only [default: %default]"),
    make_option("--output_dir", type="character", help="Output directory")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$ci_dir))     stop("--ci_dir is required")
if (is.null(opt$output_dir)) stop("--output_dir is required")
dir.create(opt$output_dir, recursive=TRUE, showWarnings=FALSE)

## ============================================================================
## Discover / validate cell types
## ============================================================================

discover_cell_types <- function(ci_dir, model) {
    subdirs <- list.dirs(ci_dir, recursive=FALSE, full.names=FALSE)
    valid <- subdirs[sapply(subdirs, function(d) {
        any(file.exists(
            file.path(ci_dir, d, model, "encode_e2g_predictions_ci_annotated.tsv.gz"),
            file.path(ci_dir, d, model, "encode_e2g_predictions.tsv.gz")
        ))
    })]
    if (length(valid) == 0) stop("No prediction files found under: ", ci_dir)
    message("Discovered ", length(valid), " cell types: ", paste(valid, collapse=", "))
    valid
}

cell_types <- if (!is.null(opt$cell_types)) {
    trimws(strsplit(opt$cell_types, ",")[[1]])
} else {
    discover_cell_types(opt$ci_dir, opt$model)
}
n_ct <- length(cell_types)

get_pred_path <- function(ct) {
    ci_path   <- file.path(opt$ci_dir, ct, opt$model, "encode_e2g_predictions_ci_annotated.tsv.gz")
    full_path <- file.path(opt$ci_dir, ct, opt$model, "encode_e2g_predictions.tsv.gz")
    if (file.exists(ci_path)) ci_path else full_path
}

## ============================================================================
## Load predictions
## ============================================================================

CI_COLS   <- c("E2G.Score.CI_low", "E2G.Score.CI_high", "E2G.Score.pct_above_threshold")
BASE_COLS <- c("chr", "start", "end", "TargetGene", "TargetGeneTSS",
               "normalizedATAC_enh", "normalizedATAC_prom",
               "RNA_pseudobulkTPM", "isSelfPromoter")

message("Loading predictions...")
preds_list <- lapply(cell_types, function(ct) {
    path <- get_pred_path(ct)
    if (!file.exists(path)) {
        warning("Prediction file not found for ", ct, ": ", path)
        return(NULL)
    }
    avail <- names(fread(path, nrows=0))
    missing_req <- setdiff(c(BASE_COLS, opt$score_col), avail)
    if (length(missing_req) > 0)
        warning("Missing columns in ", ct, ": ", paste(missing_req, collapse=", "))
    cols <- intersect(c(BASE_COLS, opt$score_col, CI_COLS), avail)
    dt   <- fread(path, select=cols)
    setnames(dt, opt$score_col, "score")
    dt[, cell_type := ct]
    dt
})
preds_list <- Filter(Negate(is.null), preds_list)
n_loaded   <- length(preds_list)
message("  Loaded ", n_loaded, "/", n_ct, " cell types")
if (n_loaded < 2) stop("Need at least 2 cell types for pairwise comparison")

has_ci <- all(CI_COLS %in% names(preds_list[[1]]))
if (!has_ci) message("  CI columns absent — CI fields will be omitted from consensus_scores")

## ============================================================================
## Step 1: Consensus element catalog
## ============================================================================

message("\nStep 1: Building consensus element catalog (merge_gap=", opt$merge_gap, "bp)...")

all_elements <- unique(rbindlist(lapply(preds_list, function(dt) dt[, .(chr, start, end)])))
setorder(all_elements, chr, start, end)

all_gr       <- makeGRangesFromDataFrame(all_elements)
consensus_gr <- reduce(all_gr, min.gapwidth = opt$merge_gap + 1L)

consensus_dt <- data.table(
    consensus_chr   = as.character(seqnames(consensus_gr)),
    consensus_start = start(consensus_gr),
    consensus_end   = end(consensus_gr)
)
consensus_dt[, consensus_id := paste0(consensus_chr, ":", consensus_start, "-", consensus_end)]

# Map each original element to its unique consensus element.
# After reduce(), every original element is contained within exactly one consensus
# element, so findOverlaps gives a clean 1-to-1 mapping.
hits <- findOverlaps(all_gr, consensus_gr, select="first")
all_elements[, consensus_id := consensus_dt$consensus_id[hits]]

# Build keyed lookup: (chr, start, end) -> consensus_id
elem_lookup <- all_elements[, .(chr, start, end, consensus_id)]
setkeyv(elem_lookup, c("chr", "start", "end"))

# Count how many cell types have an element overlapping each consensus element
ct_elem <- unique(rbindlist(lapply(preds_list, function(dt) dt[, .(chr, start, end, cell_type)])))
ct_elem <- merge(ct_elem, elem_lookup, by=c("chr", "start", "end"))
n_ct_per_consensus <- ct_elem[, .(n_cell_types = uniqueN(cell_type)), by=consensus_id]
consensus_dt <- merge(consensus_dt, n_ct_per_consensus, by="consensus_id", all.x=TRUE)
consensus_dt[is.na(n_cell_types), n_cell_types := 0L]

n_original  <- nrow(all_elements)
n_consensus <- nrow(consensus_dt)
message("  ", n_original, " original elements -> ", n_consensus, " consensus elements")

fwrite(consensus_dt[order(consensus_chr, consensus_start)],
       file.path(opt$output_dir, "consensus_elements.tsv"), sep="\t")
message("  Wrote consensus_elements.tsv")

## ============================================================================
## Step 2: Consensus score matrix
## ============================================================================

message("\nStep 2: Building consensus score matrix...")

keep_cols <- c("consensus_id", "TargetGene", "TargetGeneTSS", "cell_type",
               "score", "normalizedATAC_enh", "normalizedATAC_prom",
               "RNA_pseudobulkTPM", "isSelfPromoter", "is_above_threshold")
if (has_ci) keep_cols <- c(keep_cols, CI_COLS)

score_matrix <- rbindlist(lapply(preds_list, function(dt) {
    # Add consensus_id via keyed join (right join: all rows of dt retained)
    setkeyv(dt, c("chr", "start", "end"))
    dt2 <- elem_lookup[dt]

    dt2[, is_above_threshold := (score >= opt$threshold)]
    # isSelfPromoter comes from Python as "True"/"False" strings; coerce to logical
    if ("isSelfPromoter" %in% names(dt2) && is.character(dt2$isSelfPromoter))
        dt2[, isSelfPromoter := (isSelfPromoter == "True")]

    # When multiple original elements from the same cell type map to the same
    # (consensus_id, TargetGene), keep the highest-scoring one
    dt2 <- dt2[dt2[, .I[which.max(score)], by=.(consensus_id, TargetGene)]$V1]

    dt2[, .SD, .SDcols=intersect(keep_cols, names(dt2))]
}), fill=TRUE)

setkeyv(score_matrix, c("consensus_id", "TargetGene", "cell_type"))

fwrite(score_matrix,
       file.path(opt$output_dir, "consensus_scores.tsv.gz"), sep="\t", compress="gzip")
message("  ", nrow(score_matrix), " rows | ",
        uniqueN(score_matrix$consensus_id), " consensus elements | ",
        uniqueN(score_matrix$TargetGene), " genes")
message("  Wrote consensus_scores.tsv.gz")

## ============================================================================
## Step 3: Per-link specificity
## ============================================================================

message("\nStep 3: Computing per-link specificity...")

link_spec <- score_matrix[is_above_threshold == TRUE & isSelfPromoter == FALSE, .(
    n_cell_types_above   = .N,
    pct_cell_types_above = .N / n_ct,
    specificity_score    = 1 - .N / n_ct,   # 1.0 = unique to one cell type
    cell_types_above     = paste(sort(cell_type), collapse=",")
), by=.(consensus_id, TargetGene)]

link_spec <- merge(link_spec,
                   consensus_dt[, .(consensus_id, consensus_chr, consensus_start, consensus_end)],
                   by="consensus_id")
setorder(link_spec, consensus_chr, consensus_start, TargetGene)

fwrite(link_spec, file.path(opt$output_dir, "link_specificity.tsv"), sep="\t")
message("  ", nrow(link_spec), " links above threshold in >=1 cell type")
message("  Wrote link_specificity.tsv")

## ============================================================================
## Step 4: Per-gene summary
## ============================================================================

message("\nStep 4: Computing per-gene summary...")

gene_summary <- score_matrix[isSelfPromoter == FALSE, .(
    n_links_above_threshold = sum(is_above_threshold, na.rm=TRUE),
    max_score               = max(score, na.rm=TRUE),
    sum_score               = sum(score, na.rm=TRUE),
    mean_score              = mean(score, na.rm=TRUE),
    RNA_pseudobulkTPM       = first(RNA_pseudobulkTPM),    # gene-level, same for all elements
    normalizedATAC_prom     = first(normalizedATAC_prom)   # promoter-level, same for all elements
), by=.(TargetGene, cell_type)]

fwrite(gene_summary, file.path(opt$output_dir, "gene_summary.tsv"), sep="\t")
message("  ", nrow(gene_summary), " (gene x cell type) rows")
message("  Wrote gene_summary.tsv")

## ============================================================================
## Step 5: Pairwise metrics
## ============================================================================

n_pairs <- choose(n_ct, 2)
message("\nStep 5: Computing pairwise metrics (", n_pairs, " pairs)...")

# Pre-split by cell type for fast repeated access
sm_by_ct <- split(score_matrix, by="cell_type", keep.by=TRUE)

pairs <- combn(cell_types, 2, simplify=FALSE)

pairwise_results <- rbindlist(lapply(pairs, function(pair) {
    ct_A <- pair[1]; ct_B <- pair[2]
    sm_A <- sm_by_ct[[ct_A]]; sm_B <- sm_by_ct[[ct_B]]
    if (is.null(sm_A) || is.null(sm_B)) return(NULL)

    # Full outer join: all non-self-promoter (consensus_id, TargetGene) pairs in either cell type
    merged <- merge(
        sm_A[isSelfPromoter == FALSE, .(consensus_id, TargetGene, scoreA = score)],
        sm_B[isSelfPromoter == FALSE, .(consensus_id, TargetGene, scoreB = score)],
        by=c("consensus_id", "TargetGene"), all=TRUE
    )
    merged[is.na(scoreA), scoreA := 0]
    merged[is.na(scoreB), scoreB := 0]

    # Jaccard on above-threshold links
    above_A <- merged$scoreA >= opt$threshold
    above_B <- merged$scoreB >= opt$threshold
    nTotalA <- sum(above_A)
    nTotalB <- sum(above_B)
    nShared <- sum(above_A & above_B)
    denom   <- nTotalA + nTotalB - nShared
    jaccard <- if (denom > 0L) nShared / denom else NA_real_

    # Score correlations across all pairs (absent = 0)
    pearson       <- cor(merged$scoreA, merged$scoreB, method="pearson")
    pearson_log1p <- cor(log1p(merged$scoreA), log1p(merged$scoreB), method="pearson")
    spearman      <- cor(merged$scoreA, merged$scoreB, method="spearman")

    data.table(
        biosampleA     = ct_A,
        biosampleB     = ct_B,
        nSharedPredAwB = nShared,
        nTotalPredA    = nTotalA,
        nTotalPredB    = nTotalB,
        jaccard        = jaccard,
        Pearson        = pearson,
        Pearson_log1p  = pearson_log1p,
        Spearman       = spearman
    )
}))

fwrite(pairwise_results, file.path(opt$output_dir, "pairwise_metrics.tsv"), sep="\t")
message("  Wrote pairwise_metrics.tsv")

## ============================================================================
## Summary
## ============================================================================

message("\n=== precompute_comparisons.R complete ===")
message("Output directory: ", opt$output_dir)
message("  consensus_elements.tsv  - ", n_consensus, " elements")
message("  consensus_scores.tsv.gz - ", nrow(score_matrix), " rows")
message("  link_specificity.tsv    - ", nrow(link_spec), " links above threshold")
message("  gene_summary.tsv        - ", nrow(gene_summary), " (gene x cell type) rows")
message("  pairwise_metrics.tsv    - ", nrow(pairwise_results), " pairs")
message("\nNext steps:")
message("  - Re-run with --merge_gap 500 (or larger) to build a neighborhood-level catalog")
message("  - Use pairwise_metrics.tsv as input to plot_global_comparison.R")
message("  - Use consensus_scores.tsv.gz + link_specificity.tsv as input to gene/locus plots")
