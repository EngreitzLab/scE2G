library(data.table)
library(jsonlite)

## Inputs from snakemake
full_predictions_path <- snakemake@input$full_predictions
subsample_prediction_paths <- snakemake@input$subsample_predictions
params_json_path <- snakemake@input$params_json

model_threshold <- snakemake@params$model_threshold
ci_level <- snakemake@params$ci_level
overlap_fraction <- snakemake@params$overlap_fraction
cluster <- snakemake@params$cluster
model <- snakemake@params$model
mode <- snakemake@params$mode
n_subsamples <- snakemake@params$n_subsamples

annotated_out <- snakemake@output$annotated
summary_out <- snakemake@output$summary
gene_coverage_out <- snakemake@output$gene_coverage

## Load full predictions
message("Loading full predictions: ", full_predictions_path)
full_pred <- fread(full_predictions_path)
n_pairs <- nrow(full_pred)
message("  ", n_pairs, " E-G pairs")

## Add row index for safe re-joining after bedtools
full_pred[, .row_id := .I]

## Build element BED from full predictions (unique elements only)
full_elements <- unique(full_pred[, .(chr, start, end)])
full_elements[, element_id := paste0(chr, ":", start, "-", end)]

## Build lookup from full predictions: element_id + TargetGene -> row_id
full_pred[, element_id := paste0(chr, ":", start, "-", end)]

## Initialize score matrix (rows = E-G pairs, cols = subsamples)
score_matrix <- matrix(0.0, nrow = n_pairs, ncol = n_subsamples)

## Track genes above threshold in each subsample (for gene coverage analysis)
subsample_genes_above_threshold <- vector("list", n_subsamples)

## Helper: find matching subsample scores for all full-run pairs via bedtools intersect
## Also returns genes above threshold in the subsample
score_one_subsample <- function(sub_pred_path, subsample_idx) {
    message("Processing subsample ", subsample_idx, ": ", sub_pred_path)

    if (!file.exists(sub_pred_path)) {
        warning("Subsample predictions not found: ", sub_pred_path, " - treating as 0 for all pairs")
        return(list(scores = rep(0.0, n_pairs), genes_above_threshold = character(0)))
    }

    sub_pred <- fread(sub_pred_path)
    if (nrow(sub_pred) == 0) {
        return(list(scores = rep(0.0, n_pairs), genes_above_threshold = character(0)))
    }

    ## Get genes above threshold in this subsample
    genes_above_thresh <- unique(sub_pred[E2G.Score.qnorm >= model_threshold, TargetGene])

    sub_pred[, element_id := paste0(chr, ":", start, "-", end)]

    sub_elements <- unique(sub_pred[, .(chr, start, end)])

    ## Write BED files for bedtools intersect
    full_bed <- tempfile(fileext = ".bed")
    sub_bed <- tempfile(fileext = ".bed")
    intersect_out <- tempfile(fileext = ".bed")

    fwrite(full_elements[, .(chr, start, end, element_id)], full_bed, sep = "\t", col.names = FALSE)
    fwrite(sub_elements[, .(chr, start, end)], sub_bed, sep = "\t", col.names = FALSE)

    ## -f: min fraction of A covered by B; -r: reciprocal requirement; -wao: write all A with overlap
    cmd <- paste(
        "bedtools intersect",
        "-a", full_bed,
        "-b", sub_bed,
        "-f", overlap_fraction, "-r",
        "-wao >", intersect_out
    )
    system(cmd, ignore.stdout = FALSE, ignore.stderr = FALSE)

    ## Read intersect results
    ## Columns: chr_f, start_f, end_f, element_id_f, chr_s, start_s, end_s, overlap_len
    hits <- fread(intersect_out, header = FALSE, col.names = c(
        "chr_f", "start_f", "end_f", "element_id_f",
        "chr_s", "start_s", "end_s", "overlap_len"
    ))

    ## Keep only real overlaps (overlap_len > 0; bedtools uses -1 for no match in -wao)
    hits <- hits[overlap_len > 0]

    ## For each full element, keep the best-overlapping subsample element
    hits <- hits[, .SD[which.max(overlap_len)], by = element_id_f]

    ## Build subsample element_id for lookup
    hits[, sub_element_id := paste0(chr_s, ":", start_s, "-", end_s)]

    ## Build subsample score lookup: (sub_element_id, TargetGene) -> E2G.Score.qnorm
    sub_scores <- sub_pred[, .(sub_element_id = element_id, TargetGene, score = E2G.Score.qnorm)]

    ## Map full predictions -> matched sub element -> score
    ## Join: full_pred + element_map + sub_scores
    tmp <- copy(full_pred[, .(.row_id, element_id_f = element_id, TargetGene)])
    tmp <- merge(tmp, hits[, .(element_id_f, sub_element_id)], by = "element_id_f", all.x = TRUE)
    tmp <- merge(tmp, sub_scores, by.x = c("sub_element_id", "TargetGene"),
                 by.y = c("sub_element_id", "TargetGene"), all.x = TRUE)

    ## Sort back to original row order
    setorder(tmp, .row_id)
    scores <- tmp$score
    scores[is.na(scores)] <- 0.0

    ## Cleanup temp files
    unlink(c(full_bed, sub_bed, intersect_out))

    list(scores = scores, genes_above_threshold = genes_above_thresh)
}

## Process each subsample
for (i in seq_along(subsample_prediction_paths)) {
    result <- score_one_subsample(subsample_prediction_paths[[i]], i - 1)
    score_matrix[, i] <- result$scores
    subsample_genes_above_threshold[[i]] <- result$genes_above_threshold
}

## Compute per-pair CI statistics
ci_alpha <- (1 - ci_level) / 2

## Warn if n_subsamples is too small for the requested ci_level.
## Empirical quantiles at ci_alpha require at least floor(ci_alpha * n) >= 1,
## i.e. n >= ceil(1 / ci_alpha). Below this, CI_low == min and CI_high == max.
min_n_for_ci <- ceiling(1 / ci_alpha)
if (n_subsamples < min_n_for_ci) {
    warning(sprintf(
        "n_subsamples=%d is too small for ci_level=%.2f (requires >= %d subsamples). CI_low and CI_high will equal the subsample min/max, not a true %.0f%% CI.",
        n_subsamples, ci_level, min_n_for_ci, ci_level * 100
    ))
}

message("Computing CI statistics across ", n_subsamples, " subsamples")
row_means <- rowMeans(score_matrix)
full_pred[, E2G.Score.subsample_mean   := row_means]
full_pred[, E2G.Score.subsample_sd     := sqrt(rowSums((score_matrix - row_means)^2) / max(1L, n_subsamples - 1L))]
full_pred[, E2G.Score.CI_low           := apply(score_matrix, 1, quantile, probs = ci_alpha)]
full_pred[, E2G.Score.CI_high          := apply(score_matrix, 1, quantile, probs = 1 - ci_alpha)]
full_pred[, E2G.Score.pct_above_threshold := rowMeans(score_matrix >= model_threshold)]
full_pred[, n_subsamples_with_match    := rowSums(score_matrix > 0)]
full_pred[, E2G.Score.subsample_scores := apply(score_matrix, 1, function(x) paste(round(x, 6), collapse = ","))]

## Remove temporary helper columns before writing
full_pred[, c(".row_id", "element_id") := NULL]

## Write annotated predictions
message("Writing annotated predictions: ", annotated_out)
dir.create(dirname(annotated_out), recursive = TRUE, showWarnings = FALSE)
fwrite(full_pred, annotated_out, sep = "\t", compress = "gzip")

## Compute summary stats
params_json <- fromJSON(params_json_path)

ci_widths <- full_pred$E2G.Score.CI_high - full_pred$E2G.Score.CI_low
n_above_threshold <- sum(full_pred$E2G.Score.qnorm >= model_threshold, na.rm = TRUE)

summary_dt <- data.table(
    cluster                     = cluster,
    model                       = model,
    mode                        = mode,
    n_subsamples                = as.integer(n_subsamples),
    actual_fraction             = params_json$actual_fraction,
    n_cells                     = params_json$n_cells,
    n_fragments                 = params_json$n_fragments,
    n_cells_per_subsample       = params_json$n_cells_per_subsample,
    n_pairs_full                = n_pairs,
    n_pairs_above_threshold     = n_above_threshold,
    median_CI_width             = median(ci_widths, na.rm = TRUE),
    pct_pairs_CI_width_lt_0.1   = mean(ci_widths < 0.1, na.rm = TRUE),
    median_pct_above_threshold  = median(full_pred$E2G.Score.pct_above_threshold, na.rm = TRUE),
    median_n_subsamples_matched = median(full_pred$n_subsamples_with_match, na.rm = TRUE)
)

message("Writing summary stats: ", summary_out)
fwrite(summary_dt, summary_out, sep = "\t")

## ============================================================================
## Gene coverage analysis (bidirectional)
## ============================================================================

message("Computing gene coverage statistics")

## Genes above threshold in full run
full_genes_above_threshold <- unique(full_pred[E2G.Score.qnorm >= model_threshold, TargetGene])

## All genes in full run (regardless of threshold)
full_genes_all <- unique(full_pred$TargetGene)

## All genes above threshold in any subsample
all_subsample_genes <- unique(unlist(subsample_genes_above_threshold))

## Direction 1: For each gene in full run (above threshold), how many subsamples have it?
gene_coverage_full_to_sub <- data.table(
    TargetGene = full_genes_above_threshold,
    in_full_run = TRUE,
    n_subsamples_with_gene = sapply(full_genes_above_threshold, function(g) {
        sum(sapply(subsample_genes_above_threshold, function(sg) g %in% sg))
    })
)

## Direction 2: Genes in subsamples but not in full run (above threshold)
genes_only_in_subsamples <- setdiff(all_subsample_genes, full_genes_above_threshold)
if (length(genes_only_in_subsamples) > 0) {
    gene_coverage_sub_only <- data.table(
        TargetGene = genes_only_in_subsamples,
        in_full_run = FALSE,
        n_subsamples_with_gene = sapply(genes_only_in_subsamples, function(g) {
            sum(sapply(subsample_genes_above_threshold, function(sg) g %in% sg))
        })
    )
    gene_coverage_dt <- rbind(gene_coverage_full_to_sub, gene_coverage_sub_only)
} else {
    gene_coverage_dt <- gene_coverage_full_to_sub
}

## Add cluster info
gene_coverage_dt[, cluster := cluster]

message("Writing gene coverage stats: ", gene_coverage_out)
fwrite(gene_coverage_dt, gene_coverage_out, sep = "\t")

message("Done.")
