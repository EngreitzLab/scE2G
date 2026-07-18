library(data.table)
library(ggplot2)
library(jsonlite)
library(rmarkdown)
library(scales)
library(DT)

## Inputs from snakemake (now lists across all clusters)
annotated_paths      <- snakemake@input$annotated
summary_paths        <- snakemake@input$summary
gene_coverage_paths  <- snakemake@input$gene_coverage
params_json_paths    <- snakemake@input$params_json
subsample_stats_paths <- snakemake@input$subsample_stats
full_run_qc_stats_path <- snakemake@input$full_run_qc_stats

clusters        <- snakemake@params$clusters
model           <- snakemake@params$model
mode            <- snakemake@params$mode
model_thresholds <- snakemake@params$model_thresholds
n_subsamples    <- snakemake@params$n_subsamples
report_out      <- snakemake@output$report

## Nature color palette
COL_BLUE       <- "#006eae"
COL_BLUE_LIGHT <- "#9bcae9"
COL_RED        <- "#c5373d"
COL_YELLOW     <- "#ca9b23"
COL_GREEN      <- "#429130"
COL_GREY       <- "#6e788d"

nature_theme <- theme_classic() +
    theme(
        axis.text  = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line  = element_line(color = "black"),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        strip.text = element_text(color = "black"),
        panel.grid = element_blank()
    )

## Create color palette for clusters
n_clusters <- length(clusters)
if (n_clusters <= 8) {
    cluster_colors <- c(COL_BLUE, COL_RED, COL_GREEN, COL_YELLOW,
                        COL_GREY, "#a64691", "#49bcbc", "#f29742")[1:n_clusters]
} else {
    cluster_colors <- scales::hue_pal()(n_clusters)
}
names(cluster_colors) <- clusters

## Determine if this is a scATAC model (no UMIs)
is_scatac <- grepl("scATAC", model, ignore.case = TRUE)

## ============================================================================
## Load lightweight data (small files only)
## ============================================================================

message("Loading lightweight data")

## Load summary stats for all clusters (small)
all_summary <- rbindlist(lapply(summary_paths, fread))

## Load params JSON for all clusters (small)
all_params <- lapply(params_json_paths, fromJSON)
names(all_params) <- clusters

## Load subsample stats for all clusters (small)
all_subsample_stats <- rbindlist(lapply(seq_along(clusters), function(i) {
    dt <- fread(subsample_stats_paths[[i]])
    dt[, original_cluster := clusters[[i]]]
    dt
}))

## Load full run QC stats (small)
full_run_qc <- fread(full_run_qc_stats_path)
full_run_qc <- full_run_qc[cluster %in% clusters]

## Load gene coverage stats for all clusters (small-medium)
all_gene_coverage <- rbindlist(lapply(gene_coverage_paths, fread))

## ============================================================================
## Section 1: Basic parameters
## ============================================================================

fractions_by_cluster <- sapply(clusters, function(cl) all_params[[cl]]$actual_fraction)
if (length(unique(fractions_by_cluster)) == 1) {
    fraction_str <- sprintf("%.1f%%", fractions_by_cluster[[1]] * 100)
} else {
    fraction_str <- paste(
        sprintf("%s: %.1f%%", names(fractions_by_cluster), fractions_by_cluster * 100),
        collapse = "; "
    )
}

mode_desc <- ifelse(mode == "full", "re-run from fragments", "reuse full-run peaks, re-run from neighborhoods")

params_table <- data.table(
    Parameter = c("Cell types:", "Model name:", "Number of subsamples:", "Fraction of cells sampled:", "Run mode:"),
    Value = c(
        paste(clusters, collapse = ", "),
        model,
        as.character(n_subsamples),
        fraction_str,
        sprintf("%s (%s)", mode, mode_desc)
    )
)

## ============================================================================
## Section 2: Subsample sizes bar plots
## ============================================================================

depth_data <- rbindlist(lapply(clusters, function(cl) {
    full_row <- full_run_qc[cluster == cl]
    sub_stats <- all_subsample_stats[original_cluster == cl]

    data.table(
        cluster = cl,
        cells_full = full_row$cell_count[1],
        cells_subsample = mean(sub_stats$cell_count, na.rm = TRUE),
        fragments_full = full_row$fragments_total[1],
        fragments_subsample = mean(sub_stats$fragments_total, na.rm = TRUE),
        umis_full = if ("umi_count" %in% names(full_row)) full_row$umi_count[1] else NA_real_,
        umis_subsample = if ("umi_count" %in% names(sub_stats)) mean(sub_stats$umi_count, na.rm = TRUE) else NA_real_
    )
}))

depth_long <- melt(depth_data, id.vars = "cluster",
                   measure.vars = patterns("^cells_|^fragments_|^umis_"),
                   variable.name = "metric_type", value.name = "count")
depth_long[, metric := fcase(
    grepl("cells", metric_type), "Cells",
    grepl("fragments", metric_type), "Fragments",
    grepl("umis", metric_type), "UMIs"
)]
depth_long[, source := fifelse(grepl("_full$", metric_type), "Full cluster", "Subsample (avg)")]
depth_long[, metric := factor(metric, levels = c("Cells", "Fragments", "UMIs"))]
depth_long[, source := factor(source, levels = c("Full cluster", "Subsample (avg)"))]

make_depth_plot <- function(metric_name, threshold, y_label) {
    plot_data <- depth_long[metric == metric_name & !is.na(count)]
    if (nrow(plot_data) == 0) return(NULL)

    ggplot(plot_data, aes(x = cluster, y = count, fill = source)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
        geom_hline(yintercept = threshold, color = COL_RED, linetype = "dashed", linewidth = 0.7) +
        scale_fill_manual(values = c("Full cluster" = COL_BLUE, "Subsample (avg)" = COL_BLUE_LIGHT), name = NULL) +
        scale_y_continuous(labels = label_comma()) +
        labs(x = NULL, y = y_label, title = metric_name) +
        nature_theme +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
}

p_cells <- make_depth_plot("Cells", 100, "Number of cells")
p_fragments <- make_depth_plot("Fragments", 1e6, "Number of fragments")
p_umis <- if (!is_scatac) make_depth_plot("UMIs", 2e6, "Number of UMIs") else NULL

## ============================================================================
## MEMORY-EFFICIENT: Process each cluster one at a time
## Extract only the summary data needed for plots, then discard raw data
## ============================================================================

message("Processing clusters one at a time for memory efficiency")

## Pre-allocate lists for summary data (NOT raw data or full plots)
hexbin_data <- list()
ribbon_data_list <- list()
prediction_stability_data <- list()
element_coverage_data <- list()
cluster_metadata <- list()

## Helper: create adaptive bins with higher density around threshold
make_adaptive_breaks <- function(scores, threshold, n_threshold_bins = 50, n_other_bins = 100) {
    min_score <- min(scores, na.rm = TRUE)
    max_score <- max(scores, na.rm = TRUE)
    zone_half_width <- 0.05
    zone_low <- max(min_score, threshold - zone_half_width)
    zone_high <- min(max_score, threshold + zone_half_width)

    breaks_below <- if (min_score < zone_low) {
        n_below <- round(n_other_bins * (zone_low - min_score) / ((max_score - min_score) - (zone_high - zone_low)))
        seq(min_score, zone_low, length.out = max(2, n_below + 1))
    } else numeric(0)

    breaks_zone <- seq(zone_low, zone_high, length.out = n_threshold_bins + 1)

    breaks_above <- if (max_score > zone_high) {
        n_above <- round(n_other_bins * (max_score - zone_high) / ((max_score - min_score) - (zone_high - zone_low)))
        seq(zone_high, max_score, length.out = max(2, n_above + 1))
    } else numeric(0)

    sort(unique(c(breaks_below, breaks_zone, breaks_above)))
}

for (i in seq_along(clusters)) {
    cl <- clusters[[i]]
    threshold <- model_thresholds[[cl]]
    message(sprintf("Processing cluster %d/%d: %s", i, n_clusters, cl))

    ## Load only required columns to save memory
    pred <- fread(annotated_paths[[i]],
                  select = c("chr", "start", "end", "TargetGene",
                            "E2G.Score.qnorm", "E2G.Score.subsample_mean",
                            "E2G.Score.CI_high", "E2G.Score.CI_low",
                            "E2G.Score.pct_above_threshold", "n_subsamples_with_match"))

    n_pairs <- nrow(pred)
    max_score <- max(c(pred$E2G.Score.qnorm, pred$E2G.Score.subsample_mean), na.rm = TRUE)

    ## Store metadata
    cluster_metadata[[cl]] <- list(
        n_pairs = n_pairs,
        n_above = sum(pred$E2G.Score.qnorm >= threshold, na.rm = TRUE),
        threshold = threshold,
        max_score = max_score
    )

    ## --- Hexbin: pre-compute 2D histogram counts ---
    x_breaks <- seq(0, max_score, length.out = 81)
    y_breaks <- seq(0, max_score, length.out = 81)
    pred[, x_bin := findInterval(E2G.Score.qnorm, x_breaks, all.inside = TRUE)]
    pred[, y_bin := findInterval(E2G.Score.subsample_mean, y_breaks, all.inside = TRUE)]
    hex_counts <- pred[, .N, by = .(x_bin, y_bin)]
    hex_counts[, x_mid := x_breaks[x_bin] + diff(x_breaks)[1]/2]
    hex_counts[, y_mid := y_breaks[y_bin] + diff(y_breaks)[1]/2]
    hex_counts[, cluster := cl]
    hexbin_data[[cl]] <- hex_counts[, .(cluster, x_mid, y_mid, N)]

    ## --- Ribbon: pre-compute quantile summaries per bin ---
    breaks <- make_adaptive_breaks(pred$E2G.Score.qnorm, threshold)
    pred[, score_bin := cut(E2G.Score.qnorm, breaks = breaks, include.lowest = TRUE)]

    ribbon_summary <- pred[!is.na(score_bin), .(
        median = median(E2G.Score.subsample_mean, na.rm = TRUE),
        q025 = quantile(E2G.Score.subsample_mean, 0.025, na.rm = TRUE),
        q10 = quantile(E2G.Score.subsample_mean, 0.10, na.rm = TRUE),
        q25 = quantile(E2G.Score.subsample_mean, 0.25, na.rm = TRUE),
        q75 = quantile(E2G.Score.subsample_mean, 0.75, na.rm = TRUE),
        q90 = quantile(E2G.Score.subsample_mean, 0.90, na.rm = TRUE),
        q975 = quantile(E2G.Score.subsample_mean, 0.975, na.rm = TRUE),
        n = .N
    ), by = score_bin]
    ## Parse bin boundaries from factor labels like "[0.1,0.2)" or "(0.1,0.2]"
    ribbon_summary[, score_bin_chr := as.character(score_bin)]
    ## Remove brackets using alternation pattern
    ribbon_summary[, score_bin_clean := gsub("\\[|\\]|\\(|\\)", "", score_bin_chr)]
    ribbon_summary[, c("bin_low_str", "bin_high_str") := tstrsplit(score_bin_clean, ",", fixed = TRUE)]
    ribbon_summary[, bin_low := as.numeric(bin_low_str)]
    ribbon_summary[, bin_high := as.numeric(bin_high_str)]
    ribbon_summary[, bin_mid := (bin_low + bin_high) / 2]
    message(sprintf("  Ribbon data: %d bins, example bin label: '%s' -> clean: '%s' -> low=%s, high=%s",
                    nrow(ribbon_summary),
                    ribbon_summary$score_bin_chr[1],
                    ribbon_summary$score_bin_clean[1],
                    ribbon_summary$bin_low_str[1],
                    ribbon_summary$bin_high_str[1]))
    message(sprintf("  Ribbon data: bin_mid range [%.4f, %.4f], %d NA bin_mid values",
                    min(ribbon_summary$bin_mid, na.rm = TRUE),
                    max(ribbon_summary$bin_mid, na.rm = TRUE),
                    sum(is.na(ribbon_summary$bin_mid))))
    ribbon_summary[, c("score_bin_chr", "score_bin_clean", "bin_low_str", "bin_high_str") := NULL]
    ribbon_summary[, cluster := cl]
    ribbon_data_list[[cl]] <- ribbon_summary[order(bin_mid)]

    ## --- Prediction stability: bin and summarize ---
    pred_breaks <- seq(0, max_score + 0.001, length.out = 51)
    pred[, pred_bin := cut(E2G.Score.qnorm, breaks = pred_breaks, include.lowest = TRUE)]
    pred_summary <- pred[!is.na(pred_bin), .(
        mean_pct_above = mean(E2G.Score.pct_above_threshold, na.rm = TRUE),
        n_pairs = .N
    ), by = pred_bin]
    ## Parse bin boundaries from factor labels like "[0.1,0.2)" or "(0.1,0.2]"
    pred_summary[, pred_bin_chr := as.character(pred_bin)]
    ## Remove brackets using alternation pattern
    pred_summary[, pred_bin_clean := gsub("\\[|\\]|\\(|\\)", "", pred_bin_chr)]
    pred_summary[, c("bin_low_str", "bin_high_str") := tstrsplit(pred_bin_clean, ",", fixed = TRUE)]
    pred_summary[, bin_low := as.numeric(bin_low_str)]
    pred_summary[, bin_high := as.numeric(bin_high_str)]
    pred_summary[, bin_mid := (bin_low + bin_high) / 2]
    message(sprintf("  Pred stability: %d bins, example: '%s' -> '%s' -> low=%s, high=%s",
                    nrow(pred_summary),
                    pred_summary$pred_bin_chr[1],
                    pred_summary$pred_bin_clean[1],
                    pred_summary$bin_low_str[1],
                    pred_summary$bin_high_str[1]))
    pred_summary[, c("pred_bin_chr", "pred_bin_clean", "bin_low_str", "bin_high_str") := NULL]
    pred_summary[, cluster := cl]
    pred_summary[, threshold := threshold]
    prediction_stability_data[[cl]] <- pred_summary

    ## --- Element coverage: aggregate by element ---
    elem_cov <- pred[, .(element_id = paste0(chr, ":", start, "-", end), n_subsamples_with_match)]
    elem_cov <- elem_cov[, .(n_subsamples_covered = max(n_subsamples_with_match, na.rm = TRUE)), by = element_id]
    elem_hist <- elem_cov[, .N, by = n_subsamples_covered]
    elem_hist[, total := sum(N)]
    elem_hist[, pct := N / total]
    elem_hist[, cluster := cl]
    element_coverage_data[[cl]] <- elem_hist

    ## Clean up and force garbage collection
    rm(pred, hex_counts, ribbon_summary, pred_summary, elem_cov, elem_hist)
    gc(verbose = FALSE)
}

message("Finished processing all clusters, combining summary data")

## Combine summary data
all_hexbin <- rbindlist(hexbin_data)
all_ribbon <- rbindlist(ribbon_data_list)
all_pred_stability <- rbindlist(prediction_stability_data)
all_element_coverage <- rbindlist(element_coverage_data)

rm(hexbin_data, ribbon_data_list, prediction_stability_data, element_coverage_data)
gc(verbose = FALSE)

## ============================================================================
## Create plots from pre-computed summary data
## ============================================================================

message("Creating plots from summary data")

## --- Hexbin plots ---
score_stability_hexbin <- lapply(clusters, function(cl) {
    meta <- cluster_metadata[[cl]]
    dd <- all_hexbin[cluster == cl]

    ggplot(dd, aes(x = x_mid, y = y_mid, fill = N)) +
        geom_tile() +
        geom_abline(slope = 1, intercept = 0, color = COL_RED, linetype = "dashed", linewidth = 0.6) +
        geom_vline(xintercept = meta$threshold, color = COL_GREY, linetype = "dotted", linewidth = 0.5) +
        geom_hline(yintercept = meta$threshold, color = COL_GREY, linetype = "dotted", linewidth = 0.5) +
        scale_fill_viridis_c(trans = "log10", name = "Count", option = "viridis", na.value = "white") +
        coord_fixed(ratio = 1, xlim = c(0, meta$max_score), ylim = c(0, meta$max_score)) +
        labs(x = "Full run score", y = "Subsample mean score", title = cl,
             subtitle = sprintf("n = %s pairs", format(meta$n_pairs, big.mark = ","))) +
        nature_theme + theme(legend.position = "right")
})
names(score_stability_hexbin) <- clusters

## --- Ribbon plots ---
score_stability_ribbon <- lapply(clusters, function(cl) {
    meta <- cluster_metadata[[cl]]
    dd <- all_ribbon[cluster == cl]

    ggplot(dd, aes(x = bin_mid)) +
        geom_ribbon(aes(ymin = q025, ymax = q975), fill = COL_BLUE, alpha = 0.2) +
        geom_ribbon(aes(ymin = q10, ymax = q90), fill = COL_BLUE, alpha = 0.3) +
        geom_ribbon(aes(ymin = q25, ymax = q75), fill = COL_BLUE, alpha = 0.4) +
        geom_line(aes(y = median), color = COL_BLUE, linewidth = 0.8) +
        geom_abline(slope = 1, intercept = 0, color = COL_RED, linetype = "dashed", linewidth = 0.6) +
        geom_vline(xintercept = meta$threshold, color = COL_GREY, linetype = "dotted", linewidth = 0.5) +
        geom_hline(yintercept = meta$threshold, color = COL_GREY, linetype = "dotted", linewidth = 0.5) +
        coord_fixed(ratio = 1, xlim = c(0, meta$max_score), ylim = c(0, meta$max_score)) +
        labs(x = "Full run score (binned)", y = "Subsample mean score", title = cl,
             subtitle = "Ribbons: 50% (dark), 80%, 95% (light) intervals") +
        nature_theme
})
names(score_stability_ribbon) <- clusters

## --- Prediction stability plots ---
common_threshold <- median(unlist(model_thresholds))
zoom_range <- c(common_threshold - 0.1, common_threshold + 0.1)

p_pred_stability_full <- ggplot(all_pred_stability,
                                 aes(x = bin_mid, y = mean_pct_above, color = cluster)) +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    geom_point(aes(size = log10(n_pairs + 1)), alpha = 0.6) +
    geom_vline(xintercept = common_threshold, color = COL_GREY, linetype = "dashed", linewidth = 0.6) +
    geom_hline(yintercept = 0.5, color = COL_GREY, linetype = "dotted", linewidth = 0.4) +
    scale_color_manual(values = cluster_colors, name = "Cell type") +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    scale_size_continuous(name = "log10(n pairs)", range = c(0.5, 3), guide = "none") +
    labs(x = "Full run score", y = "Fraction of subsamples above threshold",
         title = "Prediction stability across subsamples", subtitle = "Full score range") +
    nature_theme + theme(legend.position = "right")

p_pred_stability_zoom <- ggplot(all_pred_stability,
                                 aes(x = bin_mid, y = mean_pct_above, color = cluster)) +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    geom_point(aes(size = log10(n_pairs + 1)), alpha = 0.6) +
    geom_vline(xintercept = common_threshold, color = COL_GREY, linetype = "dashed", linewidth = 0.6) +
    geom_hline(yintercept = 0.5, color = COL_GREY, linetype = "dotted", linewidth = 0.4) +
    scale_color_manual(values = cluster_colors, name = "Cell type") +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    scale_size_continuous(name = "log10(n pairs)", range = c(0.5, 3), guide = "none") +
    coord_cartesian(xlim = zoom_range) +
    labs(x = "Full run score", y = "Fraction of subsamples above threshold",
         title = "Threshold region detail",
         subtitle = sprintf("Score range: %.2f to %.2f", zoom_range[1], zoom_range[2])) +
    nature_theme + theme(legend.position = "none")

## --- Element coverage plot ---
p_element_coverage <- ggplot(all_element_coverage,
                              aes(x = factor(n_subsamples_covered), y = pct, fill = cluster)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
    scale_fill_manual(values = cluster_colors, name = "Cell type") +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = "Number of subsamples with element", y = "Percent of unique elements",
         title = "Element coverage across subsamples",
         subtitle = "0 = element not found in any subsample (peak absent)") +
    nature_theme + theme(legend.position = "bottom")

## --- Gene coverage plots ---
gene_coverage_full <- all_gene_coverage[in_full_run == TRUE]
gene_coverage_full_hist <- gene_coverage_full[, .(n_genes = .N), by = .(cluster, n_subsamples_with_gene)]
gene_coverage_full_hist[, total_genes := sum(n_genes), by = cluster]
gene_coverage_full_hist[, pct_genes := n_genes / total_genes]

p_gene_coverage_full <- if (!is_scatac) {
    ggplot(gene_coverage_full_hist,
           aes(x = factor(n_subsamples_with_gene), y = pct_genes, fill = cluster)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
        scale_fill_manual(values = cluster_colors, name = "Cell type") +
        scale_y_continuous(labels = scales::percent_format()) +
        labs(x = "Number of subsamples with gene above threshold",
             y = "Percent of genes",
             title = "Genes in full run: subsample coverage",
             subtitle = "For genes above threshold in full run, how many subsamples also have them?") +
        nature_theme + theme(legend.position = "bottom")
} else NULL

gene_coverage_sub_only <- all_gene_coverage[in_full_run == FALSE]
n_genes_sub_only <- nrow(gene_coverage_sub_only)

p_gene_coverage_sub_only <- if (!is_scatac && n_genes_sub_only > 0) {
    gene_coverage_sub_hist <- gene_coverage_sub_only[, .(n_genes = .N), by = .(cluster, n_subsamples_with_gene)]
    gene_coverage_sub_hist[, total_genes := sum(n_genes), by = cluster]
    gene_coverage_sub_hist[, pct_genes := n_genes / total_genes]

    ggplot(gene_coverage_sub_hist,
           aes(x = factor(n_subsamples_with_gene), y = pct_genes, fill = cluster)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
        scale_fill_manual(values = cluster_colors, name = "Cell type") +
        scale_y_continuous(labels = scales::percent_format()) +
        labs(x = "Number of subsamples with gene above threshold",
             y = "Percent of genes",
             title = "Genes only in subsamples: coverage distribution",
             subtitle = sprintf("Genes above threshold in subsamples but not in full run (n = %d)", n_genes_sub_only)) +
        nature_theme + theme(legend.position = "bottom")
} else NULL

## ============================================================================
## Build RMarkdown report
## ============================================================================

message("Building report")
dir.create(dirname(report_out), recursive = TRUE, showWarnings = FALSE)

plots_rds <- tempfile(fileext = ".rds")
saveRDS(list(
    params_table = params_table,
    p_cells = p_cells,
    p_fragments = p_fragments,
    p_umis = p_umis,
    score_stability_hexbin = score_stability_hexbin,
    score_stability_ribbon = score_stability_ribbon,
    p_pred_stability_full = p_pred_stability_full,
    p_pred_stability_zoom = p_pred_stability_zoom,
    p_element_coverage = p_element_coverage,
    p_gene_coverage_full = p_gene_coverage_full,
    p_gene_coverage_sub_only = p_gene_coverage_sub_only,
    n_genes_sub_only = n_genes_sub_only,
    clusters = clusters,
    is_scatac = is_scatac,
    depth_data = depth_data
), plots_rds)

## Build score stability tabs (one tab per cluster with both hexbin and ribbon)
score_stability_tabs <- ""
for (i in seq_along(clusters)) {
    cl <- clusters[[i]]
    cl_safe <- gsub("[^a-zA-Z0-9]", "_", cl)
    score_stability_tabs <- paste0(score_stability_tabs, sprintf('
### %s

```{r stability_%s, fig.width=12, fig.height=5}
score_stability_hexbin[["%s"]] + score_stability_ribbon[["%s"]] + plot_layout(ncol = 2)
```

', cl, cl_safe, cl, cl))
}

depth_plots_chunk <- '
```{r depth_plots, fig.width=12, fig.height=5}
library(patchwork)
if (is_scatac) {
    p_cells + p_fragments + plot_layout(ncol = 2, guides = "collect") & theme(legend.position = "bottom")
} else {
    p_cells + p_fragments + p_umis + plot_layout(ncol = 3, guides = "collect") & theme(legend.position = "bottom")
}
```
'

rmd_content <- paste0('---
title: "scE2G CI QC report"
output:
  html_document:
    self_contained: true
    theme: flatly
    toc: true
    toc_float: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo=FALSE, warning=FALSE, message=FALSE, fig.width=10, fig.height=5)
library(knitr)
library(patchwork)
library(DT)
objs <- readRDS("', plots_rds, '")
list2env(objs, envir=environment())
```
## Run parameters

```{r params_table}
DT::datatable(params_table, colnames = c("", ""), rownames = FALSE,
              options = list(dom = "t", ordering = FALSE))
```

## Subsample sizes per cell type

Darker bars show the full cluster depth; lighter bars show the average subsample depth.
Dashed red lines indicate minimum thresholds (cells: 100, fragments: 1M', ifelse(!is_scatac, ', UMIs: 2M', ''), ').

', depth_plots_chunk, '

### Subsample summary table

```{r depth_table}
display_cols <- c("cluster", "cells_full", "cells_subsample", "fragments_full", "fragments_subsample"',
ifelse(!is_scatac, ', "umis_full", "umis_subsample"', ''), ')
col_names <- c("Cell type", "Cells (full)", "Cells (subsample)", "Fragments (full)", "Fragments (subsample)"',
ifelse(!is_scatac, ', "UMIs (full)", "UMIs (subsample)"', ''), ')
display_df <- depth_data[, ..display_cols]
setnames(display_df, col_names)
DT::datatable(display_df, rownames = FALSE,
              options = list(scrollX = TRUE, pageLength = 25)) %>%
  DT::formatRound(columns = 2:ncol(display_df), digits = 0, mark = ",")
```

## Score stability {.tabset}

Each tab shows two plots for a cell type:

- **Left (2D histogram):** Full-run score (x) vs subsample mean score (y). Color = pair density (log scale). Red dashed = identity; dotted = threshold.
- **Right (quantile ribbons):** Ribbons show 50% (dark), 80%, 95% (light) intervals. Solid line = median. Red dashed = identity.

', score_stability_tabs, '

## {-}

## Prediction stability

For each E-G pair, what fraction of subsamples also call it above threshold? This shows how consistently pairs are classified across subsamples. Pairs near the threshold (vertical dashed line) are expected to have ~50% stability; pairs well above should approach 100%.

```{r pred_stability, fig.width=14, fig.height=6}
p_pred_stability_full + p_pred_stability_zoom + plot_layout(ncol = 2, widths = c(1.2, 1), guides = "collect") & theme(legend.position = "bottom")
```

## Element and gene coverage {.tabset}

### Element coverage

For each unique element (enhancer region), how many subsamples contain an overlapping element? Elements with 0 coverage were not called as peaks in any subsample.

```{r element_coverage, fig.width=12, fig.height=5}
p_element_coverage
```

', ifelse(!is_scatac, '
### Gene coverage: full run

For genes above threshold in the full run, how many subsamples also have them above threshold?

```{r gene_coverage_full, fig.width=12, fig.height=5}
p_gene_coverage_full
```

### Gene coverage: subsample-only

Genes that appear above threshold in at least one subsample but not in the full run. These represent genes that may be borderline expressed or have variable enhancer activity.

```{r gene_coverage_sub_only, fig.width=12, fig.height=5}
if (!is.null(p_gene_coverage_sub_only) && n_genes_sub_only > 0) {
    p_gene_coverage_sub_only
} else {
    cat("No genes found that are above threshold in subsamples but not in the full run.")
}
```

', ''), '
## {-}


')

rmd_tmp <- tempfile(fileext = ".Rmd")
writeLines(rmd_content, rmd_tmp)

render(rmd_tmp, output_file = report_out, quiet = TRUE)
unlink(c(rmd_tmp, plots_rds))

message("CI QC report written to: ", report_out)
