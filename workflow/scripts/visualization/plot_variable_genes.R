#!/usr/bin/env Rscript
#
# plot_variable_genes.R
#
# Diagnostic visualization of top and bottom variable genes from rank_variable_genes.R.
# Panels: n_links CV (bar), CI-weighted score variance (bar),
#         per-cell-type n_links boxplot, log10 TPM boxplot.
# Optionally overlays user-specified custom genes in a separate section.
#
# Filtering: genes with <--min_total_links total links OR with n_links_max < --min_max_links
#   (i.e., no cell type has more than min_max_links enhancers) are excluded.
#
# Usage:
#   Rscript plot_variable_genes.R \
#     --variable_genes results/.../comparison/variable_genes.tsv \
#     --top_n 20 --bottom_n 10 \
#     --custom_genes KDR,GATA1,HAND2,VEGFA \
#     --output results/.../comparison/variable_genes_plot.pdf

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(patchwork)
    library(optparse)
})

option_list <- list(
    make_option("--variable_genes",  type="character", help="Path to variable_genes.tsv"),
    make_option("--top_n",           type="integer",   default=20L,
                help="Most variable genes to show [default: %default]"),
    make_option("--bottom_n",        type="integer",   default=10L,
                help="Least variable genes to show [default: %default]"),
    make_option("--min_total_links", type="integer",   default=5L,
                help="Exclude genes with fewer total links across cell types [default: %default]"),
    make_option("--min_max_links",   type="integer",   default=5L,
                help="Exclude genes where no cell type has >= this many links [default: %default]"),
    make_option("--custom_genes",    type="character", default=NULL,
                help="Comma-separated gene names to plot in an additional section"),
    make_option("--output",          type="character", default=NULL,
                help="Output PDF [default: variable_genes_plot.pdf next to input]")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$variable_genes)) stop("--variable_genes is required")
if (is.null(opt$output))
    opt$output <- sub("\\.tsv$", "_plot.pdf", opt$variable_genes)

comparison_dir    <- dirname(opt$variable_genes)
gene_summary_path <- file.path(comparison_dir, "gene_summary.tsv")

## ============================================================================
## Load and filter
## ============================================================================

dt <- fread(opt$variable_genes)
has_ci <- "score_var_ci_mean" %in% names(dt)
score_col   <- if (has_ci) "score_var_ci_mean" else "score_var_mean"
score_label <- if (has_ci) "CI-weighted\nscore variance" else "Score variance\n(mean)"

dt[, total_links := n_links_mean * n_cell_types]
dt_filt <- dt[total_links >= opt$min_total_links & n_links_max >= opt$min_max_links]
message(nrow(dt_filt), " genes after filtering (", nrow(dt) - nrow(dt_filt), " removed)")

n_top    <- min(opt$top_n,    floor(nrow(dt_filt) / 2))
n_bottom <- min(opt$bottom_n, floor(nrow(dt_filt) / 2))
message("Showing top ", n_top, " and bottom ", n_bottom,
        " genes (of ", nrow(dt_filt), " filtered)")

top_genes    <- dt_filt[seq_len(n_top)]
bottom_genes <- dt_filt[seq(nrow(dt_filt) - n_bottom + 1, nrow(dt_filt))]
top_genes[,    group := "Most variable"]
bottom_genes[, group := "Least variable"]

## ============================================================================
## Custom genes
## ============================================================================

custom_genes_requested <- if (!is.null(opt$custom_genes))
    trimws(strsplit(opt$custom_genes, ",")[[1]]) else character(0)

if (length(custom_genes_requested) > 0) {
    not_found <- setdiff(custom_genes_requested, dt$TargetGene)
    if (length(not_found) > 0)
        message("Custom genes not found: ", paste(not_found, collapse=", "))
    custom_dt <- dt[TargetGene %in% custom_genes_requested]
    custom_dt[, gene_order := match(TargetGene, custom_genes_requested)]
    setorder(custom_dt, gene_order)
    custom_dt[, gene_order := NULL]
    custom_dt[, group := "Custom"]
} else {
    custom_dt <- data.table()
}

has_custom <- nrow(custom_dt) > 0

## ============================================================================
## Combined plot_dt and gene factor
## ============================================================================

plot_dt <- rbind(top_genes, bottom_genes, if (has_custom) custom_dt, fill=TRUE)
group_levels <- c("Most variable", "Least variable", if (has_custom) "Custom")
plot_dt[, group := factor(group, levels=group_levels)]

# Gene factor: most-variable at top of facet (last level = topmost row in ggplot)
gene_levels <- c(
    if (has_custom) rev(custom_dt$TargetGene),
    bottom_genes[order(-combined_rank), TargetGene],
    top_genes[order(-combined_rank),    TargetGene]
)
plot_dt[, gene := factor(TargetGene, levels=gene_levels)]

## ============================================================================
## Load gene_summary.tsv for per-cell-type distributions
## ============================================================================

has_gene_summary <- file.exists(gene_summary_path)
if (has_gene_summary) {
    gs <- fread(gene_summary_path)
    gs_plot <- gs[TargetGene %in% plot_dt$TargetGene]
    gs_plot <- merge(gs_plot, plot_dt[, .(TargetGene, gene, group)],
                     by="TargetGene", all.x=TRUE)
    gs_plot[, log10_TPM := log10(RNA_pseudobulkTPM + 1)]
} else {
    message("gene_summary.tsv not found; skipping distribution panels")
}

## ============================================================================
## Theme
## ============================================================================

group_colors <- c("Most variable"="#c5373d", "Least variable"="#6e788d", "Custom"="#006eae")

base_theme <- theme_classic() +
    theme(
        axis.text        = element_text(color="black", size=8),
        axis.ticks       = element_line(color="black"),
        axis.line        = element_line(color="black"),
        panel.grid       = element_blank(),
        legend.position  = "none",
        strip.background = element_blank(),
        strip.text       = element_blank(),
        panel.spacing    = unit(0.4, "lines")
    )

no_y <- theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank()
)

keep_y_theme <- theme(axis.text.y = element_text(face="italic"))

bar_plot <- function(data, xvar, xlabel, keep_y=FALSE) {
    p <- ggplot(data, aes(y=gene, x=.data[[xvar]], fill=group)) +
        geom_col(width=0.7) +
        scale_fill_manual(values=group_colors) +
        scale_x_continuous(expand=expansion(mult=c(0, 0.08))) +
        facet_grid(group ~ ., scales="free_y", space="free_y") +
        labs(x=xlabel, y=NULL) +
        base_theme +
        if (keep_y) keep_y_theme else no_y
    p
}

box_plot <- function(data, xvar, xlabel, keep_y=FALSE, log10_scale=FALSE) {
    p <- ggplot(data, aes(y=gene, x=.data[[xvar]], color=group)) +
        geom_boxplot(width=0.6, linewidth=0.4, outlier.size=0.8, outlier.stroke=0) +
        scale_color_manual(values=group_colors) +
        facet_grid(group ~ ., scales="free_y", space="free_y") +
        labs(x=xlabel, y=NULL) +
        base_theme + no_y
    if (log10_scale)
        p <- p + scale_x_continuous(
            trans  = "log10",
            labels = function(x) format(x, scientific=FALSE, big.mark=","),
            expand = expansion(mult=c(0.02, 0.08))
        )
    p
}

## ============================================================================
## Build panels
## ============================================================================

p_cv    <- bar_plot(plot_dt, "n_links_cv", "Enhancer count CV\n(across cell types)", keep_y=TRUE)
p_score <- bar_plot(plot_dt, score_col,    score_label,                        keep_y=FALSE)

if (has_gene_summary) {
    p_nlinks <- box_plot(gs_plot, "n_links_above_threshold", "Links above\nthreshold")
    p_tpm    <- box_plot(gs_plot, "log10_TPM", "log10(TPM + 1)", log10_scale=FALSE)
    combined <- p_cv + p_score + p_nlinks + p_tpm +
        plot_layout(nrow=1, widths=c(1.5, 1.5, 1.5, 1.5))
    plot_width <- 10
} else {
    combined <- p_cv + p_score + plot_layout(nrow=1, widths=c(1.5, 1.5))
    plot_width <- 7
}

## ============================================================================
## Save
## ============================================================================

dir.create(dirname(opt$output), recursive=TRUE, showWarnings=FALSE)
n_genes     <- n_top + n_bottom + nrow(custom_dt)
plot_height <- max(3, n_genes * 0.18 + 2)
ggsave(opt$output, combined, width=plot_width, height=plot_height)
message("Wrote: ", opt$output)
