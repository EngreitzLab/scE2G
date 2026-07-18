#!/usr/bin/env Rscript
#
# plot_expression_vs_enhancers.R
#
# Scatter plots relating gene expression changes to enhancer landscape changes
# between two cell types. One point per gene.
#
# Panels produced:
#   1. log2FC(TPM) vs delta max E2G score
#   2. log2FC(TPM) vs delta n_links_above_threshold
#   3. log2FC(TPM) vs delta normalizedATAC_prom (promoter accessibility change)
#
# Genes are colored by expression direction: up / down / stable.
#
# Usage:
#   Rscript plot_expression_vs_enhancers.R \
#     --gene_summary results/.../comparison/gene_summary.tsv \
#     --cell_type_A d0_ipsc \
#     --cell_type_B d4_ec \
#     --output results/.../comparison/expression_vs_enhancers_d0_vs_d4.pdf

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(optparse)
    library(scales)
})

## ============================================================================
## CLI
## ============================================================================

option_list <- list(
    make_option("--gene_summary",     type="character", help="Path to gene_summary.tsv"),
    make_option("--cell_type_A",      type="character", help="Reference cell type"),
    make_option("--cell_type_B",      type="character", help="Comparison cell type"),
    make_option("--output",           type="character", help="Output PDF path"),
    make_option("--min_tpm",          type="numeric",   default=1.0,
                help="Min TPM in at least one cell type (filters unexpressed genes) [default: %default]"),
    make_option("--log2fc_threshold", type="numeric",   default=1.0,
                help="log2FC threshold to call up/down [default: %default]"),
    make_option("--tpm_pseudocount",  type="numeric",   default=0.1,
                help="Pseudocount added before log2FC [default: %default]"),
    make_option("--highlight_genes",  type="character", default=NULL,
                help="Comma-separated gene names to label on plots"),
    make_option("--width",            type="numeric",   default=10),
    make_option("--height",           type="numeric",   default=4)
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$gene_summary)) stop("--gene_summary is required")
if (is.null(opt$cell_type_A))  stop("--cell_type_A is required")
if (is.null(opt$cell_type_B))  stop("--cell_type_B is required")
if (is.null(opt$output))       stop("--output is required")

dir.create(dirname(opt$output), recursive=TRUE, showWarnings=FALSE)

highlight_genes <- if (!is.null(opt$highlight_genes)) {
    trimws(strsplit(opt$highlight_genes, ",")[[1]])
} else character(0)

## ============================================================================
## Theme and colors
## ============================================================================

COL_UP     <- "#c5373d"   # red — up in B
COL_DOWN   <- "#006eae"   # blue — down in B
COL_STABLE <- "#96a0b3"   # grey — stable

nature_theme <- theme_classic() +
    theme(
        axis.text   = element_text(color="black"),
        axis.ticks  = element_line(color="black"),
        axis.line   = element_line(color="black"),
        panel.grid  = element_blank(),
        legend.text = element_text(color="black"),
        legend.title= element_text(color="black"),
        strip.text  = element_text(color="black"),
        strip.background = element_blank()
    )

## ============================================================================
## Load and reshape
## ============================================================================

gs <- fread(opt$gene_summary)

ct_A <- opt$cell_type_A; ct_B <- opt$cell_type_B
for (ct in c(ct_A, ct_B)) {
    if (!ct %in% gs$cell_type) stop("Cell type '", ct, "' not found in gene_summary.tsv")
}

# Pivot to wide: one row per gene with A and B columns
wide_A <- gs[cell_type == ct_A, .(
    TargetGene,
    n_links_A  = n_links_above_threshold,
    max_score_A = max_score,
    sum_score_A = sum_score,
    TPM_A       = RNA_pseudobulkTPM,
    ATAC_prom_A = normalizedATAC_prom
)]
wide_B <- gs[cell_type == ct_B, .(
    TargetGene,
    n_links_B   = n_links_above_threshold,
    max_score_B = max_score,
    sum_score_B = sum_score,
    TPM_B        = RNA_pseudobulkTPM,
    ATAC_prom_B  = normalizedATAC_prom
)]

df <- merge(wide_A, wide_B, by="TargetGene", all=FALSE)

# Filter: require min TPM in at least one cell type
df <- df[pmax(TPM_A, TPM_B, na.rm=TRUE) >= opt$min_tpm]
message("Genes after TPM filter (>= ", opt$min_tpm, "): ", nrow(df))

## ============================================================================
## Compute change metrics
## ============================================================================

df[, log2FC_TPM       := log2((TPM_B + opt$tpm_pseudocount) / (TPM_A + opt$tpm_pseudocount))]
df[, delta_max_score  := max_score_B  - max_score_A]
df[, delta_n_links    := n_links_B    - n_links_A]
df[, delta_ATAC_prom  := ATAC_prom_B  - ATAC_prom_A]

# Expression direction
df[, direction := fcase(
    log2FC_TPM >=  opt$log2fc_threshold, "up",
    log2FC_TPM <= -opt$log2fc_threshold, "down",
    default = "stable"
)]
df[, direction := factor(direction, levels=c("up", "stable", "down"))]

n_up     <- sum(df$direction == "up")
n_down   <- sum(df$direction == "down")
n_stable <- sum(df$direction == "stable")
message("Expression direction: ", n_up, " up, ", n_stable, " stable, ", n_down, " down")

color_scale <- scale_color_manual(
    values = c(up=COL_UP, stable=COL_STABLE, down=COL_DOWN),
    labels = c(
        up    = sprintf("Up in %s (n=%d)", ct_B, n_up),
        stable= sprintf("Stable (n=%d)", n_stable),
        down  = sprintf("Down in %s (n=%d)", ct_B, n_down)
    ),
    name = NULL
)

## ============================================================================
## Shared scatter helper
## ============================================================================

make_scatter <- function(y_col, y_lab, subtitle=NULL) {
    # Downsample stable genes for overplotting
    df_plot <- rbind(
        df[direction != "stable"],
        df[direction == "stable"][sample(.N, min(.N, 2000))]
    )
    df_plot <- df_plot[order(direction == "stable")]  # stable points plotted first

    p <- ggplot(df_plot, aes(x=log2FC_TPM, y=get(y_col), color=direction)) +
        geom_point(size=0.6, alpha=0.5, stroke=0) +
        geom_hline(yintercept=0, color="black", linewidth=0.4, linetype="dashed") +
        geom_vline(xintercept=c(-opt$log2fc_threshold, opt$log2fc_threshold),
                   color="black", linewidth=0.3, linetype="dotted") +
        color_scale +
        labs(
            x        = sprintf("log2(TPM %s / TPM %s)", ct_B, ct_A),
            y        = y_lab,
            subtitle = subtitle
        ) +
        nature_theme +
        theme(legend.position="bottom",
              legend.key.size=unit(0.4, "lines"))

    # Label highlighted genes
    if (length(highlight_genes) > 0) {
        df_hl <- df_plot[TargetGene %in% highlight_genes]
        if (nrow(df_hl) > 0) {
            if (!requireNamespace("ggrepel", quietly=TRUE)) {
                p <- p + geom_text(
                    data=df_hl,
                    aes(label=TargetGene),
                    size=2.5, color="black", hjust=-0.1
                )
            } else {
                library(ggrepel)
                p <- p + ggrepel::geom_text_repel(
                    data=df_hl,
                    aes(label=TargetGene),
                    size=2.5, color="black",
                    max.overlaps=20, segment.size=0.3
                )
            }
        }
    }
    p
}

## ============================================================================
## Build panels
## ============================================================================

p1 <- make_scatter("delta_max_score", "Delta max E2G score (B - A)")
p2 <- make_scatter("delta_n_links",   "Delta links above threshold (B - A)")
p3 <- make_scatter("delta_ATAC_prom", "Delta promoter ATAC (B - A)")

## ============================================================================
## Combine and save
## ============================================================================

if (!requireNamespace("patchwork", quietly=TRUE)) {
    stop("patchwork is required: install.packages('patchwork')")
}
library(patchwork)

combined <- (p1 | p2 | p3) +
    plot_annotation(
        title   = sprintf("%s vs %s", ct_A, ct_B),
        theme   = theme(plot.title = element_text(size=11, color="black"))
    ) +
    plot_layout(guides="collect") & theme(legend.position="bottom")

ggsave(opt$output, combined, width=opt$width, height=opt$height)
message("Wrote: ", opt$output)

## ============================================================================
## Print summary statistics
## ============================================================================

message("\nCorrelation summary (", ct_A, " vs ", ct_B, "):")
message("  cor(log2FC_TPM, delta_max_score) = ",
        round(cor(df$log2FC_TPM, df$delta_max_score, use="complete"), 3))
message("  cor(log2FC_TPM, delta_n_links)   = ",
        round(cor(df$log2FC_TPM, df$delta_n_links,   use="complete"), 3))
message("  cor(log2FC_TPM, delta_ATAC_prom) = ",
        round(cor(df$log2FC_TPM, df$delta_ATAC_prom, use="complete"), 3))
