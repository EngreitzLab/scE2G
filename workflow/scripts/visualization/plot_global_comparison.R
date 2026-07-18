#!/usr/bin/env Rscript
#
# plot_global_comparison.R
#
# Double-triangle heatmap comparing cell types by E2G prediction similarity.
# Upper triangle = Jaccard (above-threshold link overlap).
# Lower triangle = Pearson correlation of scores (all links, absent = 0).
# Cell types clustered by average-linkage on Jaccard x Pearson product.
#
# Usage:
#   Rscript plot_global_comparison.R \
#     --pairwise_metrics results/.../comparison/pairwise_metrics.tsv \
#     --output results/.../comparison/global_comparison.pdf
#
# Optional --sample_key TSV (columns: biosample, display_name, dataset) adds
# display names and a dataset color bar alongside the heatmap.

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
    make_option("--pairwise_metrics", type="character", help="Path to pairwise_metrics.tsv"),
    make_option("--output",           type="character", help="Output PDF path"),
    make_option("--sample_key",       type="character", default=NULL,
                help="Optional TSV (biosample, display_name, dataset) for display names and color bar"),
    make_option("--lower_metric",     type="character", default="Pearson",
                help="Metric for lower triangle: Pearson, Pearson_log1p, Spearman [default: %default]"),
    make_option("--upper_metric",     type="character", default="jaccard",
                help="Metric for upper triangle [default: %default]"),
    make_option("--lower_limits",     type="character", default="0,0.6",
                help="Color limits for lower triangle 'min,max' [default: %default]"),
    make_option("--upper_limits",     type="character", default="0,0.6",
                help="Color limits for upper triangle 'min,max' [default: %default]"),
    make_option("--width",            type="numeric",   default=NULL,
                help="Heatmap width in inches (default: auto from n cell types)"),
    make_option("--height",           type="numeric",   default=NULL,
                help="Heatmap height in inches (default: auto from n cell types)")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$pairwise_metrics)) stop("--pairwise_metrics is required")
if (is.null(opt$output))           stop("--output is required")

parse_limits <- function(s) as.numeric(strsplit(s, ",")[[1]])
lower_limits <- parse_limits(opt$lower_limits)
upper_limits <- parse_limits(opt$upper_limits)

dir.create(dirname(opt$output), recursive=TRUE, showWarnings=FALSE)

## ============================================================================
## Color palettes (Nature)
## ============================================================================

blues   <- c("#c5e5fb", "#9bcae9", "#5496ce", "#006eae", "#00488d", "#002359")
reds    <- c("#f6ceca", "#e9a0a5", "#dc6464", "#c5373d", "#9b241c", "#730c0d")
diag_color  <- "#1c2a43"
absent_color <- "#ffffff"

map_to_colors <- function(values, gradient, limits) {
    ramp   <- colorRamp(gradient)
    scaled <- rescale(values, to=c(0, 1), from=limits)
    scaled <- squish(scaled, range=c(0, 1))
    cols   <- rgb(ramp(scaled), maxColorValue=255)
    cols[is.na(values)] <- absent_color
    cols
}

nature_theme <- theme_minimal() +
    theme(
        axis.text       = element_text(color="black", size=10),
        axis.title      = element_blank(),
        axis.text.x     = element_text(angle=60, hjust=1),
        legend.text     = element_text(color="black", size=9),
        legend.title    = element_text(color="black", size=9),
        panel.grid      = element_blank()
    )

make_colorbar <- function(gradient, limits, title) {
    cb <- data.table(y=seq(limits[1], limits[2], length.out=200), x=1)
    cb[, fill_color := map_to_colors(y, gradient, limits)]
    ggplot(cb, aes(x=x, y=y, fill=fill_color)) +
        geom_tile() +
        scale_fill_identity() +
        scale_y_continuous(
            breaks = c(limits[1], mean(limits), limits[2]),
            labels = c(limits[1], mean(limits), limits[2]),
            expand = c(0, 0)
        ) +
        labs(x=NULL, y=NULL, title=title) +
        theme_minimal() +
        theme(
            axis.text.x  = element_blank(),
            axis.text.y  = element_text(size=8, color="black"),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(color="black"),
            plot.title   = element_text(size=8, color="black", hjust=0.5),
            panel.grid   = element_blank()
        )
}

## ============================================================================
## Load data
## ============================================================================

df <- fread(opt$pairwise_metrics)

# Validate metric columns exist
for (col in c(opt$lower_metric, opt$upper_metric)) {
    if (!col %in% names(df)) stop("Column '", col, "' not found in pairwise_metrics.tsv")
}

all_samples <- unique(c(df$biosampleA, df$biosampleB))

# Apply display names from sample key if provided
display_names <- setNames(all_samples, all_samples)  # identity by default
dataset_map   <- NULL

if (!is.null(opt$sample_key)) {
    sk <- fread(opt$sample_key)
    req <- c("biosample", "display_name")
    if (!all(req %in% names(sk))) stop("sample_key must have columns: biosample, display_name")
    dn <- setNames(sk$display_name, sk$biosample)
    display_names[names(dn)] <- dn
    if ("dataset" %in% names(sk)) {
        dataset_map <- setNames(sk$dataset, sk$biosample)
    }
}

## ============================================================================
## Mirror to symmetric matrix and cluster
## ============================================================================

# Add symmetric counterparts
df_sym <- rbind(
    df[, .(biosampleA, biosampleB,
           lower = get(opt$lower_metric),
           upper = get(opt$upper_metric))],
    df[, .(biosampleA = biosampleB, biosampleB = biosampleA,
           lower = get(opt$lower_metric),
           upper = get(opt$upper_metric))]
)
df_sym[, product := lower * upper]

# Wide matrix for clustering
wide <- dcast(df_sym, biosampleA ~ biosampleB, value.var="product")
mat  <- as.matrix(wide[, -1])
rownames(mat) <- wide$biosampleA
mat[is.na(mat)] <- 0

clust <- hclust(dist(mat), method="average")
ordered_samples <- rownames(mat)[clust$order]
ordered_names   <- display_names[ordered_samples]

# Add diagonal rows
df_diag <- data.table(
    biosampleA = all_samples, biosampleB = all_samples,
    lower = NA_real_, upper = NA_real_, product = NA_real_
)
df_all <- rbind(df_sym, df_diag)

# Assign display names and factor levels
df_all[, nameA := factor(display_names[biosampleA], levels=ordered_names)]
df_all[, nameB := factor(display_names[biosampleB], levels=ordered_names)]

# Classify each cell as upper triangle, lower triangle, or diagonal
df_all[, tri := fcase(
    as.integer(nameA) > as.integer(nameB), "upper",
    as.integer(nameA) < as.integer(nameB), "lower",
    default = "diag"
)]

# Map colors
df_all[tri == "lower", fill_color := map_to_colors(lower, blues, lower_limits)]
df_all[tri == "upper", fill_color := map_to_colors(upper, reds,  upper_limits)]
df_all[tri == "diag",  fill_color := diag_color]

## ============================================================================
## Heatmap
## ============================================================================

hm <- ggplot(df_all, aes(x=nameA, y=nameB, fill=fill_color)) +
    geom_tile(color="white", linewidth=0.3) +
    scale_fill_identity() +
    coord_fixed() +
    labs(x=NULL, y=NULL) +
    nature_theme

## ============================================================================
## Combine heatmap with color bar legends (and optional dataset bar)
## ============================================================================

library(patchwork)

cb_upper <- make_colorbar(reds,  upper_limits, opt$upper_metric)
cb_lower <- make_colorbar(blues, lower_limits, opt$lower_metric)
legend_col <- cb_upper / cb_lower

# Dynamic sizing: each cell type contributes ~0.45 inches; minimum 3 inches
n_samples  <- length(all_samples)
hm_size    <- opt$width  %||% max(3, n_samples * 0.45 + 1.0)
hm_height  <- opt$height %||% hm_size
cb_width   <- 0.5  # narrow color bar column

if (!is.null(dataset_map)) {
    sk_plot <- data.table(
        biosample    = ordered_samples,
        display_name = factor(ordered_names, levels=ordered_names),
        dataset      = dataset_map[ordered_samples]
    )

    ds_levels  <- unique(sk_plot$dataset)
    ds_palette <- c("#006eae", "#c5373d", "#429130", "#ca9b23",
                    "#6e788d", "#a64691", "#49bcbc", "#f29742")
    ds_cols <- setNames(ds_palette[seq_along(ds_levels)], ds_levels)

    ds_bar <- ggplot(sk_plot, aes(x=1, y=display_name, fill=dataset)) +
        geom_tile() +
        scale_fill_manual(values=ds_cols, name="Dataset") +
        theme_minimal() +
        theme(
            axis.text    = element_blank(),
            axis.title   = element_blank(),
            axis.ticks   = element_blank(),
            legend.text  = element_text(size=8, color="black"),
            legend.title = element_text(size=8, color="black"),
            panel.grid   = element_blank()
        )

    combined <- hm + ds_bar + legend_col +
        plot_layout(ncol=3, widths=c(hm_size, 0.3, cb_width))
    ggsave(opt$output, combined, width=hm_size + 0.3 + cb_width, height=hm_height)
} else {
    combined <- hm + legend_col + plot_layout(ncol=2, widths=c(hm_size, cb_width))
    ggsave(opt$output, combined, width=hm_size + cb_width, height=hm_height)
}

message("Wrote: ", opt$output)
