#!/usr/bin/env Rscript
#
# plot_gene_comparison.R
#
# Bubble matrix comparing the enhancer landscape for a single gene across cell types.
# Rows = consensus elements (ordered by genomic position).
# Columns = cell types (ordered by --cell_types or alphabetically).
# Dot size = normalizedATAC_enh (accessibility at element).
# Dot color = E2G score (purple gradient; grey if absent/not called).
# Dot border = black if above threshold, absent if not called.
# Top strip = RNA_pseudobulkTPM and normalizedATAC_prom per cell type.
#
# Usage:
#   Rscript plot_gene_comparison.R \
#     --gene GATA2 \
#     --comparison_dir results/.../comparison \
#     --threshold 0.177 \
#     --output results/.../comparison/gene_plots/GATA2.pdf

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
    make_option("--gene",            type="character", help="Target gene name"),
    make_option("--comparison_dir",  type="character", help="Directory with precompute outputs"),
    make_option("--output",          type="character", default=NULL,
                help="Output PDF path (default: {comparison_dir}/gene_plots/{gene}.pdf)"),
    make_option("--threshold",       type="numeric",   default=0.177),
    make_option("--cell_types",      type="character", default=NULL,
                help="Comma-separated ordered cell types (default: alphabetical from data)"),
    make_option("--max_elements",    type="integer",   default=75L,
                help="Max elements to show; keeps highest-scoring [default: %default]"),
    make_option("--width",           type="numeric",   default=NULL),
    make_option("--height",          type="numeric",   default=NULL),
    make_option("--bw_dir",          type="character", default=NULL,
                help="Directory containing {cell_type}/ATAC_norm.bw files for coverage tracks"),
    make_option("--locus",           type="character", default=NULL,
                help="Override locus window as chrX:start-end (default: span all above-threshold links)")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$gene))           stop("--gene is required")
if (is.null(opt$comparison_dir)) stop("--comparison_dir is required")

scores_path <- file.path(opt$comparison_dir, "consensus_scores.tsv.gz")
if (!file.exists(scores_path)) stop("consensus_scores.tsv.gz not found in: ", opt$comparison_dir)

if (is.null(opt$output)) {
    out_dir <- file.path(opt$comparison_dir, "gene_plots")
    dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
    opt$output <- file.path(out_dir, paste0(opt$gene, ".pdf"))
} else {
    dir.create(dirname(opt$output), recursive=TRUE, showWarnings=FALSE)
}

## ============================================================================
## Nature theme and palettes
## ============================================================================

purples      <- c("#e9d3ea", "#d3a9ce", "#b778b3", "#a64791", "#792374", "#430b4e")
absent_color <- "#c5c1a5"   # stone — element not called in this cell type
below_color  <- "#e5e5e9"   # light grey — called but below threshold

`%||%` <- function(a, b) if (!is.null(a)) a else b

nature_theme <- theme_classic() +
    theme(
        axis.text   = element_text(color="black"),
        axis.ticks  = element_line(color="black"),
        axis.line   = element_line(color="black"),
        strip.text  = element_text(color="black"),
        panel.grid  = element_blank(),
        legend.text = element_text(color="black"),
        legend.title= element_text(color="black")
    )

## ============================================================================
## Load data for this gene
## ============================================================================

message("Loading data for gene: ", opt$gene)

# Efficiently extract only rows for this gene using awk
cmd <- sprintf("zcat '%s' | awk -F'\\t' 'NR==1 || $2==\"%s\"'", scores_path, opt$gene)
dt  <- fread(cmd=cmd)

if (nrow(dt) == 0) stop("Gene '", opt$gene, "' not found in consensus_scores.tsv.gz")

# Exclude self-promoter rows
if ("isSelfPromoter" %in% names(dt)) {
    dt <- dt[isSelfPromoter == FALSE]
}

if (nrow(dt) == 0) stop("No non-self-promoter predictions found for gene '", opt$gene, "'")

message("  Found ", uniqueN(dt$consensus_id), " elements across ",
        uniqueN(dt$cell_type), " cell types")

## ============================================================================
## Build complete (element x cell_type) grid
## ============================================================================

cell_types <- if (!is.null(opt$cell_types)) {
    trimws(strsplit(opt$cell_types, ",")[[1]])
} else {
    sort(unique(dt$cell_type))
}

all_elements <- unique(dt[, .(consensus_id)])

# Limit to top elements if too many
if (nrow(all_elements) > opt$max_elements) {
    top_elements <- dt[, .(max_score = max(score, na.rm=TRUE)), by=consensus_id
                       ][order(-max_score)][seq_len(opt$max_elements), consensus_id]
    all_elements <- data.table(consensus_id = top_elements)
    message("  Showing top ", opt$max_elements, " elements by max score across cell types")
}

# Complete grid
grid <- CJ(consensus_id = all_elements$consensus_id, cell_type = cell_types)
grid <- merge(grid, dt, by=c("consensus_id", "cell_type"), all.x=TRUE)

# Mark absent elements (not called as a peak in this cell type)
grid[, is_absent := is.na(score)]
grid[is_absent == TRUE,  score              := 0]
grid[is_absent == TRUE,  normalizedATAC_enh := 0]
grid[is_absent == TRUE,  is_above_threshold := FALSE]

# Add element coordinates for genomic ordering
elem_coords <- unique(dt[, .(consensus_id)])
elem_coords[, consensus_chr   := sub(":.*", "", consensus_id)]
elem_coords[, consensus_start := as.integer(sub(".*:(\\d+)-.*", "\\1", consensus_id))]
elem_coords[, consensus_end   := as.integer(sub(".*-", "", consensus_id))]
elem_coords[, elem_mid := (consensus_start + consensus_end) / 2L]
setorder(elem_coords, consensus_chr, consensus_start)
elem_coords[, element_label := paste0("E", seq_len(.N))]
elem_coords[, element_label := factor(element_label, levels=rev(element_label))]  # top = most upstream

grid <- merge(grid, elem_coords[, .(consensus_id, element_label, elem_mid)], by="consensus_id")
grid[, cell_type := factor(cell_type, levels=cell_types)]

# Categorize each point for color encoding
grid[, point_type := fcase(
    is_absent          == TRUE, "absent",
    is_above_threshold == TRUE, "above",
    default = "below"
)]
grid[, point_type := factor(point_type, levels = c("absent", "below", "above"))]

## ============================================================================
## Get TSS position for reference line
## ============================================================================

tss_pos <- dt[!is.na(TargetGeneTSS), TargetGeneTSS[1]]
tss_label <- if (!is.na(tss_pos)) {
    # Find closest element by midpoint for x-axis reference
    elem_coords[which.min(abs(elem_mid - tss_pos)), element_label]
} else NULL

## ============================================================================
## Bubble matrix
## ============================================================================

# Scale ATAC size: use rank-based scaling to handle skewed distributions
atac_max <- max(grid$normalizedATAC_enh, na.rm=TRUE)
if (atac_max == 0) atac_max <- 1

# Score color scale: grey for absent/below, purple for above
score_range <- range(grid[point_type == "above", score], na.rm=TRUE)
if (length(score_range) == 0 || all(is.na(score_range))) score_range <- c(0, 1)

p_bubbles <- ggplot(grid, aes(x=cell_type, y=element_label)) +
    # Absent: small open circle — element not called as a peak in this cell type
    geom_point(
        data = grid[point_type == "absent"],
        aes(color = point_type),
        size=0.5, shape=1, stroke=0.4
    ) +
    # Below threshold: grey filled circle — element called but E2G score below threshold
    geom_point(
        data = grid[point_type == "below"],
        aes(size = normalizedATAC_enh, color = point_type),
        shape=21, fill=below_color, stroke=0.3
    ) +
    # Above threshold: colored by score, size by ATAC, black border
    geom_point(
        data = grid[point_type == "above"],
        aes(size = normalizedATAC_enh, fill = score, color = point_type),
        shape=21, stroke=0.5
    ) +
    scale_fill_gradientn(
        colors  = purples,
        limits  = score_range,
        oob     = squish,
        name    = "E2G score",
        guide   = guide_colorbar(barwidth=0.8, barheight=4, order=1)
    ) +
    scale_size_continuous(
        range = c(0.8, 6),
        limits = c(0, atac_max),
        name   = "ATAC signal\n(normalized)",
        guide  = guide_legend(
            override.aes = list(fill="grey40", color="black", shape=21),
            order = 2
        )
    ) +
    scale_color_manual(
        values = c("absent" = absent_color, "below" = "grey60", "above" = "black"),
        labels = c(
            "absent" = "Not called as peak",
            "below"  = "Peak called; score below threshold\n(size = ATAC signal)",
            "above"  = "Peak called; score above threshold\n(size = ATAC signal)"
        ),
        name  = "Element status",
        guide = guide_legend(
            override.aes = list(
                shape  = c(1,           21,          21),
                fill   = c(NA,          below_color, purples[4]),
                size   = c(2,           3,           3),
                stroke = c(0.4,         0.3,         0.5)
            ),
            order = 3
        )
    ) +
    labs(x=NULL, y=NULL) +
    nature_theme +
    theme(
        axis.text.x     = element_text(angle=45, hjust=1),
        legend.position = "right"
    )

## ============================================================================
## Top strip: TPM and ATAC_prom per cell type
## ============================================================================

# Get gene-level metrics (same for all elements, take first)
gene_meta <- dt[, .(
    TPM       = first(RNA_pseudobulkTPM),
    ATAC_prom = first(normalizedATAC_prom)
), by=cell_type]

# Fill missing cell types with 0
missing_ct <- setdiff(cell_types, gene_meta$cell_type)
if (length(missing_ct) > 0) {
    gene_meta <- rbind(gene_meta,
                       data.table(cell_type=missing_ct, TPM=0, ATAC_prom=0))
}
gene_meta[, cell_type := factor(cell_type, levels=cell_types)]

strip_data <- melt(gene_meta, id.vars="cell_type",
                   measure.vars=c("TPM", "ATAC_prom"),
                   variable.name="metric", value.name="value")
strip_data[, metric := factor(metric,
    levels=c("TPM", "ATAC_prom"),
    labels=c("RNA TPM", "ATAC prom."))]

p_strip <- ggplot(strip_data, aes(x=cell_type, y=value)) +
    geom_col(width=0.7, fill="#6e788d") +
    facet_wrap(~metric, scales="free_y", nrow=1) +
    labs(x=NULL, y=NULL, title=opt$gene) +
    nature_theme +
    theme(
        axis.text.x     = element_text(angle=45, hjust=1),
        strip.background= element_blank(),
        plot.title      = element_text(face="italic", size=11, color="black")
    )

## ============================================================================
## Locus plot
## ============================================================================

## Parse coordinates for ALL elements of this gene (not filtered to top-N)
all_elem_coords <- unique(dt[, .(consensus_id)])
all_elem_coords[, consensus_chr   := sub(":.*", "", consensus_id)]
all_elem_coords[, consensus_start := as.integer(sub(".*:(\\d+)-.*", "\\1", consensus_id))]
all_elem_coords[, consensus_end   := as.integer(sub(".*-", "", consensus_id))]
all_elem_coords[, elem_mid        := (consensus_start + consensus_end) / 2L]

## Window: span all above-threshold links + TSS, with small padding (or user-supplied coords)
if (!is.null(opt$locus)) {
    locus_chr <- sub(":.*", "", opt$locus)
    coords    <- sub(".*:", "", opt$locus)
    window_lo <- as.integer(sub("-.*", "", coords))
    window_hi <- as.integer(sub(".*-", "", coords))
} else {
    above_thresh_coords <- all_elem_coords[
        consensus_id %in% unique(dt[is_above_threshold == TRUE, consensus_id])
    ]
    locus_chr <- above_thresh_coords$consensus_chr[1]
    window_lo <- min(above_thresh_coords$consensus_start, if (!is.na(tss_pos)) tss_pos else Inf,  na.rm=TRUE)
    window_hi <- max(above_thresh_coords$consensus_end,   if (!is.na(tss_pos)) tss_pos else -Inf, na.rm=TRUE)
    padding   <- max((window_hi - window_lo) * 0.05, 2000L)
    window_lo <- window_lo - padding
    window_hi <- window_hi + padding
}

## Elements in window
window_elems <- all_elem_coords[
    consensus_chr == locus_chr & consensus_start >= window_lo & consensus_end <= window_hi
]

## Full (element x cell_type) grid for window, with absent fill
locus_grid <- CJ(consensus_id = window_elems$consensus_id, cell_type = cell_types)
locus_grid <- merge(locus_grid, dt, by=c("consensus_id", "cell_type"), all.x=TRUE)
locus_grid[, is_absent := is.na(score)]
locus_grid[is_absent == TRUE, score              := 0]
locus_grid[is_absent == TRUE, normalizedATAC_enh := 0]
locus_grid[is_absent == TRUE, is_above_threshold := FALSE]
locus_grid[, point_type := factor(fcase(
    is_absent          == TRUE, "absent",
    is_above_threshold == TRUE, "above",
    default = "below"
), levels = c("absent", "below", "above"))]
locus_grid <- merge(locus_grid, window_elems, by="consensus_id")
locus_grid[, cell_type := factor(cell_type, levels=cell_types)]
locus_grid[, is_top_n  := consensus_id %in% all_elements$consensus_id]

## Normalize ATAC 0-1 per cell type and normalize x coords 0-1 within window
locus_grid[, atac_max := max(normalizedATAC_enh, na.rm=TRUE), by=cell_type]
locus_grid[, atac_norm    := fifelse(atac_max > 0, normalizedATAC_enh / atac_max, 0)]
locus_grid[, x_start_norm := (consensus_start - window_lo) / (window_hi - window_lo)]
locus_grid[, x_end_norm   := (consensus_end   - window_lo) / (window_hi - window_lo)]
locus_grid[, x_mid_norm   := (elem_mid        - window_lo) / (window_hi - window_lo)]
locus_grid[, atac_max := NULL]

tss_x_norm <- if (!is.na(tss_pos)) (tss_pos - window_lo) / (window_hi - window_lo) else NA_real_

## x-axis labels: convert 0-1 back to genomic kb
x_label_fn <- function(x) {
    pos_bp <- x * (window_hi - window_lo) + window_lo
    paste0(format(round(pos_bp / 1e3), big.mark=",", scientific=FALSE), " kb")
}

## Load bigwig coverage tracks if --bw_dir provided
bw_coverage <- NULL
if (!is.null(opt$bw_dir)) {
    suppressPackageStartupMessages(library(rtracklayer))
    query_gr <- GRanges(locus_chr, IRanges(as.integer(window_lo), as.integer(window_hi)))
    bw_list <- lapply(cell_types, function(ct) {
        bw_path <- file.path(opt$bw_dir, ct, "ATAC_norm.bw")
        if (!file.exists(bw_path)) {
            message("  No bigwig found for cell type: ", ct)
            return(NULL)
        }
        tryCatch({
            gr <- import.bw(bw_path, which=query_gr, as="GRanges")
            dt <- as.data.table(as.data.frame(gr)[, c("seqnames","start","end","score")])
            dt[, cell_type := ct]
            dt
        }, error = function(e) {
            message("  Failed to load bigwig for ", ct, ": ", conditionMessage(e))
            NULL
        })
    })
    bw_coverage <- rbindlist(bw_list[!sapply(bw_list, is.null)])
    if (nrow(bw_coverage) > 0) {
        bw_coverage[, cell_type := factor(cell_type, levels=cell_types)]
        bw_coverage[, cov_max := max(score, na.rm=TRUE), by=cell_type]
        bw_coverage[, cov_norm := fifelse(cov_max > 0, score / cov_max, 0)]
        bw_coverage[, cov_max := NULL]
        bw_coverage[, x_start_norm := (start - window_lo) / (window_hi - window_lo)]
        bw_coverage[, x_end_norm   := (end   - window_lo) / (window_hi - window_lo)]
        message("  Loaded bigwig coverage for ", uniqueN(bw_coverage$cell_type), " cell types")
    }
}

## Build locus plot: ATAC bigwig coverage + arcs only (no element bars)
p_locus <- ggplot()

## Bigwig coverage track in Nature blue
if (!is.null(bw_coverage) && nrow(bw_coverage) > 0) {
    p_locus <- p_locus +
        geom_rect(
            data = bw_coverage,
            aes(xmin=x_start_norm, xmax=x_end_norm, ymin=0, ymax=cov_norm),
            fill="#00488d", color=NA, inherit.aes=FALSE
        )
}

## Arcs: above-threshold links only; purple = cell-type variable, grey = shared across all cell types
locus_arc <- locus_grid[point_type == "above" & !is.na(tss_x_norm)]
if (nrow(locus_arc) > 0) {
    locus_arc[, tss_x := tss_x_norm]
    score_lim <- range(locus_arc$score, na.rm=TRUE)

    ## Classify element as variable (above threshold in fewer than all shown cell types) vs shared
    n_ct_shown <- length(cell_types)
    element_n_above <- locus_grid[, .(n_above = sum(is_above_threshold == TRUE, na.rm=TRUE)),
                                   by=consensus_id]
    locus_arc <- merge(locus_arc, element_n_above, by="consensus_id", all.x=TRUE)
    locus_arc[, arc_type := factor(
        fifelse(n_above < n_ct_shown, "variable", "shared"),
        levels = c("shared", "variable")
    )]

    arc_left  <- locus_arc[x_mid_norm <  tss_x_norm]
    arc_right <- locus_arc[x_mid_norm >= tss_x_norm]

    if (nrow(arc_left) > 0) {
        p_locus <- p_locus + geom_curve(
            data = arc_left,
            aes(x=x_mid_norm, xend=tss_x, y=0, yend=0, color=arc_type, linewidth=score),
            curvature=-0.3, ncp=10
        )
    }
    if (nrow(arc_right) > 0) {
        p_locus <- p_locus + geom_curve(
            data = arc_right,
            aes(x=x_mid_norm, xend=tss_x, y=0, yend=0, color=arc_type, linewidth=score),
            curvature=0.3, ncp=10
        )
    }

    p_locus <- p_locus +
        scale_color_manual(
            values = c("shared"="#96a0b3", "variable"="#792374"),
            labels = c("shared"="Shared (all cell types)", "variable"="Cell-type variable"),
            name   = "Link type",
            guide  = guide_legend(order=1)
        ) +
        scale_linewidth_continuous(
            range  = c(0.2, 0.6),
            limits = score_lim,
            breaks = pretty(score_lim, n=3),
            name   = "E2G score",
            guide  = guide_legend(
                override.aes = list(color="black", linewidth=c(0.2, 0.4, 0.6)),
                order=2
            )
        )
}

## TSS line (gene label shown in gene track above, not per-facet)
if (!is.na(tss_x_norm)) {
    p_locus <- p_locus +
        geom_vline(xintercept=tss_x_norm, linetype="dashed", color="black", linewidth=0.25)
}

p_locus <- p_locus +
    facet_wrap(~cell_type, ncol=1, strip.position="right") +
    scale_x_continuous(labels=x_label_fn, limits=c(0, 1), expand=c(0.01, 0)) +
    labs(x=locus_chr, y="ATAC signal (norm.)",
         title=sprintf("%s locus - top %d of %d elements outlined",
                       opt$gene, nrow(all_elements), uniqueN(dt$consensus_id))) +
    scale_y_continuous(breaks=c(0, 0.5, 1)) +
    nature_theme +
    theme(
        strip.text.y.right = element_text(angle=0, hjust=0),
        strip.background   = element_blank(),
        panel.spacing      = unit(0.1, "lines"),
        legend.position    = "right"
    )

## ============================================================================
## Distance from TSS plot
## ============================================================================

has_tss <- !is.na(tss_pos)

if (has_tss) {
    dist_data <- elem_coords[consensus_id %in% all_elements$consensus_id]
    dist_data[, dist_kb := (elem_mid - tss_pos) / 1000]

    p_distance <- ggplot(dist_data, aes(x = dist_kb, y = element_label)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_segment(aes(xend = 0, yend = element_label), color = "grey70", linewidth = 0.3) +
        geom_point(shape = 21, fill = "#96a0b3", color = "black", size = 2, stroke = 0.4) +
        labs(x = "Distance from TSS (kb)", y = NULL) +
        nature_theme +
        theme(
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y  = element_blank()
        )
}

## ============================================================================
## Combine and save
## ============================================================================

if (!requireNamespace("patchwork", quietly=TRUE)) {
    stop("patchwork is required: install.packages('patchwork')")
}
library(patchwork)

n_elements  <- uniqueN(grid$element_label)
n_ct        <- length(cell_types)
plot_height <- opt$height %||% max(5, n_elements * 0.25 + 3)

if (has_tss) {
    p_bubbles_main <- p_bubbles + theme(
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank()
    )
    combined <- plot_spacer() + p_strip + p_distance + p_bubbles_main +
        plot_layout(ncol = 2, heights = c(1, 4), widths = c(1, n_ct))
    plot_width <- opt$width %||% max(7, n_ct * 1.2 + 5)
} else {
    combined   <- p_strip / p_bubbles + plot_layout(heights = c(1, 4))
    plot_width <- opt$width %||% max(5, n_ct * 1.2 + 3)
}

ggsave(opt$output, combined, width=plot_width, height=plot_height)
message("Wrote: ", opt$output)

## Gene track panel: TSS marker + gene name, assembled above the locus facets
if (!is.na(tss_x_norm)) {
    p_gene_track <- ggplot() +
        annotate("segment", x=0, xend=1, y=0, yend=0, color="grey70", linewidth=0.4) +
        annotate("point",   x=tss_x_norm, y=0, shape=25, fill="black", color="black", size=2.5) +
        annotate("text",    x=tss_x_norm, y=0.7, label=opt$gene,
                 fontface="italic", size=3, hjust=0.5, color="black") +
        scale_x_continuous(limits=c(0, 1), expand=c(0.01, 0)) +
        scale_y_continuous(limits=c(-0.3, 1.2)) +
        labs(x=NULL, y=" ") +   # blank y label to match p_locus left-axis width
        theme_classic() +
        theme(
            axis.text    = element_blank(),
            axis.ticks   = element_blank(),
            axis.line    = element_blank(),
            axis.title.y = element_text(color=NA),  # invisible spacer
            panel.grid   = element_blank(),
            plot.margin  = margin(4, 5, 0, 5)
        )
    p_locus_out <- p_gene_track / p_locus +
        plot_layout(heights=c(1, length(cell_types) * 5))
} else {
    p_locus_out <- p_locus
}

locus_output <- sub("\\.pdf$", "_locus.pdf", opt$output)
locus_height <- length(cell_types) * 1.0 + 1.8
ggsave(locus_output, p_locus_out, width=plot_width, height=locus_height)
message("Wrote: ", locus_output)
