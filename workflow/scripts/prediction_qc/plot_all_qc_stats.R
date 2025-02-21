suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
})

assign_model_colors <- function(model_names) {
  # Define specific model colors
  predefined_colors <- c(
    "multiome_powerlaw_v2" = "#792374",
    "scATAC_powerlaw_v2" = "#006479",
    "multiome_megamap_v2" = "#00488d",
	"multiome_powerlaw_v3" = "#792374",
    "scATAC_powerlaw_v3" = "#006479",
    "multiome_megamap_v3" = "#00488d"
  )

  # Define available fallback colors
  fallback_colors <- c("#429130", "#c5373d", "#e96a00", "#ca9b23")
  
  # Create an empty named vector
  model_colors <- setNames(character(length(model_names)), model_names)
  
  # Assign predefined colors if present
  for (model in names(predefined_colors)) {
    if (model %in% model_names) {
      model_colors[model] <- predefined_colors[model]
    }
  }
  
  # Get models that do not have predefined colors
  remaining_models <- model_names[model_colors == ""]
  
  # Assign fallback colors
  if (length(remaining_models) > 0) {
    assigned_fallbacks <- head(fallback_colors, length(remaining_models))
    model_colors[remaining_models] <- assigned_fallbacks
  }
  
  # Assign random colors if needed
  uncolored_models <- names(model_colors[model_colors == ""])
  if (length(uncolored_models) > 0) {
    random_colors <- colors()[sample(length(colors()), length(uncolored_models))]
    model_colors[uncolored_models] <- random_colors
  }
  
  return(model_colors)
}

plot_dataset_stats <- function(stats, label_key, out_file) {
	# 1) frag/cell by log10(cells)
	n_label <- paste0("N = ", nrow(stats))
	x_label <- paste0(label_key["cell_count"], "\n", n_label)

	s1 <- ggplot(stats, aes(x = cell_count, y = frag_per_cell)) +
		geom_point(shape = 16, size = 3, color = "#1c2a43") +
		labs(x = x_label, y = label_key["frag_per_cell"]) +
		scale_x_log10() +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), axis.ticks = element_line(color = "#000000"),
			aspect.ratio = 1, legend.position = "None")

	# 2) umi/cell by log10(cells)
	s2 <- ggplot(stats, aes(x = cell_count, y = umi_per_cell)) +
		geom_point(shape = 16, size = 3, color = "#1c2a43") +
		labs(x = x_label, y = label_key["umi_per_cell"]) +
		scale_x_log10() +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), axis.ticks = element_line(color = "#000000"),
			aspect.ratio = 1, legend.position = "None")
	
	# combine, arrange, and save
	plot_list <- list(s1, s2)
	combined_plot <- ggarrange(plotlist = plot_list, nrow =1, ncol = 2)

	ht <- 4
	wd <- 7
	ggsave(out_file, combined_plot, width = wd, height = ht, device = "pdf")

}

plot_scatter_set <- function(stats, x_vars, x_thresh, y_vars, colors, label_key, out_file) {
  	plot_list <- list()
  
 	 for (y in y_vars) {
    	for (i in seq_along(x_vars)) {
			x <- x_vars[i]
			x_line <- x_thresh[i]
			x_max <- max(x_line + 1, max(stats[[x]]) + 1)
			
			p <- ggplot(stats, aes(x = !!sym(x), y = !!sym(y))) +
				geom_vline(xintercept = x_line, linetype = "dashed", color = "#96a0b3") +
				geom_point(aes(color = model_name), shape = 16, size = 3) +
				labs(x = label_key[x], y = label_key[y]) +
				scale_color_manual(values = colors) +
				scale_x_log10(limits = c(1e-10, x_max)) +
				theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), axis.ticks = element_line(color = "#000000"),
					aspect.ratio = 1, legend.position = "bottom")

      	plot_list[[paste(y, x)]] <- p
    	}	
  	}
  
	# Arrange plots in a grid, add a shared legend, and save
	combined_plot <- ggarrange(plotlist = plot_list, nrow = length(y_vars), ncol = length(x_vars), common.legend = TRUE, legend = "bottom")

	ht <- length(y_vars) * 3 + 1
	wd <- length(x_vars) * 3 + 1
	ggsave(out_file, combined_plot, width = wd, height = ht, device = "pdf")
}

plot_violin_metrics <- function(stats, label_key, colors, out_file) {
  
	plots <- list()  # probably ~9
	
	for (i in seq_along(label_key)) {
		metric_col <- names(label_key)[i]
		axis_label <- label_key[metric_col]

		p <- ggplot(stats, aes(y = model_name, x = !!sym(metric_col), fill = model_name)) +
		geom_violin(trim = TRUE, scale = "count", linewidth = 0, alpha = 0.7, ) + 
		geom_jitter(aes(color = model_name), height = 0.2, width = 0, shape = 16, size = 2, alpha = 1) +  
		scale_fill_manual(values = colors) + scale_color_manual(values = colors) + 
		labs(x = axis_label, y = "") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), axis.ticks = element_line(color = "#000000"),
			axis.text.y=element_blank(), legend.position = "none")
		
		plots[[i]] <- p  # Add plot to the list
	}
	
	# Combine all plots into a grid with common legend
	ncol <- 3
	nrow <- ceiling(length(label_key) / ncol)
	
	combined_plot <- ggarrange(plotlist = plots, ncol = 3, nrow = nrow, common.legend = TRUE,  legend = "bottom")

	ht <- nrow * 3 + 1
	wd <- ncol * 3 + 1
	ggsave(out_file, plot = combined_plot, width = wd, height = ht, device = "pdf")
  }

##### MAIN

## read in and combine stats files
stats_files <- snakemake@input$stats_files %>% strsplit(" ") %>% unlist()
stats <- lapply(stats_files, fread) %>% rbindlist() %>% as_tibble()
fwrite(stats, snakemake@output$all_stats, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

model_names <- stats$model_name %>% unique()
cp <- assign_model_colors(model_names)

## make plots
# some params
label_key <- c(fragments_total = "# unique ATAC fragments in cluster", cell_count = "# cells in cluster", umi_count = "# RNA UMIs in cluster",
	frag_per_cell = "Mean unique ATAC fragments per cell", umi_per_cell = "Mean RNA UMIs per cell",
	n_enh_gene_links = "# enhancer-gene links", n_enh_elements = "# unique enhancers",
	mean_genes_per_enh = "Mean # genes per enhancer", mean_enh_per_gene = "Mean # enhancers per gene",
	n_genes_with_enh = "# genes with 1+ enhancer", n_genes_active_promoter = "# genes with accessible promoter", n_genes_not_expressed = "# genes below TPM threshold",
	mean_dist_to_tss = "Mean distance to TSS (bp)", mean_enh_width = "Mean width of enhancer element (bp)")

# output files
ds_stats_out <- snakemake@output$ds_stats
enh_stats_out <- snakemake@output$enh_stats
eg_stats_out <- snakemake@output$eg_stats
gene_stats_out <- snakemake@output$gene_stats
dist_stats_out <- snakemake@output$dist_size_stats
all_distributions_out <- snakemake@output$all_distributions

enh_vars <- c("n_enh_gene_links", "n_enh_elements")
eg_vars <- c("mean_genes_per_enh", "mean_enh_per_gene")
gene_vars <- c("n_genes_with_enh", "n_genes_active_promoter")
dist_vars <- c("mean_dist_to_tss", "mean_enh_width")

if (sum(stats$umi_count) > 0) { # plot metrics by UMI count
	x_vars = c("fragments_total", "cell_count", "umi_count")
	gene_vars <- c(gene_vars, "n_genes_not_expressed")
	
} else {
	x_vars = c("fragments_total", "cell_count")
}

# plot stats about datasets
stats_ds <- dplyr::select(stats, cluster, fragments_total, cell_count, umi_count) %>% distinct() %>%
	mutate(frag_per_cell = fragments_total / cell_count,
		umi_per_cell = umi_count / cell_count)
plot_dataset_stats(stats = stats_ds, label_key = label_key, out_file = ds_stats_out)

# plot scatter plots of each metric versus sequencing depth stats
x_thresh = c(2e6, 100, 1e6)
plot_scatter_set(stats = stats, x_vars = x_vars, x_thresh = x_thresh, y_vars = enh_vars, colors = cp, label_key = label_key, out_file = enh_stats_out)
plot_scatter_set(stats = stats, x_vars = x_vars, x_thresh = x_thresh, y_vars = eg_vars, colors = cp, label_key = label_key, out_file = eg_stats_out)
plot_scatter_set(stats = stats, x_vars = x_vars, x_thresh = x_thresh, y_vars = gene_vars, colors = cp, label_key = label_key, out_file = gene_stats_out)
plot_scatter_set(stats = stats, x_vars = x_vars, x_thresh = x_thresh, y_vars = dist_vars, colors = cp, label_key = label_key, out_file = dist_stats_out)

# violin plots of each metric (each model a diff color) 
violin_label_key <- label_key[6:length(label_key)]
plot_violin_metrics(stats = stats, label_key = violin_label_key, colors = cp, out_file = all_distributions_out)









