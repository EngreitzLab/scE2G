suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
  library(ggplot2)
})

get_stats_from_pred <- function(pred_full, pred_thresholded, score_column) {

	# stats from full predictions
	n_genes_active_promoter  <- pred_full %>% dplyr::select(TargetGene) %>% distinct() %>% nrow()
	n_genes_low_TPM <- pred_full %>% dplyr::filter(!!sym(score_column) == 0) %>%
		dplyr::select(TargetGene) %>% distinct() %>% nrow()

	# stats from thresholded predictions
	enh_gene <- pred_thresholded %>% dplyr::select(chr, start, end, TargetGene)

	# get stats
	n_enh <- enh_gene %>% dplyr::select(-c(TargetGene)) %>% distinct() %>% nrow()
	n_gene <- enh_gene %>% dplyr::select(TargetGene) %>% distinct() %>% nrow()
	n_links <- enh_gene %>% nrow()
	
	n_gene_per_enh <- enh_gene %>% group_by(chr, start, end) %>%
		tally() %>% pull(n) %>% mean()
	n_enh_per_gene <- enh_gene %>% group_by(TargetGene) %>% 
		tally() %>% pull(n) %>% mean()

	if (!("distanceToTSS" %in% colnames(pred_thresholded))) {
		pred_thresholded$distanceToTSS <- pred_thresholded$distance
	}
	mean_distance <- pred_thresholded %>% pull(distanceToTSS) %>% mean()

	mean_size <- enh_gene %>% dplyr::select(chr, start, end) %>% distinct() %>%
		mutate(width = end - start) %>%
		pull(width) %>% mean()

	res <- c(n_enh_elements = n_enh, n_genes_with_enh = n_gene, n_enh_gene_links = n_links,
		mean_genes_per_enh = n_gene_per_enh, mean_enh_per_gene = n_enh_per_gene,
		mean_dist_to_tss = mean_distance, mean_enh_width = mean_size,
		n_genes_active_promoter = n_genes_active_promoter, n_genes_not_expressed = n_genes_low_TPM) %>%
		as.list()%>% data.frame()

	return(res)
}



##### RUN

pred_full <- fread(snakemake@input$pred_full)
pred_thresholded <- fread(snakemake@input$pred_thresholded) %>% dplyr::filter(class != "promoter")
score_column <- snakemake@params$score_column
out_file <- snakemake@output$stats

res <- get_stats_from_pred(pred_full, pred_thresholded, score_column)

res$cluster <- snakemake@wildcards$cluster
res$model_name <- snakemake@wildcards$model_name
res$fragments_total <- as.numeric(readLines(snakemake@input$fragment_count))
res$cell_count <- as.numeric(readLines(snakemake@input$cell_count))

umi_path <- snakemake@input$umi_count
if (file.info(umi_path)$isdir) {
	res$umi_count = 0
} else {
	res$umi_count <- as.numeric(readLines(umi_path))
}

fwrite(res, out_file, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
