## Create Activity-only feature table for an ABC sample

# required packages and functions
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# load inputs
abc_gene_list <- fread(snakemake@input$abc_gene_list, sep = "\t")
abc_element_list <- fread(snakemake@input$abc_element_list, sep = "\t")
pred <- fread(snakemake@input$prediction_file)
this_cluster <- snakemake@wildcards$cluster
tpm_threshold <- snakemake@params$tpm_threshold

genes_out <- snakemake@output$gene_list
elements_out <- snakemake@output$element_list

# all pred columns: chr	start	end	name	class	TargetGene	TargetGeneTSS	TargetGeneIsExpressed	TargetGeneEnsembl_ID
# isSelfPromoter	CellType	distance	normalizedATAC_prom	ABC.Score	ABC.Score.1	numCandidateEnhGene	numTSSEnhGene	numNearbyEnhancers
# ubiqExpressed	RNA_meanLogNorm	RNA_pseudobulkTPM	RNA_percentCellsDetected	Kendall	ARC.E2G.Score	E2G.Score	E2G.Score.qnorm	E2G.Score.qnorm.ignoreTPM


## update gene list
# ABC gene list columns: chr	start	end	name	score	strand	Ensembl_ID	gene_type	symbol	tss	Expression	is_ue	cellType
# ATAC.tagAlign.sort.gz.readCount	...	PromoterActivityQuantile

gene_columns <- c("normalizedATAC_prom", "ubiqExpressed")

pred_genes <- dplyr::select(pred, TargetGene, TargetGeneEnsembl_ID, any_of(gene_columns), CellType) %>%
	distinct() %>%
	rename(name = TargetGene, Ensembl_ID = TargetGeneEnsembl_ID)

gene_list <- abc_gene_list %>%
	dplyr::select(-c(is_ue, Expression, cellType), -contains("RPKM")) %>%
	mutate(removed_by_promoter_activity = !(name %in% pred_genes$name)) %>% # refers to promoter activity quantile
	left_join(pred_genes, by = c("name", "Ensembl_ID"))

if ("RNA_pseudobulkTPM" %in% colnames(pred)) { # we have RNA data
	RNA_columns <- c("RNA_meanLogNorm", "RNA_pseudobulkTPM", "RNA_percentCellsDetected")
	gex <- fread(snakemake@input$gene_expr_file) %>%
		select(name = TargetGene, any_of(RNA_columns))

	gene_list <- left_join(gene_list, gex, by = "name") %>% 
		mutate(removed_by_TPM = (RNA_pseudobulkTPM < tpm_threshold),
		considered_for_predictions = !removed_by_promoter_activity & !removed_by_TPM)
} else {
	gene_list <- mutate(gene_list, considered_for_predictions = !removed_by_promoter_activity)
}
fwrite(gene_list, genes_out, sep = "\t", col.names = TRUE, row.names = FALSE,, quote = FALSE)

## update enhancer list
# ABC enhancer list coluns: chr	start	end	ATAC.tagAlign.sort.gz.readCount	...	ATAC.RPMnonpromoter.quantile	activity_base	activity_base_no_qnorm name class
element_columns <- c("numNearbyEnhancers", "sumNearbyEnhancers", "normalizedATAC_enh")

pred_elements <- dplyr::select(pred, chr, start, end, any_of(element_columns)) %>%
	distinct()

element_list <- left_join(abc_element_list, pred_elements, by = c("chr", "start", "end"))
fwrite(element_list, elements_out, sep = "\t", col.names = TRUE, row.names = FALSE,, quote = FALSE)


