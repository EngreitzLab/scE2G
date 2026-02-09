# required packages
library(data.table)
library(dplyr)

# inputs from snakemake
model_dirs = snakemake@params$model_dirs %>% strsplit(" ") %>% unlist()

for (i in 1:length(model_dirs)){
	ft = fread(file.path(model_dirs[i],"feature_table.tsv"))
    if (i==1){df = ft} else {df = rbind(df, ft)}
}

# for sc-E2G pipeline
if (("ARC.E2G.Score" %in% df$feature) | ("Kendall" %in% df$feature)){
	ARC_rows = data.frame(c("RNA_meanLogNorm", "RNA_pseudobulkTPM", "RNA_percentCellsDetected", "Kendall", "ARC.E2G.Score", "ABC.Score", "normalizedATAC_enh"),
		c("mean_log_normalized_rna", "RnaPseudobulkTPM", "RnaDetectedPercent", "Kendall", "ARC.E2G.Score", "ABC.Score", "normalized_atac_enh"),
		c(NA, NA, NA, NA, NA, NA, NA),
		c("mean", "mean", "mean", "max", "mean", "sum", "mean"),
		c(0, 0, 0, 0, 0, 0, 0),
		c("Mean log normalized RNA expression", "RNA pseudobulk TPM", "RNA percent cells detected", "Kendall correlation", "ARC-E2G score", "ABC score", "ATAC signal at E"));

	colnames(ARC_rows) = colnames(df)

	# handle ABC.Score & powerlaw.Score redundancy
	if ("ABC.Score" %in% df$feature) {
		ARC_rows <- dplyr::filter(ARC_rows, feature != "ABC.Score")
	}
	
	df = rbind(df, ARC_rows)
}

df = dplyr::distinct(df)
fwrite(df, file = snakemake@output$biosample_features, sep = "\t")
