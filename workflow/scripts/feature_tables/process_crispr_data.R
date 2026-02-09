## Create Activity-only feature table for an ABC sample

# required packages and functions
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# load inputs
df <- fread(snakemake@input$crispr_features)
genes <- fread(snakemake@params$genes)
colnames(genes) <- c('chr', 'start', 'end', 'gene', 'score', 'strand', 'Ensembl_ID', 'gene_type')

# rename some column names for later utility
df = dplyr::rename(df, 'chr'='chrom', 'start'='chromStart', 'end'='chromEnd', 'TargetGene'='measuredGeneSymbol')

# drop elements where Regulated=NA
df = dplyr::filter(df, Regulated!="NA")

# filter to elements with target gene in gene universe (should be redundant with applying filter_genes=tss_universe in feature overlap)
df = dplyr::filter(df, TargetGene %in% genes$gene)

# save output to file
fwrite(df, file = snakemake@output$processed, sep = "\t", na = "NA", quote = FALSE)
