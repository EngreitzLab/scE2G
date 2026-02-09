## Create Activity-only feature table for an ABC sample

# required packages and functions
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# load feature config file
config <- fread(snakemake@input$feature_table_file)
# load ABC table and compute denominator
abc <- fread(snakemake@input$abc)
abc$ABC.Numerator = abc$ABC.Score.Numerator
abc$ABC.Denominator = abc$ABC.Score/abc$ABC.Numerator

# get a list of all input feature names from feature table
input_features = c(config$input_col, config$second_input)
input_features = unique(na.omit(input_features)) # remove NAs, keep unique

# core columns from ABC for output
core_cols <- c(
  "chr", "start", "end", "name", "class", "TargetGene", "TargetGeneTSS",
  "TargetGeneIsExpressed", "TargetGeneEnsembl_ID", "isSelfPromoter",
  "CellType", "distance"
)
# filter to existing core columns that are present in ABC or that will be added as features later
core_cols <- core_cols[core_cols %in% colnames(abc)]
core_cols <- core_cols[!core_cols %in% input_features]
output <- select(abc, all_of(core_cols))

# assemble feature table:

# add features present in ABC table (eg. activity_base)
output <- abc %>%
  select(name, TargetGene, any_of(input_features)) %>%
  left_join(output, ., by = c("name", "TargetGene"))

# add number of candidate enhancers between enhancers and genes if required
if (file_test("-f", snakemake@input$NumCandidateEnhGene)){ # this is a file and not a directory
	NumCandidateEnhGene <- fread(snakemake@input$NumCandidateEnhGene,
	header = TRUE,
	col.names = c("name", "TargetGene", "numCandidateEnhGene")
	)
	output <- NumCandidateEnhGene %>%
  		select(name, TargetGene, any_of(input_features)) %>%
  		left_join(output, ., by = c("name", "TargetGene"))
}

# # add number of protein-coding TSSs between enhancer and target if required
if (file_test("-f", snakemake@input$NumTSSEnhGene)){ # this is a file and not a directory
	NumTSSEnhGene <- fread(snakemake@input$NumTSSEnhGene,
	header = TRUE,
	col.names = c("name", "TargetGene", "numTSSEnhGene")
	)
	output <- NumTSSEnhGene %>%
	select(name, TargetGene, any_of(input_features)) %>%
	left_join(output, ., by = c("name", "TargetGene"))
}

# add number of enhancers wtihin 5kb of target if required
if (file_test("-f", snakemake@input$NumEnhancersEG5kb)){ # this is a file and not a directory
	NumEnhancersEG5kb <- fread(snakemake@input$NumEnhancersEG5kb,
	header = FALSE, sep = "\t",
	col.names = c("name", "numNearbyEnhancers")
	)
	output <- NumEnhancersEG5kb %>%
	select(name, any_of(input_features)) %>%
	left_join(output, ., by = "name")
}

# add sum of activity of enhancers wtihin 5kb of target if required
if (file_test("-f", snakemake@input$SumEnhancersEG5kb)){ # this is a file and not a directory
	SumEnhancersEG5kb <- fread(snakemake@input$SumEnhancersEG5kb,
	header = FALSE, sep = "\t",
	col.names = c("name", "sumNearbyEnhancers")
	)
	output <- SumEnhancersEG5kb %>%
	select(name, any_of(input_features)) %>%
	left_join(output, ., by = "name")
}

# add ubiquitously expressed and P2 gene info
ubiqExprGenes <- fread(snakemake@input$geneClass)

output <- ubiqExprGenes %>%
  select(TargetGene, any_of(input_features)) %>%
  left_join(output, ., by = "TargetGene")

# save output to file
fwrite(output, file = snakemake@output[[1]], sep = "\t", na = "NA", quote = FALSE)
