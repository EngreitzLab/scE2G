## Create Activity-only feature table for an ABC sample

# required packages and functions
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  source(file.path(snakemake@scriptdir, "get_fill_values.R"))
})

# load feature config file
config <- fread(snakemake@input$feature_table_file)
# load ABC table
df <- fread(snakemake@input$plus_external_features)

# confirm that all required input feature columns are present
input_features = c(config$input_col, config$second_input)
input_features = unique(input_features)
input_features <- na.omit(input_features[nzchar(input_features)])

if (any(!(input_features %in% colnames(df)))) {
    missing = input_features[!(input_features %in% colnames(df))]
    message(missing)
    stop("Above required features are not present.")
}

# calculate interaction terms and name them correctly
intx <- dplyr::filter(config, !is.na(second_input), nzchar(second_input))
for (i in seq_len(nrow(intx))) {
  df[[intx$feature[i]]] <- df[[intx$input_col[i]]] * df[[intx$second_input[i]]]
}

# rename single features to final names
single <- dplyr::filter(config, !(feature %in% intx$feature))
for (i in seq_len(nrow(single))){
  names(df)[names(df) == single$input_col[i]] <- single$feature[i]
}

# fill NAs 
fill_values <- get_fill_values(df, config = config)
df <- replace_na(df, replace = fill_values)

# reorder columns for output
core_cols <- c(
  "chr", "start", "end", "name", "class", "TargetGene", "TargetGeneTSS",
  "TargetGeneIsExpressed", "TargetGeneEnsembl_ID", "isSelfPromoter",
  "CellType", "distance", "ABC.Score"
)
core_cols <- core_cols[!core_cols %in% input_features] # remove core_cols that were renamed
output <- select(df, all_of(core_cols), all_of(config$feature))

# save output to file
fwrite(df, file = snakemake@output[[1]], sep = "\t", na = "NA", quote = FALSE)
