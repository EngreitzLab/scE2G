## Overlap an E-G features table with CRISPR data

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(GenomicRanges)
  source(file.path(snakemake@scriptdir, "get_fill_values.R"))
})

## Define functions --------------------------------------------------------------------------------

# function to merge features with crispr data
merge_feature_to_crispr <- function(crispr, features, feature_score_cols, agg_fun, fill_value) {
  # create GRanges for crispr and features E-G pairs
  crispr_gr <- with(crispr, GRanges(
    seqnames = paste0(chrom, ":", measuredGeneSymbol),
    ranges = IRanges(chromStart, chromEnd)
  ))

  feat_gr <- with(features, GRanges(
    seqnames = paste0(chr, ":", TargetGene),
    ranges = IRanges(start, end)
  ))

  # set same seqlevels for both GRanges objects to avoid warnings
  seqlevels_all_pairs <- as.character(unique(c(seqnames(crispr_gr), seqnames(feat_gr))))
  seqlevels(crispr_gr) <- seqlevels_all_pairs
  seqlevels(feat_gr) <- seqlevels_all_pairs

  # find overlaps between crispr and features
  ovl <- findOverlaps(crispr_gr, feat_gr)

  # merge rows of the crispr and features table based on overlapping E-G pairs
  merged <- cbind(crispr[queryHits(ovl)], features[subjectHits(ovl), feature_score_cols, with = FALSE])

  # aggregate feature scores, where multiple features overlapped one crispr enhancer
  agg_cols <- setdiff(colnames(merged), feature_score_cols)
  merged <- aggregate_features(merged,
    feature_score_cols = feature_score_cols,
    agg_cols = agg_cols, agg_fun = agg_fun[feature_score_cols]
  )

  # report how many crispr pairs overlap >= 1 features
  perc_overlap <- nrow(merged) / nrow(crispr)
  message("Features overlapping crispr: ", round(perc_overlap * 100, digits = 2), "%")

  # get pairs from crispr table that are missing from features
  missing <- crispr[setdiff(seq_len(nrow(crispr)), queryHits(ovl)), ]

  # fill in missing value for features in missing pairs
  for (i in feature_score_cols) {
    missing[[i]] <- fill_value
  }

  # only return crispr pairs that overlap features
  output <- merged

  # sort output by cre position
  output <- output[order(dataset, chrom, chromStart, chromEnd, measuredGeneSymbol), ]
  return(list(merged=output, missing=missing))
}

# function to aggregate multiple overlaps
aggregate_features <- function(merged, feature_score_cols, agg_cols, agg_fun) {
  # aggregate feature columns as specified by agg_fun
  agg_list <- mapply(
    FUN = function(feat_col, agg_func, agg_cols) {
      agg_func <- get(agg_func)
      merged[, setNames(.(agg_func(get(feat_col))), feat_col), by = agg_cols]
    }, feat_col = feature_score_cols, agg_func = agg_fun, MoreArgs = list(agg_cols = agg_cols),
    SIMPLIFY = FALSE
  )

  # merge all the aggregates together to make collapsed data.frame
  output <- as.data.table(Reduce(function(df1, df2) merge(df1, df2, by = agg_cols), agg_list))
}

## Overlap features with CRISPR data ---------------------------------------------------------------

# load feature table
features <- fread(snakemake@input$features)

# load crispri data and only retain relevant columns
crispr <- fread(snakemake@input$crispr)
crispr <- select(crispr, -c(pair_uid, merged_uid, merged_start, merged_end))

# load feature config file and only retain entries for features in input data
config <- fread(snakemake@input$feature_table_file)
config <- filter(config, feature %in% colnames(features))

# load tss annotations
tss <- fread(snakemake@input$tss, col.names = c("chr", "start", "end", "name", "score", "strand", "Ensembl_ID", "gene_type"))

# create vector with aggregation functions for each feature
agg_funs <- deframe(distinct(select(config, feature, aggregate_function)))

# filter out any CRISPR E-G pairs involving genes not part of the TSS universe or feature genes
filter_genes = "tss_universe"
if (filter_genes== "tss_universe") { 
  missing_genes <- setdiff(crispr$measuredGeneSymbol, tss$name)
  message("Removing CRISPR data for ", length(missing_genes), " genes not part of TSS universe")
  crispr <- filter(crispr, !measuredGeneSymbol %in% missing_genes)
} else if (filter_genes == "feature_genes") {
  missing_genes <- setdiff(crispr$measuredGeneSymbol, features$TargetGene)
  message("Removing CRISPR data for ", length(missing_genes), " genes not part of feature genes")
  crispr <- filter(crispr, !measuredGeneSymbol %in% missing_genes)
}

# overlap feature table with CRISPR data
df_list <- merge_feature_to_crispr(crispr,
  feature = features,
  feature_score_cols = unique(config$feature),
  agg_fun = agg_funs, fill_value = NA_real_
)
output = df_list$merged
missing = df_list$missing

# add TSS coordinates from TSS annotations
output <- tss %>%
  select(measuredGeneSymbol = name, startTSS_ref = start, endTSS_ref = end) %>%
  left_join(output, ., by = "measuredGeneSymbol")
missing <- tss %>%
  select(measuredGeneSymbol = name, startTSS_ref = start, endTSS_ref = end) %>%
  left_join(missing, ., by = "measuredGeneSymbol")

# re-calculate distance to tss based on crispr data for missing values
if ("distanceToTSS" %in% colnames(output)) {
  output <- output %>%
    mutate(enh_center = (chromStart + chromEnd) / 2) %>%
    mutate(distanceToTSS = case_when(
      is.na(distanceToTSS) & is.na(startTSS_ref) ~ abs(enh_center - (startTSS + endTSS) / 2),
      is.na(distanceToTSS) ~ abs(enh_center - (startTSS_ref + endTSS_ref) / 2),
      TRUE ~ distanceToTSS
    )) %>%
    select(-c(enh_center))
}

# remove columns added to compute missing distance to TSS
output <- select(output, -c(startTSS_ref, endTSS_ref))

# replace NAs with fill values if specified
if (snakemake@wildcards$nafill == "NAfilled") {
  # get fill values for each feature
  fill_values <- get_fill_values(features, config = config)

  # replace NAs with fill values
  output <- replace_na(output, replace = fill_values)
}

# write output to file
fwrite(output, file = snakemake@output$features, sep = "\t")
fwrite(missing, file = snakemake@output$missing, sep= "\t" )

