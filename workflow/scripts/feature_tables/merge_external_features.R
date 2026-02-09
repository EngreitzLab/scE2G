## Merge features with ABC table

# save.image("merge_abc.rda")
# stop()

# required packages
library(data.table)
library(tidyverse)
library(GenomicRanges)

## Define functions --------------------------------------------------------------------------------

# function to overlap features to abc table
overlap_feature_with_abc <- function(abc, feature, feature_score_cols, agg_fun, fill_value) {
  
  # convert fill value to numeric
  fill_value <- suppressWarnings(as.numeric(fill_value))
  
  # only retain relevant columns from feature
  feature <- select(feature, chr, start, end, TargetGene, all_of(feature_score_cols))
  
  # create GRanges for abc  and feature E-G pairs
  abc_gr <- with(abc, GRanges(seqnames = paste0(chr, ":", TargetGene),
                              ranges = IRanges(start, end)))
  
  feat_gr <- with(feature, GRanges(seqnames = paste0(chr, ":", TargetGene),
                                   ranges = IRanges(start, end)))
  
  # set same seqlevels for both GRanges objects to avoid warnings
  seqlevels_all_pairs <- as.character(unique(c(seqnames(abc_gr), seqnames(feat_gr))))
  seqlevels(abc_gr) <- seqlevels_all_pairs
  seqlevels(feat_gr) <- seqlevels_all_pairs
  
  # find overlaps between abc and features
  ovl <- findOverlaps(abc_gr, feat_gr)
  
  # merge abc table with features
  merged <- cbind(abc[queryHits(ovl)], feature[subjectHits(ovl), feature_score_cols, with = FALSE])
  
  # aggregate feature scores, where multiple features overlapped one abc enhancer
  agg_cols <- setdiff(colnames(merged), feature_score_cols)
  merged <- aggregate_features(merged, feature_score_cols = feature_score_cols,
                               agg_cols = agg_cols, agg_fun = agg_fun)
  
  # report how many abc pairs overlap >= 1 features
  perc_overlap <- nrow(merged) / nrow(abc)
  message("Features overlapping ABC: ", round(perc_overlap * 100, digits = 2), "%")
  
  # get pairs from abc table that are missing from features 
  missing <- abc[setdiff(seq_len(nrow(abc)), queryHits(ovl)), ]
  
  if (nrow(missing)>0){
    # fill in missing value for features in missing pairs
    for (i in feature_score_cols) {
        missing[[i]] <- fill_value
    }
     # combine merged and missing pairs to create output
    output <- rbind(merged, missing)
  } else {
    output = merged
  }

  # sort output by cre position
  output <- output[order(chr, start, end, TargetGene), ]
  
  return(output)
  
}

# function to aggregate multiple overlaps
aggregate_features <- function(merged, feature_score_cols, agg_cols, agg_fun) {
  
  # aggregate feature columns as specified by agg_fun
  agg_list <- mapply(FUN = function(feat_col, agg_func, agg_cols) {
    agg_func <- get(agg_func)
    merged[, setNames(.(agg_func(get(feat_col))), feat_col), by = agg_cols]
  }, feat_col = feature_score_cols, agg_func = agg_fun, MoreArgs = list(agg_cols = agg_cols),
  SIMPLIFY = FALSE)
  
  # merge all the aggregates together to make collapsed data.frame
  output <- as.data.table(Reduce(function(df1, df2) merge(df1, df2, by = agg_cols), agg_list))
  
  return(output)
  
}

# function to merge features by target gene column
merge_feature_by_gene <- function(abc, feature, feature_score_cols){
    feature = dplyr::select(feature, TargetGene, all_of(feature_score_cols))
    abc = dplyr::left_join(abc, feature, by="TargetGene")

    n_genes_abc = length(unique(abc$TargetGene))
    n_represented = length(intersect(unique(abc$TargetGene), unique(feature$TargetGene)))
    pct_overlap <-n_represented/n_genes_abc
    message("Genes overlapping ABC: ", round(pct_overlap * 100, digits = 2), "%")

    return(abc)
}

## Process features --------------------------------------------------------------------------------

# inputs from snakemake
abc <- fread(snakemake@input$predictions_extended)
features <- fread(snakemake@input$feature_table_file, header = TRUE)
ext_ft = fread(snakemake@input$external_features_config, header=TRUE)

# filter to relevant data....
input_features = c(features$input_col, features$second_input)
input_features = unique(na.omit(input_features)) # remove NAs, keep unique
needed = setdiff(input_features, colnames(abc))
ext_ft = dplyr::filter(ext_ft, input_col %in% needed)
unique_sources = na.omit(unique(ext_ft$source_file))

# currently requires overlap of genomic coordinates
if (length(unique_sources>0)){
    for (i in 1:length(unique_sources)){
        message("Processing: ", unique_sources[i])
        this_ext = dplyr::filter(ext_ft, source_file == unique_sources[i])
        this_source = fread(unique_sources[i])
        # rename features to input_col
        colnames(this_source)[match(this_ext$source_col, colnames(this_source))] = this_ext$input_col
        
        if ("overlap" %in% this_ext$join_by){ # join by overlapping genomic coordinates of enhancer and matching target gene
            abc = overlap_feature_with_abc(abc, feature=this_source, feature_score_cols = this_ext$input_col, agg_fun = this_ext$aggregate_function, fill_value = NA_real_)
        } else if (all(this_ext$join_by == "TargetGene")) { # features just matched to target gene
        abc = merge_feature_by_gene(abc, feature=this_source, feature_score_cols = this_ext$input_col)
        } else {
            message("Features must be added by overlap or TargetGene.")
        }
    }
}

# write output to file
fwrite(abc, file = snakemake@output$plus_external_features, sep = "\t")
