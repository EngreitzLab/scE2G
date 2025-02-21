# required packages
library(data.table)
library(dplyr)

# inputs from snakemake
new_feature = snakemake@input[1]
dataset_config = fread(snakemake@params$dataset_config)
ds = snakemake@wildcards$sample
e2g_path = snakemake@params$e2g_path

# if specified, load external_feature_config and make paths absolute
datasets_this = dplyr::filter(dataset_config, cluster==ds)
if ("external_features_config" %in% colnames(datasets_this)){
    if (nchar(datasets_this$external_features_config[1])>5) {
        efc = fread(datasets_this$external_features_config[1])
        # make paths absolute
        for (i in 1:nrow(efc)){
            potential_file = file.path(e2g_path, efc$source_file[i])
            if (file.exists(potential_file)){
                efc$source_file[i] = potential_file
            }
        }
    }
}
# if it doesn't exist, initialize empty df
if (!(exists('efc') && !is.null(efc))){
    message("Initializing empty efc")
    efc =  data.frame(input_col = character(), source_col = character(), aggregate_function = character(), join_by = character(), source_file = character())
}
col_names = c("input_col", "source_col", "aggregate_function", "join_by", "source_file")

ARC_rows = data.frame(input_col =c("mean_log_normalized_rna", "RnaPseudobulkTPM", "RnaDetectedPercent", "Kendall", "ARC.E2G.Score"),
		source_col=c("mean_log_normalized_rna", "RnaPseudobulkTPM", "RnaDetectedPercent", "Kendall", "ARC.E2G.Score"),
		aggregate_function=c("mean", "mean", "mean", "max", "mean"))
ARC_rows$join_by="overlap"
ARC_rows$source_file = new_feature

Kendall_row = data.frame("Kendall", "Kendall", "max", "overlap", new_feature); names(Kendall_row) = col_names

if (grepl("Pairs.Kendall.tsv.gz", new_feature)){
    efc = rbind(efc, Kendall_row)
}
if (grepl("EnhancerPredictionsAllPutative_ARC.tsv.gz", new_feature)){
    efc = rbind(efc, ARC_rows)
}

fwrite(efc, file = snakemake@output$external_features_config, sep = "\t")
