# required packages
library(data.table)
library(dplyr)

# inputs from snakemake
dataset_config = fread(snakemake@input$dataset_config)
ds = snakemake@wildcards$dataset
e2g_path = snakemake@params$e2g_path

# make external features config
datasets_this = dplyr::filter(dataset_config, biosample==ds)
if ("external_features_config" %in% colnames(datasets_this)){
    if (!(is.na(datasets_this$external_features_config[1]))) {
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
}
if (!(exists('efc') && !is.null(efc))){
    efc =  data.frame(input_col = character(), source_col = character(), aggregate_function = character(), join_by = character(), source_file = character())
}
fwrite(efc, file = snakemake@output$external_features_config, sep = "\t")
