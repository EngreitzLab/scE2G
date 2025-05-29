
# merge features with crispr data
rule overlap_features_crispr_apply:
	input:
		prediction_file = os.path.join(RESULTS_DIR, "{cluster}", "{model_name}", "encode_e2g_predictions.tsv.gz"),
		crispr = config['crispr_dataset'],
		feature_table_file = os.path.join(RESULTS_DIR, "{cluster}", "feature_table.tsv"),
		tss = config['gene_TSS500']
	params:
		fill_value_script = os.path.join(SCRIPTS_DIR, "model_application", "get_fill_values.R")
	output: 
		features = os.path.join(RESULTS_DIR, "{cluster}", "{model_name}", "EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_{nafill}.tsv.gz"),
	conda:
		"../envs/sc_e2g.yml" 
	resources:
		mem_mb=32*1000
	script:
		"../scripts/model_application/merge_features_with_crispr_data_apply.R"

# calculate performance metrics 
rule crispr_benchmarking:
	input:
		crispr_features = expand(os.path.join(RESULTS_DIR, "{cluster}", "{model_name}", "EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz"), zip, cluster=BIOSAMPLE_DF["cluster"], model_name=BIOSAMPLE_DF["model_dir_base"]),
	output:
		comp_table = os.path.join(RESULTS_DIR, "crispr_benchmarking_performance_summary.tsv")
	params:
		model_names = BIOSAMPLE_DF["model_dir_base"].tolist(),
		model_thresh = BIOSAMPLE_DF["model_threshold"].tolist(),
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/sc_e2g.yml" 
	resources:
		mem_mb=64*1000
	shell: 
		""" 
		python {params.scripts_dir}/prediction_qc/benchmark_performance.py \
			--crispr_features "{input.crispr_features}" \
			--output_file {output.comp_table} \
			--model_thresholds "{params.model_thresh}" \
			--model_names "{params.model_names}" \
		"""

rule run_e2g_qnorm:
	input:
		final_features = os.path.join(RESULTS_DIR, "{cluster}", "genomewide_features.tsv.gz"),
	params:
		epsilon = 0.01,
		feature_table_file = lambda wildcards: encode_e2g.get_feature_table_file(wildcards.cluster, wildcards.model_name),
		trained_model = lambda wildcards: encode_e2g.get_trained_model(wildcards.cluster, wildcards.model_name),
		model_dir = lambda wildcards: encode_e2g._get_model_dir_from_wildcards(wildcards.cluster, wildcards.model_name, BIOSAMPLE_DF),
		tpm_threshold =  lambda wildcards: encode_e2g.get_tpm_threshold(wildcards.cluster, wildcards.model_name, BIOSAMPLE_DF),
		crispr_benchmarking = config["benchmark_performance"],
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/sc_e2g.yml"
	resources:
		mem_mb=encode_e2g.ABC.determine_mem_mb
	output: 
		prediction_file = os.path.join(RESULTS_DIR, "{cluster}", "{model_name}", "encode_e2g_predictions.tsv.gz")
	shell: 
		""" 
		python {params.scripts_dir}/model_application/run_e2g_cv.py \
			--predictions {input.final_features} \
			--feature_table_file {params.feature_table_file} \
			--epsilon {params.epsilon} \
			--tpm_threshold {params.tpm_threshold} \
			--trained_model {params.trained_model} \
			--model_dir {params.model_dir} \
			--crispr_benchmarking {params.crispr_benchmarking} \
			--output_file {output.prediction_file}
		"""
def get_gex_file(wildcards):
	with checkpoints.features_required.get(sample=wildcards.cluster).output.to_generate.open() as f:
		val = f.read().strip()
		if val == "Kendall" or val == "ARC":
			return os.path.join(RESULTS_DIR, wildcards.cluster, "Kendall", "gene_expression_metrics.tsv.gz")
		else:
			return RESULTS_DIR

rule element_and_gene_summaries:
	input:
		abc_gene_list = os.path.join(RESULTS_DIR, "{cluster}", "Neighborhoods", "GeneList.txt"),
		abc_element_list = os.path.join(RESULTS_DIR, "{cluster}", "Neighborhoods", "EnhancerList.txt"),
		prediction_file = os.path.join(RESULTS_DIR, "{cluster}", "{model_name}", "encode_e2g_predictions.tsv.gz"),
		gene_expr_file = get_gex_file
	params:
		tpm_threshold = lambda wildcards: encode_e2g.get_tpm_threshold(wildcards.cluster, wildcards.model_name, BIOSAMPLE_DF)
	conda:
		"../envs/sc_e2g.yml"
	resources:
		mem_mb=encode_e2g.ABC.determine_mem_mb
	output: 
		gene_list = os.path.join(RESULTS_DIR, "{cluster}", "{model_name}", "scE2G_gene_list.tsv.gz"),
		element_list = os.path.join(RESULTS_DIR, "{cluster}", "{model_name}", "scE2G_element_list.tsv.gz")
	script:
		"../scripts/prediction_qc/generate_element_gene_lists.R"

def get_count_file(wildcards, metric):
	with checkpoints.features_required.get(sample=wildcards.cluster).output.to_generate.open() as f:
		val = f.read().strip()
		if val == "Kendall" or val == "ARC":
			return os.path.join(RESULTS_DIR, wildcards.cluster, f"{metric}.txt")
		else:
			return RESULTS_DIR

rule get_stats_per_model_per_cluster:
	input:
		pred_full = os.path.join(RESULTS_DIR, "{cluster}", "{model_name}", "encode_e2g_predictions.tsv.gz"),
		pred_thresholded = os.path.join(RESULTS_DIR, "{cluster}", "{model_name}", "encode_e2g_predictions_threshold{threshold}.tsv.gz"),
		fragment_count = os.path.join(RESULTS_DIR, "{cluster}", "fragment_count.txt"),
		cell_count = lambda wildcards: get_count_file(wildcards, "cell_count"),
		umi_count = lambda wildcards: get_count_file(wildcards, "umi_count")
	params:
		score_column = "E2G.Score.qnorm"
	conda:
		"../envs/sc_e2g.yml"
	resources:
		mem_mb=32000
	output: 
		stats = os.path.join(RESULTS_DIR, "{cluster}", "{model_name}", "encode_e2g_predictions_threshold{threshold}_stats.tsv")
	script:
		"../scripts/prediction_qc/get_stats_per_cluster.R"

biosample_model_threshold = list(zip(BIOSAMPLE_DF["biosample"], BIOSAMPLE_DF["model_dir_base"], BIOSAMPLE_DF["model_threshold"]))
rule plot_stats:
	input:
		stats_files = [os.path.join(RESULTS_DIR, biosample, model_name, f"encode_e2g_predictions_threshold{threshold}_stats.tsv") for biosample,model_name,threshold in biosample_model_threshold]
	params:
		score_column = "E2G.Score.qnorm"
	conda:
		"../envs/sc_e2g.yml"
	resources:
		mem_mb=encode_e2g.ABC.determine_mem_mb
	output: 
		all_stats = os.path.join(RESULTS_DIR, "qc_plots", "all_qc_stats.tsv"),
		ds_stats = os.path.join(RESULTS_DIR, "qc_plots", "dataset_metrics.pdf"),
		enh_stats = os.path.join(RESULTS_DIR, "qc_plots", "enhancer_metrics.pdf"),
		eg_stats = os.path.join(RESULTS_DIR, "qc_plots", "eg_metrics.pdf"),
		gene_stats  = os.path.join(RESULTS_DIR, "qc_plots", "gene_metrics.pdf"),
		dist_size_stats = os.path.join(RESULTS_DIR, "qc_plots", "distance_and_size_metrics.pdf"),
		all_distributions = os.path.join(RESULTS_DIR, "qc_plots", "qc_metric_distributions.pdf")
	script:
		"../scripts/prediction_qc/plot_all_qc_stats.R"

