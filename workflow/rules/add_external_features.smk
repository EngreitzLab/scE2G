# Feature requirement detection via upfront configuration evaluation
# (Replaces checkpoint features_required)

def features_to_generate(wildcards):
	"""Return ARC file path if needed, otherwise RESULTS_DIR.

	Uses upfront-computed BIOSAMPLE_NEEDS_ARC dict instead of checkpoint.
	"""
	if BIOSAMPLE_NEEDS_ARC.get(wildcards.sample, False):
		return os.path.join(RESULTS_DIR, wildcards.sample, "ARC", "EnhancerPredictionsAllPutative_ARC.tsv.gz")
	else:
		return RESULTS_DIR

# activate generation of Kendall/ARC and format external_features_config
rule make_external_features_config:
	input:
		feature_inputs = features_to_generate
	params:
		dataset_config = config["cell_clusters"],
		e2g_path = WORKFLOW_DIR
	output:
		external_features_config = os.path.join(RESULTS_DIR, "{sample}", "external_features_config.tsv")
	conda:
		"../envs/sc_e2g.yml"
	resources:
		mem_mb=8*1000
	script:
		"../scripts/feature_computation/format_external_features_config_sc.R"
