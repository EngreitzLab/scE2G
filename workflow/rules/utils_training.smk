## Import training configuration for ENCODE_rE2G
def get_e2g_training_config(config, encode_re2g_dir):
	"""
		This function reads the ENCODE_rE2G configuration file, updates certain paths with values from 
		the provided `config`, and returns the updated configuration.

		:param config: Configuration of sc-E2G.
		:param encode_re2g_dir: Path to the ENCODE_rE2G directory.
	"""

	e2g_config_file = os.path.join(encode_re2g_dir, "config/config_training.yaml")
	with open(e2g_config_file, 'r') as stream:
		e2g_config = yaml.safe_load(stream)

	# Update ENCODE_rE2G configuration
	e2g_config["E2G_DIR_PATH"] = encode_re2g_dir
	e2g_config["dataset_config"] = config["ABC_BIOSAMPLES"]
	e2g_config["results_dir"] = config["results_dir"]
	e2g_config["model_config"] = config["model_config"]
	e2g_config["run_feature_analysis"] = config["run_feature_analysis"]
	e2g_config["reference_features"] = config["reference_features"]
	e2g_config["final_score_col"] = config["final_score_col"]

	# If files specified in scE2G, update ENCODE_rE2G config
	if "gene_TSS500" in config:
		e2g_config["gene_TSS500"] = config["gene_TSS500"]
	if "genes" in config:
		e2g_config["genes"] = config["genes"]
	if "crispr_dataset" in config:
		e2g_config["crispr_dataset"] = config["crispr_dataset"]

	return e2g_config


