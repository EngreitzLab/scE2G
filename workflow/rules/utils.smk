import glob
import numpy as np

MAX_MEM_MB = 250 * 1000  # 250GB

## Update paths in the config obj to absolute path
def make_paths_absolute(obj, base_path):
	"""
	Use absolute paths to be compatible with github submodules
	Recursively go through the dictionary and convert relative paths to absolute paths.
	"""
	if isinstance(obj, dict):
		for key, value in obj.items():
			obj[key] = make_paths_absolute(value, base_path)
	elif isinstance(obj, str):
		# We assume all strings are paths. If converting the string
		# to an absolute path results in a valid file, then the str was a path
		new_file = os.path.join(base_path, obj)
		if os.path.exists(new_file):
			return new_file
	return obj

def determine_mem_mb(wildcards, input, attempt, min_gb=8):
	"""Memory resource calculator for snakemake rules."""
	input_size_mb = input.size_mb
	if ".gz" in str(input):
		input_size_mb *= 8  # assume gz compressed the file <= 8x
	attempt_multiplier = 2 ** (attempt - 1)  # Double memory for each retry
	mem_to_use_mb = attempt_multiplier * max(4 * input_size_mb, min_gb * 1000)
	return min(mem_to_use_mb, MAX_MEM_MB)

## Convert cell cluster config to biosample config for the ABC pipeline
## Each cell cluster is treated as a distinct biosample within the ABC framework
def make_biosample_config(cluster_config, biosample_config, results_dir):
	"""
	This function transforms a cell cluster configuration table into a biosample configuration table.
	It is designed for use in the ABC pipeline. In this context, each cell cluster is treated as an individual biosample.
	Additionally, the function sets up necessary file paths and directories for the biosample data.

	:param cluster_config: Path to the input cell cluster configuration file.
	:param biosample_config: Path for the output biosample configuration file.
	:param results_dir: Directory where result files are stored.
	"""

	df = pd.read_csv(cluster_config, sep='\t')
	
	# Assign each cell cluster as a biosample and prepare respective data paths
	df['biosample'] = df['cluster']
	df['DHS'] = ''
	df['ATAC'] = df['cluster'].apply(
		lambda cluster: os.path.join(results_dir, cluster, "tagAlign", "tagAlign.sort.gz")
	)
	df['H3K27ac'] = ''
	df['default_accessibility_feature'] = 'ATAC'

	output_dir = os.path.dirname(biosample_config)

	# Create the output directory if it doesn't exist
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	df.to_csv(biosample_config, sep='\t', index=False)


## Import configuration for ABC submodule
def get_abc_config(config, abc_dir):
	"""
	This function reads the ABC configuration file, updates certain paths with values from
	the provided scE2G `config`, and returns the updated configuration.

	:param config: Configuration of scE2G.
	:param abc_dir: Path to the ABC directory.
	"""
	abc_config_file = os.path.join(abc_dir, "config/config.yaml")
	with open(abc_config_file, 'r') as stream:
		abc_config = yaml.safe_load(stream)

	abc_config["ABC_DIR_PATH"] = abc_dir
	abc_config["results_dir"] = config["results_dir"]

	if "dataset_config" in config:
		abc_config["biosamplesTable"] = config["dataset_config"]
	else:
		abc_config["biosamplesTable"] = config["ABC_BIOSAMPLES"]

	if "gene_TSS500" in config:
		abc_config["ref"]["genome_tss"] = config["gene_TSS500"]
	if "genes" in config:
		abc_config["ref"]["genes"] = config["genes"]
	if "chr_sizes" in config:
		abc_config["ref"]["chrom_sizes"] = config["chr_sizes"]
	if "regions_blocklist" in config:
		abc_config["ref"]["regions_blocklist"] = config["regions_blocklist"]
	if "macs2_genomesize" in config:
		abc_config["params_macs"]["genome_size"] = config["macs2_genomesize"]

	return abc_config

# update scE2G config to have consistent gene reference files from ABC
def update_scE2G_config(config, abc_config, abc_dir):
	if "crispr_dataset" not in config:
		# crispr_dataset is now in scE2G resources, not ABC
		pass

	if "gene_TSS500" not in config:
		# Make path absolute if it's relative to ABC directory
		genome_tss = abc_config["ref"]["genome_tss"]
		if not os.path.isabs(genome_tss):
			genome_tss = os.path.join(abc_dir, genome_tss)
		config["gene_TSS500"] = genome_tss

	if "chr_sizes" not in config:
		# Make path absolute if it's relative to ABC directory
		chrom_sizes = abc_config["ref"]["chrom_sizes"]
		if not os.path.isabs(chrom_sizes):
			chrom_sizes = os.path.join(abc_dir, chrom_sizes)
		config["chr_sizes"] = chrom_sizes

	return config

## Functions consolidated from ENCODE_rE2G

def make_accessibility_file_df(biosample_df, biosample_activities):
	"""Create a dataframe with one row per accessibility file per biosample."""
	df = biosample_df[["biosample", "ATAC", "DHS"]].copy()
	df["single_access_file"] = ""
	df["access_base"] = ""
	df["access_simple_id"] = ""

	new_rows = []
	for index, row in df.iterrows():
		counter = 1
		this_biosample = row["biosample"]
		default_accessibility = biosample_activities[this_biosample]
		for this_file in row[default_accessibility].split(","):
			new_row = row.copy()
			new_row["single_access_file"] = this_file
			new_row["access_base"] = os.path.splitext(os.path.basename(this_file))[0]
			new_row["access_simple_id"] = default_accessibility + f"_id{counter}"

			if (os.path.splitext(os.path.basename(this_file))[1]==".bam") or ("tagAlign.gz" in os.path.basename(this_file)):
				new_rows.append(new_row)
				counter += 1

	new_df = pd.DataFrame(new_rows)
	return new_df

def expand_biosample_df(biosample_df):
	"""Expand biosample dataframe with model directories and thresholds."""
	if "model_dir" not in biosample_df.columns:
		biosample_df['model_dir'] = np.nan
	biosample_df['model_dir_base'] = ''
	biosample_df['model_threshold'] = float(0)
	biosample_df["tpm_threshold"] = float(0)

	new_rows = []
	for index, row in biosample_df.iterrows():
		if pd.isna(row['model_dir']):
			biosample_df.loc[index, "model_dir"] = os.path.normpath(_get_biosample_model_dir_from_row(row))

		for model_dir in biosample_df.loc[index, 'model_dir'].split(','):
			new_row = row.copy()
			new_row['model_dir'] = model_dir
			new_row['model_dir_base'] = os.path.basename(model_dir)
			new_rows.append(new_row)
	new_df = pd.DataFrame(new_rows)
	new_df['model_threshold'] = [float(get_model_threshold(this_biosample, this_model, new_df)) for this_biosample, this_model in zip(new_df['biosample'], new_df['model_dir_base'])]
	new_df['tpm_threshold'] = [float(get_tpm_threshold(this_biosample, this_model, new_df)) for this_biosample, this_model in zip(new_df['biosample'], new_df['model_dir_base'])]

	return new_df

def _validate_model_dir(potential_dir):
	"""Validate that a model directory contains required files."""
	files = os.listdir(potential_dir)
	if "model.pkl" not in files:
		raise Exception("model.pkl is not provided in the specified model directory")
	if "feature_table.tsv" not in files:
		raise Exception("feature_table.tsv is not provided in the specified model directory")
	if sum(s.startswith("threshold_") for s in files) == 0:
		raise Exception("A threshold is not provided in the specified model directory")
	elif sum(s.startswith("threshold_") for s in files) > 1:
		raise Exception("More than one threshold is provided in the specified model directory")

def _get_model_dir_from_wildcards(biosample, model_name, biosample_df=None):
	"""Get model directory path from biosample and model name."""
	if biosample_df is None:
		df = BIOSAMPLE_DF
	else:
		df = biosample_df
	model_dir = df.loc[(df["biosample"]==biosample) & (df["model_dir_base"]==model_name), "model_dir"]
	return model_dir.values[0]

def _get_biosample_model_dir_from_row(row):
	"""Infer model directory from biosample row configuration."""
	if "model_dir" in BIOSAMPLE_DF.columns:
		if pd.notna(row["model_dir"]):
			_validate_model_dir(row["model_dir"])
			return row["model_dir"]

	access_type = row["default_accessibility_feature"].lower()
	input_params = [access_type]
	if not pd.isna(row["H3K27ac"]):
		input_params.append("h3k27ac")

	hic_file = row["HiC_file"]
	if pd.isna(hic_file):
		input_params.append("powerlaw")
	elif row["HiC_type"] == "avg":
		input_params.append("avg_hic")
	elif os.path.basename(hic_file) == os.path.basename(config["MEGAMAP_HIC_FILE"]):
		input_params.append("megamap")
	elif row["HiC_type"] == "hic":
		input_params.append("intact_hic")

	model_folder = os.path.join(MODEL_DIR, "_".join(input_params))
	if not os.path.exists(model_folder):
		raise Exception(f"{model_folder} not found. Model with input params not supported")
	return model_folder

def get_feature_table_file(biosample, model_name):
	"""Get path to feature table file for a biosample/model combination."""
	model_dir = _get_model_dir_from_wildcards(biosample, model_name)
	return os.path.join(model_dir, "feature_table.tsv")

def get_trained_model(biosample, model_name):
	"""Get path to trained model file for a biosample/model combination."""
	model_dir = _get_model_dir_from_wildcards(biosample, model_name)
	return os.path.join(model_dir, "model.pkl")

def get_model_threshold(biosample, model_name, biosample_df=None):
	"""Get model score threshold for a biosample/model combination."""
	model_dir = _get_model_dir_from_wildcards(biosample, model_name, biosample_df)
	threshold_files = glob.glob(os.path.join(model_dir, 'score_threshold_*'))
	assert len(threshold_files) == 1, "Should have exactly 1 threshold file in directory"
	threshold_file = os.path.basename(threshold_files[0])
	return threshold_file.split("threshold_")[1]

def get_tpm_threshold(biosample, model_name, biosample_df=None):
	"""Get TPM threshold for a biosample/model combination."""
	model_dir = _get_model_dir_from_wildcards(biosample, model_name, biosample_df)

	tpm_files = glob.glob(os.path.join(model_dir, 'tpm_threshold_*'))
	if len(tpm_files) == 0:
		return 0
	else:
		tpm_file = os.path.basename(tpm_files[0])
	return tpm_file.split("threshold_")[1]


# decide whether ARC-E2G should use "ABC.Score" or "powerlaw.Score"
def get_abc_score_col(cluster):
	row = CELL_CLUSTER_DF.loc[CELL_CLUSTER_DF["cluster"] == cluster].iloc[0]
	if pd.isnull(row["HiC_type"]):
		return "powerlaw.Score"
	else:
		return "ABC.Score"
		


