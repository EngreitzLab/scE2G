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
	# Memory resource calculator for snakemake rules
	input_size_mb = input.size_mb
	if ".gz" in str(input):
		input_size_mb *= 8  # assume gz compressesed the file <= 8x
	attempt_multiplier = 2 ** (attempt - 1)  # Double memory for each retry
	mem_to_use_mb = attempt_multiplier *  max(4 * input_size_mb, min_gb * 1000)
	return min(mem_to_use_mb, MAX_MEM_MB)

#cluster	rna_matrix_file	atac_frag_file	model_dir	HiC_file	HiC_type	HiC_resolution	alt_TSS	alt_genes
def get_cell_cluster_config(cell_type_config, model_dir, n_iter, cell_cluster_config, scratch_dir):
	if os.path.isfile(cell_cluster_config):
		res = pd.read_table(cell_cluster_config)
	else:
		# column values
		L2 = cell_type_config['L2'].unique().tolist()
		L3 = cell_type_config['L3'].unique().tolist()
		iter_values = list(range(1, n_iter + 1))
		type_values = ["real", "permute"]
		level_values = ["L1", "L2", "L3"]

		rows = []

		for i in iter_values:
			#rows.append((i, "real", "L1", "all"))
			#rows.extend((i, t, "L2", identity) for t, identity in itertools.product(type_values, L2))
			rows.extend((i, "permute", "L3_from_L1", identity) for identity in L3)
			rows.extend((i, "permute", "L3_from_L2", identity) for identity in L3)
			rows.extend((i, "real", "L3", identity) for identity in L3)

		# convert to df and populate necessary columns
		res = pd.DataFrame(rows, columns=["iter", "type", "level", "identity"])
		res["cluster"] = res.apply(lambda row: f"iter{row['iter']}.{row['type']}.{row['level']}.{row['identity']}", axis=1)
		res["atac_frag_file"] = res.apply(lambda row: os.path.join(scratch_dir, f"iter{row['iter']}.{row['type']}.{row['level']}.{row['identity']}", "atac_fragments.tsv.gz"), axis=1)
		res["rna_matrix_file"] = res.apply(lambda row: os.path.join(scratch_dir, f"iter{row['iter']}.{row['type']}.{row['level']}.{row['identity']}", "rna_matrix.csv.gz"), axis=1)
		res["model_dir"] = model_dir
		res[["HiC_file", "HiC_type", "HiC_resolution", "alt_TSS", "alt_genes"]] = ''

		# make output directory if necesary & save
		output_dir = os.path.dirname(cell_cluster_config)
		if not os.path.exists(output_dir):
			os.makedirs(output_dir)
		res.to_csv(cell_cluster_config, sep='\t', index=False)

	return res


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
	if os.path.isfile(biosample_config):
		df = pd.read_table(cluster_config)
	else:
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


## Import configuration for ENCODE_rE2G
def get_e2g_config(config, encode_re2g_dir):
	"""
		This function reads the ENCODE_rE2G configuration file, updates certain paths with values from 
		the provided `config`, and returns the updated configuration.

		:param config: Configuration of scE2G.
		:param encode_re2g_dir: Path to the ENCODE_rE2G directory.
	"""

	e2g_config_file = os.path.join(encode_re2g_dir, "config/config.yaml")
	with open(e2g_config_file, 'r') as stream:
		e2g_config = yaml.safe_load(stream)

	# Update ENCODE_rE2G configuration
	e2g_config["E2G_DIR_PATH"] = encode_re2g_dir
	e2g_config["ABC_BIOSAMPLES"] = config["ABC_BIOSAMPLES"]
	e2g_config["IGV_dir"] = IGV_DIR
	e2g_config["results_dir"] = RESULTS_DIR
	e2g_config["model_dir"] = config["model_dir"]
	e2g_config["final_score_col"] = config["final_score_col"]
	e2g_config["feature_table"] = FEATURE_TABLE

	# If files specified in scE2G, update ENCODE_rE2G
	if "gene_TSS500" in config:
		e2g_config["gene_TSS500"] = config["gene_TSS500"]
	if "genes" in config:
		e2g_config["genes"] = config["genes"]
	if "crispr_dataset" in config:
		e2g_config["crispr_dataset"] = config["crispr_dataset"]
	
	return e2g_config

# update scE2G config to have consistent gene reference files to E2G 
def update_scE2G_config(config, e2g_config, encode_re2g_dir):
	if "crispr_dataset" not in config:
		config["crispr_dataset"] = os.path.join(encode_re2g_dir, e2g_config["crispr_dataset"])
	
	if "gene_TSS500" not in config:
		config["gene_TSS500"] = os.path.join(encode_re2g_dir, e2g_config["gene_TSS500"])

	if "chr_sizes" not in config:
		config["chr_sizes"] = os.path.join(encode_re2g_dir, e2g_config["chr_sizes"])

	return config


# decide whether ARC-E2G should use "ABC.Score" or "powerlaw.Score"
def get_abc_score_col(cluster):
	row = CELL_CLUSTER_DF.loc[CELL_CLUSTER_DF["cluster"] == cluster].iloc[0]
	if pd.isnull(row["HiC_type"]):
		return "powerlaw.Score"
	else:
		return "ABC.Score"
		


