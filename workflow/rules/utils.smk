import filecmp
import tempfile


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

## Convert cell cluster config to biosample config for the ABC pipeline
## Each cell cluster is treated as a distinct biosample within the ABC framework
def make_biosample_config(cluster_config, biosample_config, results_dir):
	"""
	This function transforms a cell cluster configuration table into a biosample configuration table.
	It is designed for use in the ABC pipeline. In this context, each cell cluster is treated as an individual biosample.
	Additionally, the function sets up necessary file paths and directories for the biosample data.

	This function is called at parse time, so it runs on every single Snakemake
	invocation -- including `--dry-run`, `--unlock`, and each remote job's
	re-parse on a compute node. The biosample config it produces is a declared
	input of downstream rules (e.g. ENCODE_rE2G's make_biosample_feature_table),
	so rewriting it unconditionally would refresh its mtime on every invocation
	and make Snakemake consider the whole downstream DAG out of date, even when
	nothing changed. It is therefore written "if changed": the table is rendered
	to a temporary file and only moved into place when it differs from what is
	already there, leaving the mtime alone for no-op invocations while still
	updating (and correctly triggering reruns) whenever the cluster config
	really does change.

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

	# Render to a temporary file in the destination directory, so the bytes
	# compared below are produced by exactly the same code path that would
	# write the final file, and so the move into place is atomic.
	tmp_fd, tmp_path = tempfile.mkstemp(
		dir=output_dir,
		prefix=os.path.basename(biosample_config) + ".",
		suffix=".tmp"
	)
	os.close(tmp_fd)
	try:
		# mkstemp() is deliberately 0600; restore the permissions a plain
		# open() would have given the file so group-shared results stay readable.
		umask = os.umask(0)
		os.umask(umask)
		os.chmod(tmp_path, 0o666 & ~umask)

		df.to_csv(tmp_path, sep='\t', index=False)

		# Unchanged content: leave the existing file (and its mtime) untouched.
		if os.path.exists(biosample_config) and filecmp.cmp(tmp_path, biosample_config, shallow=False):
			return

		# First run, or the content really changed: publish atomically.
		os.replace(tmp_path, biosample_config)
		tmp_path = None
	finally:
		if tmp_path is not None and os.path.exists(tmp_path):
			os.remove(tmp_path)


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
	e2g_config["results_dir"] = config["results_dir"]
	e2g_config["model_dir"] = config["model_dir"]
	e2g_config["final_score_col"] = config["final_score_col"]

	# If files specified in scE2G, update ENCODE_rE2G
	if "gene_TSS500" in config:
		e2g_config["gene_TSS500"] = config["gene_TSS500"]
	if "genes" in config:
		e2g_config["genes"] = config["genes"]
	if "crispr_dataset" in config:
		e2g_config["crispr_dataset"] = config["crispr_dataset"]
	if "chr_sizes" in config:
		e2g_config["chr_sizes"] = config["chr_sizes"]
	if "regions_blocklist" in config:
		e2g_config["regions_blocklist"] = config["regions_blocklist"]
	if "macs2_genomesize" in config:
		e2g_config["macs2_genomesize"] = config["macs2_genomesize"]
	
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
		


