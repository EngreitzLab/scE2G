checkpoint features_required:
	input:
		feature_table_file = os.path.join(RESULTS_DIR, "{sample}", "feature_table.tsv")
	output:
		to_generate = os.path.join(RESULTS_DIR, "{sample}", "to_generate.txt") # file with "Kendall" "ARC" or "Neither"
	run:
		Kendall = False
		ARC = False
		with open(input.feature_table_file, "r") as f:
			for line in f:
				columns = line.strip().split("\t")
				if (("Kendall" in columns[1]) or ("Kendall" in columns[2])):
					Kendall = True
				if (("ARC.E2G.Score" in columns[1]) or ("ARC.E2G.Score" in columns[2])):
					ARC = True
		final_val = "Neither"
		if ARC or Kendall:
			final_val = "ARC"
		with open(output.to_generate, "w") as out:
			out.write(final_val)


# return file paths for features to generate
def features_to_generate(wildcards):
	with checkpoints.features_required.get(sample=wildcards.sample).output.to_generate.open() as f:
		val = f.read().strip()
		if val == "Kendall":
			return os.path.join(RESULTS_DIR, "{sample}", "Kendall", "Pairs.Kendall.tsv.gz")
		elif val == "ARC":
			return os.path.join(RESULTS_DIR, "{sample}", "ARC", "EnhancerPredictionsAllPutative_ARC.tsv.gz")
		else:
			return RESULTS_DIR

# activate generation of Kendall/ARC and format external_features_config
rule make_external_features_config:
	input:
		feature_inputs = features_to_generate
	params:
		dataset_config = config["cell_clusters"],
		e2g_path = config["encode_re2g_dir"]
	output:
		external_features_config = os.path.join(RESULTS_DIR, "{sample}", "external_features_config.tsv")
	conda:
		"../envs/sc_e2g.yml"
	resources:
		mem_mb=8*1000
	script:
		"../scripts/feature_computation/format_external_features_config_sc.R"
