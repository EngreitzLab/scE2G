from functools import partial
		
checkpoint basic_features_required: # (vs scE2G feature_required)
	input:
		feature_table_file = os.path.join(RESULTS_DIR, "{sample}", "feature_table.tsv")
	output:
		numCandidateEnhGene = os.path.join(RESULTS_DIR, "{sample}", "new_features", "generate_numCandidateEnhGene.txt"), # file with NumCandidateEnhGene, NumTSSEnhGene, NumEnhancersEGXkb, SumEnhancersEGXkb
		numTSSEnhGene = os.path.join(RESULTS_DIR, "{sample}", "new_features", "generate_numTSSEnhGene.txt"),
		nearbyEnhancers = os.path.join(RESULTS_DIR, "{sample}", "new_features", "generate_nearbyEnhancers.txt"),
	run:
		req_numCandidateEnhGene = "False"
		req_numTSSEnhGene = "False"
		req_nearbyEnhancers = "False"
		with open(input.feature_table_file, "r") as f:
			for line in f:
				columns = line.strip().split("\t")
				if (("numCandidateEnhGene" in columns[1]) or ("numCandidateEnhGene" in columns[2])):
					req_numCandidateEnhGene = "True"
				if (("numTSSEnhGene" in columns[1]) or ("numTSSEnhGene" in columns[2])):
					req_numTSSEnhGene = "True"
				if (("numNearbyEnhancers" in columns[1])or  ("numNearbyEnhancers" in columns[2])):
					req_nearbyEnhancers = "True"
				if (("sumNearbyEnhancers" in columns[1])or  ("sumNearbyEnhancers" in columns[2])):
					req_nearbyEnhancers = "True"

		with open(output.numCandidateEnhGene, "w") as out:
			out.write(req_numCandidateEnhGene)
		with open(output.numTSSEnhGene, "w") as out:
			out.write(req_numTSSEnhGene)
		with open(output.nearbyEnhancers, "w") as out:
			out.write(req_nearbyEnhancers) 


# return file paths for features to generate
def get_numCandidateEnhGene_file(wildcards):
	with checkpoints.basic_features_required.get(sample=wildcards.biosample).output.numCandidateEnhGene.open() as f:
		val = f.read().strip()
		if val == "True":
			return os.path.join(RESULTS_DIR, wildcards.biosample, "new_features", "NumCandidateEnhGene.tsv")
		else:
			return RESULTS_DIR

def get_numTSSEnhGene_file(wildcards):
	with checkpoints.basic_features_required.get(sample=wildcards.biosample).output.numTSSEnhGene.open() as f:
		val = f.read().strip()
		if val == "True":
			return os.path.join(RESULTS_DIR, wildcards.biosample, "new_features", "NumTSSEnhGene.tsv")
		else:
			return RESULTS_DIR

def get_numNearbyEnhancers_file(wildcards):
	with checkpoints.basic_features_required.get(sample=wildcards.biosample).output.nearbyEnhancers.open() as f:
		val = f.read().strip()
		if val == "True":
			return os.path.join(RESULTS_DIR, wildcards.biosample, "new_features", "NumEnhancersEG5kb.txt")
		else:
			return RESULTS_DIR

def get_sumNearbyEnhancers_file(wildcards):
	with checkpoints.basic_features_required.get(sample=wildcards.biosample).output.nearbyEnhancers.open() as f:
		val = f.read().strip()
		if val == "True":
			return os.path.join(RESULTS_DIR, wildcards.biosample, "new_features", "SumEnhancersEG5kb.txt")
		else:
			return RESULTS_DIR

# generate feature "numCandidateEnhGene"
rule generate_num_candidate_enh_gene:
	input:
		abc_predictions = lambda wildcards: os.path.join(ABC_BIOSAMPLES_DIR[wildcards.biosample], "Predictions", "EnhancerPredictionsAllPutative.tsv.gz"),
	params:
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=partial(determine_mem_mb, min_gb=16)
	output:
		NumCandidateEnhGene = os.path.join(RESULTS_DIR, "{biosample}", "new_features", "NumCandidateEnhGene.tsv")
	shell: 
		""" 
		python {params.scripts_dir}/feature_tables/gen_num_candidate_enh_gene.py \
			--abc_predictions {input.abc_predictions} \
			--out_file {output.NumCandidateEnhGene}
		"""

# generate feature "numTSSEnhGene"
rule generate_num_tss_enh_gene:
	input:
		abc_predictions = lambda wildcards: os.path.join(ABC_BIOSAMPLES_DIR[wildcards.biosample], "Predictions", "EnhancerPredictionsAllPutative.tsv.gz"),
	params:
		gene_TSS500 = config['gene_TSS500'],
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=partial(determine_mem_mb, min_gb=32)
	output:
		numTSSEnhGene = os.path.join(RESULTS_DIR, "{biosample}", "new_features", "NumTSSEnhGene.tsv"),
		extendedEnhancerRegions = temp(os.path.join(RESULTS_DIR, "{biosample}",  "new_features", "extendedEnhancerRegions.txt")),
		enhancerTSSInt = temp(os.path.join(RESULTS_DIR, "{biosample}", "new_features", "extendedEnhancerRegions_TSS_int.tsv.gz"))
	shell: 
		""" 
		python {params.scripts_dir}/feature_tables/gen_num_tss_enh_gene.py \
			--abc_predictions {input.abc_predictions} \
			--ref_gene_tss {params.gene_TSS500} \
			--extended_enhancers {output.extendedEnhancerRegions} \
			--enhancer_tss_int {output.enhancerTSSInt} \
			--out_file {output.numTSSEnhGene}
		"""

# generate features "numNearbyEnhancers" and "sumNearbyEnhancers"
rule generate_num_sum_enhancers:
	input:
		abc_predictions = lambda wildcards: os.path.join(ABC_BIOSAMPLES_DIR[wildcards.biosample], "Predictions", "EnhancerPredictionsAllPutative.tsv.gz"),
		enhancer_list = lambda wildcards: os.path.join(ABC_BIOSAMPLES_DIR[wildcards.biosample], "Neighborhoods", "EnhancerList.txt"),
	params:
		chr_sizes = config['chr_sizes'],
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=partial(determine_mem_mb, min_gb=16)
	output: 
		NumEnhancersEG = os.path.join(RESULTS_DIR, "{biosample}", "new_features", "NumEnhancersEG{kb}kb.txt"),
		SumEnhancersEG = os.path.join(RESULTS_DIR, "{biosample}", "new_features", "SumEnhancersEG{kb}kb.txt"),
		enhMidpoint = temp(os.path.join(RESULTS_DIR, "{biosample}", "new_features", "enhancerMidpoint_{kb}kb.txt")),
		enhExpanded = temp(os.path.join(RESULTS_DIR, "{biosample}", "new_features", "enhancerMidpoint_exp{kb}kb.txt")),
		predSlim = temp(os.path.join(RESULTS_DIR, "{biosample}", "new_features", "EnhancerPredictionsAllPutative_{kb}kb.slim.txt")),
		enhPredInt = temp(os.path.join(RESULTS_DIR, "{biosample}", "new_features",  "enhancerExp{kb}kb_intPred.txt")),
	shell: 
		""" 
		python {params.scripts_dir}/feature_tables/gen_num_sum_nearby_enhancers.py \
			--abc_predictions {input.abc_predictions} \
			--enhancer_list {input.enhancer_list} \
			--distance_threshold_kb {wildcards.kb} \
			--chr_sizes {params.chr_sizes} \
			--enh_midpoint {output.enhMidpoint} \
			--enh_expanded {output.enhExpanded} \
			--pred_slim {output.predSlim} \
			--enh_pred_int {output.enhPredInt} \
			--out_num {output.NumEnhancersEG} \
			--out_sum {output.SumEnhancersEG}

		"""

# create activity-only feature table
rule activity_only_features:
	input:
		feature_table_file = os.path.join(RESULTS_DIR, "{biosample}", "feature_table.tsv"),
		abc = lambda wildcards: os.path.join(ABC_BIOSAMPLES_DIR[wildcards.biosample], "Predictions", "EnhancerPredictionsAllPutative.tsv.gz"),
		NumCandidateEnhGene = get_numCandidateEnhGene_file,
		NumTSSEnhGene = get_numTSSEnhGene_file,
		NumEnhancersEG5kb = get_numNearbyEnhancers_file,
		SumEnhancersEG5kb = get_sumNearbyEnhancers_file,
		geneClasses = config["gene_classes"]
	output: 
		predictions_extended = os.path.join(RESULTS_DIR, "{biosample}", "ActivityOnly_features.tsv.gz")
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=determine_mem_mb
	script:
		"../scripts/feature_tables/activity_only_features.R"

# add external features
if config["final_score_col"] == "E2G.Score.qnorm": # if sc-E2G
	min_mem = 32
else:
	min_mem = 8

rule add_external_features:
	input:
		predictions_extended = os.path.join(RESULTS_DIR, "{biosample}", "ActivityOnly_features.tsv.gz"),
		feature_table_file = os.path.join(RESULTS_DIR, "{biosample}", "feature_table.tsv"),
		external_features_config = ancient(os.path.join(RESULTS_DIR, "{biosample}", "external_features_config.tsv"))
	output:
		plus_external_features = os.path.join(RESULTS_DIR, "{biosample}",  "ActivityOnly_plus_external_features.tsv.gz")
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=partial(determine_mem_mb, min_gb=min_mem)  
	script:
		"../scripts/feature_tables/merge_external_features.R"

# compute interaction or squared terms, fill NAs, rename features to finals, fill nas
rule gen_final_features:
	input:
		plus_external_features = os.path.join(RESULTS_DIR, "{biosample}", "ActivityOnly_plus_external_features.tsv.gz"),
		feature_table_file = os.path.join(RESULTS_DIR, "{biosample}", "feature_table.tsv")
	output:
		final_features = os.path.join(RESULTS_DIR, "{biosample}", "genomewide_features.tsv.gz")
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=determine_mem_mb
	script:
		"../scripts/feature_tables/gen_final_features.R"
