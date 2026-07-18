from functools import partial

# Feature requirement detection via upfront configuration evaluation
# (Replaces checkpoint basic_features_required)

def get_numCandidateEnhGene_file(wildcards):
	"""Return feature file path if needed, otherwise RESULTS_DIR."""
	if BIOSAMPLE_NEEDS_NUMCANDIDATEENHGENE.get(wildcards.biosample, False):
		return os.path.join(RESULTS_DIR, wildcards.biosample, "new_features", "NumCandidateEnhGene.tsv")
	return RESULTS_DIR

def get_numTSSEnhGene_file(wildcards):
	"""Return feature file path if needed, otherwise RESULTS_DIR."""
	if BIOSAMPLE_NEEDS_NUMTSSENHGENE.get(wildcards.biosample, False):
		return os.path.join(RESULTS_DIR, wildcards.biosample, "new_features", "NumTSSEnhGene.tsv")
	return RESULTS_DIR

def get_numNearbyEnhancers_file(wildcards):
	"""Return feature file path if needed, otherwise RESULTS_DIR."""
	if BIOSAMPLE_NEEDS_NEARBYENHANCERS.get(wildcards.biosample, False):
		return os.path.join(RESULTS_DIR, wildcards.biosample, "new_features", "NumEnhancersEG5kb.txt")
	return RESULTS_DIR

def get_sumNearbyEnhancers_file(wildcards):
	"""Return feature file path if needed, otherwise RESULTS_DIR."""
	if BIOSAMPLE_NEEDS_NEARBYENHANCERS.get(wildcards.biosample, False):
		return os.path.join(RESULTS_DIR, wildcards.biosample, "new_features", "SumEnhancersEG5kb.txt")
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
