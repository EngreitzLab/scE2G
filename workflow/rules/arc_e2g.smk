## integrate Kendall and ABC into ARC-E2G
rule arc_e2g:
	input:
		abc_predictions = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"Predictions", 
				"EnhancerPredictionsAllPutative.tsv.gz"
			),
		kendall_predictions = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"Kendall", 
				"Pairs.Kendall.tsv.gz")
	params:
	 	abc_score_col = lambda wildcards: get_abc_score_col(wildcards.cluster)
	output:
		arc_predictions = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"ARC",
				"EnhancerPredictionsAllPutative_ARC.tsv.gz")
	resources: 
		mem_mb=64*1000
	conda:
		"../envs/sc_e2g.yml"
	script:
		"../scripts/feature_computation/compute_arc_e2g.R"
