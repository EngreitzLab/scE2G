## Make enhancer-gene pairs for computing Kendall correlation
## using non-extended peaks
rule make_kendall_pairs:
	input:
		narrowPeak = os.path.join(
			RESULTS_DIR, 
			"{cluster}", 
			"Peaks", 
			"macs2_peaks.narrowPeak.sorted"
		),
		allPutative = os.path.join(
			RESULTS_DIR, 
			"{cluster}", 
			"Predictions", 
			"EnhancerPredictionsAllPutative.tsv.gz"
		)
	output:
		kendallPairs = os.path.join(
			RESULTS_DIR, 
			"{cluster}", 
			"Kendall", 
			"Pairs.tsv.gz"
		)
	resources:
		mem_mb=32*1000
	conda:
		"../envs/sc_e2g.yml"
	script:
		"../scripts/feature_computation/make_kendall_pairs.R"


