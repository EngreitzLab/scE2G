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
			"Paires.tsv.gz"
		)

	conda:
		"../envs/sc_e2g.yml"
	script:
		"../scripts/make_kendall_pairs.R"

