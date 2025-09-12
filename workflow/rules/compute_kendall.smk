from functools import partial

## compute Kendall correlation
rule compute_kendall:
	input:
		kendall_pairs_path = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"Kendall", 
				"Pairs.tsv.gz"
			),
		atac_matrix = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"Kendall", 
				"atac_matrix.rds"
			),
		rna_matrix = 
			lambda wildcards: CELL_CLUSTER_DF.loc[wildcards.cluster, "rna_matrix_file"],
	params:
		gene_gtf = config["gene_annotations"],
		abc_genes = config['gene_TSS500']
	output:
		kendall_predictions = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"Kendall", 
				"Pairs.Kendall.tsv.gz"),
		umi_count = (os.path.join(RESULTS_DIR, "{cluster}", "umi_count.txt")),
		cell_count = (os.path.join(RESULTS_DIR, "{cluster}", "cell_count.txt")),
		all_gex = os.path.join(RESULTS_DIR, "{cluster}", "Kendall", "gene_expression_metrics.tsv.gz")
	resources: 
		mem_mb=partial(encode_e2g.ABC.determine_mem_mb, min_gb=63),
		runtime=lambda wildcards, attempt: attempt*12*60
	threads: 
		config["threads"]
	conda:
		"../envs/sc_e2g.yml"
	script:
		"../scripts/feature_computation/compute_kendall.R"
