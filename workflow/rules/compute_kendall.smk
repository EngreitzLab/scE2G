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
			lambda wildcards: CELL_CLUSTER_DF.loc[wildcards.cluster, "rna_matrix_file"]
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
		umi_count = temp(os.path.join(RESULTS_DIR, "{cluster}", "umi_count.txt")) 
	resources: 
		mem_mb=200*1000
	conda:
		"../envs/sc_e2g.yml"
	script:
		"../scripts/feature_computation/compute_kendall.R"
