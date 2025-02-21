print(CELL_CLUSTER_DF["rna_matrix_file"].tolist())
## generate single-cell atac-seq matrix
rule generate_atac_matrix:
	input:
		kendall_pairs_path = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"Kendall", 
				"Pairs.tsv.gz"
			),
		atac_frag_path = 
			lambda wildcards: CELL_CLUSTER_DF.loc[wildcards.cluster, "atac_frag_file"],
		rna_matrix_path = 
			lambda wildcards: CELL_CLUSTER_DF.loc[wildcards.cluster, "rna_matrix_file"]
	output:
		atac_matrix_path = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"Kendall", 
				"atac_matrix.rds"
			)
	resources:
		mem_mb=128*1000
	conda:
		"../envs/sc_e2g.yml"
	script:
		"../scripts/feature_computation/generate_atac_matrix.R"
