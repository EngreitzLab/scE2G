## get list of cells defining this cluster
def get_cell_barcode_file(RNA_filt):
	if RNA_filt:
		return RESULTS_DIR
	else:
		return os.path.join(RESULTS_DIR, "{cluster}", "Kendall", "cell_barcodes.txt")

rule get_cell_barcodes:
	input:
		frag_file = get_processed_fragment_file
	resources:
		mem_mb = encode_e2g.ABC.determine_mem_mb,
	output:
		cell_barcodes = os.path.join(RESULTS_DIR, "{cluster}", "Kendall", "cell_barcodes.txt")
	shell:
		"""
		zcat {input.frag_file} | cut -f 4 | awk '!seen[$0]++' > {output.cell_barcodes}
		"""

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
			lambda wildcards: CELL_CLUSTER_DF.loc[wildcards.cluster, "rna_matrix_file"],
		cell_barcodes_path = get_cell_barcode_file(config["RNA_matrix_filtered"])
	output:
		atac_matrix_path = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"Kendall", 
				"atac_matrix.rds"
			)
	params:
		max_cell_count = config['max_cell_count']
	resources:
		mem_mb=encode_e2g.ABC.determine_mem_mb
	conda:
		"../envs/sc_e2g.yml"
	script:
		"../scripts/feature_computation/generate_atac_matrix.R"
