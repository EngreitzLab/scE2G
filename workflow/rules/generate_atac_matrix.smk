## get list of cells defining this cluster
def get_cell_barcode_file(RNA_filt):
	if RNA_filt:
		# A pre-filtered RNA matrix needs no barcode list, so this branch has no
		# real file to point at and returns the results directory as a placeholder.
		# generate_atac_matrix.R guards the input with `file_test("-f", ...)`, which
		# is FALSE for a directory, so it is never read -- but Snakemake still tracks
		# its mtime, and a directory's mtime changes whenever an entry is added or
		# removed anywhere directly inside it. Any rule that creates a new
		# subdirectory under RESULTS_DIR therefore invalidates generate_atac_matrix
		# for EVERY cluster, cascading into compute_kendall and arc_e2g and
		# recomputing hours of work for byte-identical output.
		#
		# ancient() keeps the dependency while telling Snakemake to ignore the
		# timestamp, which is exactly right for an input the script never reads. A
		# missing output still triggers the job normally.
		return ancient(RESULTS_DIR)
	else:
		# A real file here, and its timestamp genuinely should trigger a rerun.
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
