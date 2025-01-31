## Convert fragment file to tagAlign file
rule frag_to_tagAlign:
	input:
		frag_file = 
			lambda wildcards: CELL_CLUSTER_DF.loc[wildcards.cluster, "atac_frag_file"]
	output:
		tagAlign_sort_file = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"tagAlign",
				"tagAlign.sort.gz"
			),
		fragment_count = 
			temp(os.path.join(
				RESULTS_DIR,
				"{cluster}",
				"fragment_count.txt"
			)),
		cell_count = 
			temp(os.path.join(
				RESULTS_DIR,
				"{cluster}",
				"cell_count.txt"
			))
		processed_fragments = temp(os.path.join(RESULTS_DIR, "{cluster}", "atac_fragments.chr_filtered.sorted.tsv.gz"))
	params:
		chrSizes = encode_e2g_config["chr_sizes"]
	conda:
		"../envs/sc_e2g.yml"
	threads: 8
	resources:
		mem_mb=determine_mem_mb,
		runtime=720*2,
		temp_dir = 
			os.path.join(
				RESULTS_DIR,
				"tmp"
		)
	shell:
		"""
		# Make, sort and compress tagAlign file from fragment file
		LC_ALL=C 
		
		# save chr-filtered & sorted fragment file for later use
		chr_list=$(cut -f1 {params.chrSizes} | sort | uniq)
		zcat {input.frag_file} | sed '/^#/d' | awk -v chr_list="$chr_list" '$1 ~ chr_list' | \ # remove header, filter chromosomes
			sort -k 1,1V -k 2,2n -k3,3n --parallel {threads} -T {resources.temp_dir} | \ # sort
			gzip > {output.processed_fragments}

		zcat {output.processed_fragments} | \
			awk -v OFS='\t' '{{mid=int(($2+$3)/2); print $1,$2,mid,"N",1000,"+"; print $1,mid+1,$3,"N",1000,"-"}}' | \
			bgzip -c > {output.tagAlign_sort_file}  

		# get fragment & cell count (of unfiltered fragments)
		zcat {input.frag_file} | awk '$1 !~ /_/' | wc -l > {output.fragment_count}
		zcat {input.frag_file} | awk '$1 !~ /_/' | cut -f4 | sort -u | wc -l > {output.cell_count}

		# Index the tagAlign file
		tabix -p bed {output.tagAlign_sort_file}
		"""

## create bigwig from fragment file
rule frag_to_bigWig:
	input:
		processed_fragments = temp(os.path.join(RESULTS_DIR, "{cluster}", "atac_fragments.chr_filtered.sorted.tsv.gz"))
	params:
		chrSizes = encode_e2g_config["chr_sizes"]
	output:
		bigWig_file = os.path.join(IGV_DIR, "{cluster}", "ATAC.bw"),
		bedGraph_file = temp(os.path.join(IGV_DIR, "{cluster}", "ATAC.bg"))
	resources:
		mem_mb=determine_mem_mb,
		runtime_hr=24
	threads: 16
	conda: 
		"../envs/sc_e2g.yml"
	shell:
		"""
			LC_ALL=C
			zcat {input.processed_fragments} | \
				bedtools genomecov -bg -i stdin -g {params.chrSizes} > {output.bedGraph_file}
				#sort -k1,1 -k2,2n --parallel={threads} > {output.bedGraph_file} # redundant now
			bedGraphToBigWig {output.bedGraph_file} {params.chrSizes} {output.bigWig_file}
		"""

rule frag_to_norm_bigWig:
	input:
		processed_fragments = temp(os.path.join(RESULTS_DIR, "{cluster}", "atac_fragments.chr_filtered.sorted.tsv.gz")),
		fragment_count = os.path.join(RESULTS_DIR, "{cluster}", "fragment_count.txt")
	params:
		chrSizes = encode_e2g_config["chr_sizes"]
	output:
		bigWig_file = os.path.join(IGV_DIR, "{cluster}", "ATAC_norm.bw"),
		bedGraph_file = temp(os.path.join(IGV_DIR, "{cluster}", "ATAC_norm.bg"))
	resources:
		mem_mb=determine_mem_mb,
		runtime_hr=24
	threads: 16
	conda: 
		"../envs/sc_e2g.yml"
	shell:
		"""
			LC_ALL=C
			frag_count=$(<{input.fragment_count})
			scale_factor=$(echo "1000000 / $frag_count" | bc -l)

			zcat {input.processed_fragments} | \
				bedtools genomecov -bg -i stdin -g {params.chrSizes} -scale $scale_factor > {output.bedGraph_file}
				#sort -k1,1 -k2,2n --parallel={threads} > {output.bedGraph_file} # redundant now
			bedGraphToBigWig {output.bedGraph_file} {params.chrSizes} {output.bigWig_file}
		"""
