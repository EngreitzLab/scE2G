def get_processed_fragment_file(wildcards):
	if config["fragments_preprocessed"]:
		fp = CELL_CLUSTER_DF.loc[wildcards.cluster, "atac_frag_file"]
	else:
		fp = os.path.join(RESULTS_DIR, wildcards.cluster, "fragments_filtered.tsv.gz")
	return fp

## Convert fragment file to tagAlign file
rule frag_to_tagAlign:
	input:
		frag_file = get_processed_fragment_file
	output:
		tagAlign_sort_file = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"tagAlign",
				"tagAlign.sort.gz"
			)
	conda:
		"../envs/sc_e2g.yml"
	threads: 8
	resources:
		mem_mb=determine_mem_mb,
		runtime=720*2,
		temp_dir = 
			os.path.join(
				RESULTS_DIR, "tmp"
		)
	shell:
		"""
		# Make, sort and compress tagAlign file from fragment file
		export BUFFER_SIZE=$(awk -v mem_mb={resources.mem_mb} -v threads={threads} 'BEGIN {{ result = mem_mb/threads/2; print int(result) }}')

		LC_ALL=C
		zcat {input.frag_file} | \
			awk -v OFS='\t' '{{mid=int(($2+$3)/2); print $1,$2,mid,"N",1000,"+"; print $1,mid+1,$3,"N",1000,"-"}}' | \
			sort -k 1,1V -k 2,2n -k3,3n --parallel {threads} -T {resources.temp_dir} -S $BUFFER_SIZE | \
		bgzip -c > {output.tagAlign_sort_file}  

		# Index the tagAlign file
		tabix -p bed {output.tagAlign_sort_file}
		"""

## get fragment cell count
if config["fragments_preprocessed"]:
	rule get_fragment_count:
		input:
			frag_file = lambda wildcards: CELL_CLUSTER_DF.loc[wildcards.cluster, "atac_frag_file"]
		resources: mem_mb=determine_mem_mb
		output:
			fragment_count = (os.path.join(RESULTS_DIR, "{cluster}", "fragment_count.txt")),
		shell:
			"""
				zcat {input.frag_file} | wc -l > {output.fragment_count}
			"""
else:
	rule process_fragment_file:
		input:
			frag_file = lambda wildcards: CELL_CLUSTER_DF.loc[wildcards.cluster, "atac_frag_file"]
		params:
			chrSizes = config["chr_sizes"]
		output:
			fragment_count = (os.path.join(RESULTS_DIR, "{cluster}", "fragment_count.txt")),
			fragments_filtered = (os.path.join(RESULTS_DIR, "{cluster}", "fragments_filtered.tsv.gz"))
		threads: 8
		resources:
			mem_mb=determine_mem_mb,
			runtime=720*2,
			temp_dir = 
				os.path.join(
					RESULTS_DIR,
					"tmp"
			)
		conda: 
			"../envs/sc_e2g.yml"
		shell:
			"""
			LC_ALL=C 
			# get fragment & cell count
			awk 'NR==FNR {{keep[$1]; next}} $1 in keep' {params.chrSizes} <(zcat {input.frag_file})  | bgzip > {output.fragments_filtered}
			tabix -p bed {output.fragments_filtered}

			zcat {output.fragments_filtered} | wc -l > {output.fragment_count}
			"""

## create bigwig from fragment file
rule frag_to_bigWig:
	input:
		frag_file = get_processed_fragment_file
	params:
		chrSizes = config["chr_sizes"]
	output:
		bigWig_file = os.path.join(IGV_DIR, "{cluster}", "ATAC.bw"),
		bedGraph_file = temp(os.path.join(IGV_DIR, "{cluster}", "ATAC.bg"))
	resources:
		mem_mb=determine_mem_mb,
		runtime_hr=24,
		temp_dir = os.path.join(RESULTS_DIR, "tmp")
	threads: 16
	conda: 
		"../envs/sc_e2g.yml"
	shell:
		"""
			LC_ALL=C
			export BUFFER_SIZE=$(awk -v mem_mb={resources.mem_mb} -v threads={threads} 'BEGIN {{ result = mem_mb/threads/2; print int(result) }}')
			zcat {input.frag_file} | \
				bedtools genomecov -bg -i stdin -g {params.chrSizes} | \
				sort -k1,1 -k2,2n --parallel={threads} -S $BUFFER_SIZE -T {resources.temp_dir} > {output.bedGraph_file}
			bedGraphToBigWig {output.bedGraph_file} {params.chrSizes} {output.bigWig_file}
		"""

rule frag_to_norm_bigWig:
	input:
		frag_file = get_processed_fragment_file,
		fragment_count = os.path.join(RESULTS_DIR, "{cluster}", "fragment_count.txt")
	params:
		chrSizes = config["chr_sizes"]
	output:
		bigWig_file = os.path.join(IGV_DIR, "{cluster}", "ATAC_norm.bw"),
		bedGraph_file = temp(os.path.join(IGV_DIR, "{cluster}", "ATAC_norm.bg"))
	resources:
		mem_mb=determine_mem_mb,
		runtime_hr=24,
		temp_dir=os.path.join(RESULTS_DIR, "tmp")
	threads: 16
	conda: 
		"../envs/sc_e2g.yml"
	shell:
		"""
			LC_ALL=C
			export BUFFER_SIZE=$(awk -v mem_mb={resources.mem_mb} -v threads={threads} 'BEGIN {{ result = mem_mb/threads/2; print int(result) }}')
			frag_count=$(<{input.fragment_count})
			scale_factor=$(awk "BEGIN {{print 1000000 / $frag_count}}")

			zcat {input.frag_file} | \
				bedtools genomecov -bg -i stdin -g {params.chrSizes} -scale $scale_factor| \
				sort -k1,1 -k2,2n --parallel={threads} -S $BUFFER_SIZE -T {resources.temp_dir} > {output.bedGraph_file}
			bedGraphToBigWig {output.bedGraph_file} {params.chrSizes} {output.bigWig_file}
		"""
