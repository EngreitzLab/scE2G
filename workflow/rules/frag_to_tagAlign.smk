## Fragment processing rules
## Handles conversion of fragment files to tagAlign format, bigWig generation, and indexing

def get_processed_fragment_file(wildcards):
	"""Return path to processed fragment file based on preprocessing mode."""
	if config.get("fragments_preprocessed", False):
		return CELL_CLUSTER_DF.loc[wildcards.cluster, "atac_frag_file"]
	else:
		return os.path.join(RESULTS_DIR, wildcards.cluster, "fragments_filtered.tsv.gz")


def get_processed_fragment_file_index(wildcards):
	"""Return path to index file for processed fragments."""
	if config.get("fragments_preprocessed", False):
		# When preprocessed, use ensure_fragment_index rule to create/link index
		return os.path.join(RESULTS_DIR, wildcards.cluster, "input_fragments.tsv.gz.tbi")
	else:
		return os.path.join(RESULTS_DIR, wildcards.cluster, "fragments_filtered.tsv.gz.tbi")


## Rule to ensure index exists for preprocessed fragment files
## Creates index on-the-fly if user didn't provide one
rule ensure_fragment_index:
	input:
		frag_file = lambda wildcards: CELL_CLUSTER_DF.loc[wildcards.cluster, "atac_frag_file"]
	output:
		index_file = os.path.join(RESULTS_DIR, "{cluster}", "input_fragments.tsv.gz.tbi")
	conda:
		"../envs/sc_e2g.yml"
	resources:
		mem_mb=ABC.determine_mem_mb
	shell:
		"""
		# Check if index exists alongside input file
		if [ -f "{input.frag_file}.tbi" ]; then
			ln -sf {input.frag_file}.tbi {output.index_file}
		else
			# Create index on-the-fly
			tabix -p bed {input.frag_file}
			ln -sf {input.frag_file}.tbi {output.index_file}
		fi
		"""


## Filter fragment file to valid chromosomes (only when fragments_preprocessed=False)
rule process_fragment_file:
	input:
		frag_file = lambda wildcards: CELL_CLUSTER_DF.loc[wildcards.cluster, "atac_frag_file"]
	params:
		chrSizes = config["chr_sizes"]
	output:
		fragments_filtered = os.path.join(RESULTS_DIR, "{cluster}", "fragments_filtered.tsv.gz"),
		fragments_filtered_index = os.path.join(RESULTS_DIR, "{cluster}", "fragments_filtered.tsv.gz.tbi")
	threads: 8
	resources:
		mem_mb=ABC.determine_mem_mb,
		runtime=720*2,
		temp_dir = os.path.join(RESULTS_DIR, "tmp")
	conda:
		"../envs/sc_e2g.yml"
	shell:
		"""
		LC_ALL=C
		# Filter fragments to valid chromosomes and compress
		awk 'NR==FNR {{keep[$1]; next}} $1 in keep' {params.chrSizes} <(zcat {input.frag_file}) | \
			bgzip > {output.fragments_filtered}

		# Create tabix index
		tabix -p bed {output.fragments_filtered}
		"""


## Convert fragment file to tagAlign file
## Also produces fragment count as additional output
rule frag_to_tagAlign:
	input:
		frag_file = get_processed_fragment_file,
		frag_file_index = get_processed_fragment_file_index
	output:
		tagAlign_sort_file = os.path.join(RESULTS_DIR, "{cluster}", "tagAlign", "tagAlign.sort.gz"),
		tagAlign_sort_index = os.path.join(RESULTS_DIR, "{cluster}", "tagAlign", "tagAlign.sort.gz.tbi"),
		fragment_count = os.path.join(RESULTS_DIR, "{cluster}", "fragment_count.txt")
	conda:
		"../envs/sc_e2g.yml"
	threads: 8
	resources:
		mem_mb=ABC.determine_mem_mb,
		runtime=720*2,
		temp_dir = os.path.join(RESULTS_DIR, "tmp")
	shell:
		"""
		# Calculate buffer size for sort
		export BUFFER_SIZE=$(awk -v mem_mb={resources.mem_mb} -v threads={threads} 'BEGIN {{ result = mem_mb/threads/2; print int(result) }}')

		LC_ALL=C
		# Count fragments while converting to tagAlign format
		# tee to wc -l counts input lines, then awk converts to tagAlign (2 rows per fragment)
		zcat {input.frag_file} | \
			tee >(wc -l > {output.fragment_count}) | \
			awk -v OFS='\t' '{{mid=int(($2+$3)/2); print $1,$2,mid,"N",1000,"+"; print $1,mid+1,$3,"N",1000,"-"}}' | \
			sort -k 1,1V -k 2,2n -k3,3n --parallel {threads} -T {resources.temp_dir} -S $BUFFER_SIZE | \
			bgzip -c > {output.tagAlign_sort_file}

		# Index the tagAlign file
		tabix -p bed {output.tagAlign_sort_file}
		"""


## Create bigWig coverage tracks from fragment file
## Generates both unnormalized and CPM-normalized versions in a single pass
rule frag_to_bigwig:
	input:
		frag_file = get_processed_fragment_file,
		fragment_count = os.path.join(RESULTS_DIR, "{cluster}", "fragment_count.txt")
	params:
		chrSizes = config["chr_sizes"]
	output:
		bigWig_unnorm = os.path.join(IGV_DIR, "{cluster}", "ATAC.bw"),
		bigWig_norm = os.path.join(IGV_DIR, "{cluster}", "ATAC_norm.bw"),
		bedGraph_file = temp(os.path.join(IGV_DIR, "{cluster}", "ATAC.bg"))
	resources:
		mem_mb=ABC.determine_mem_mb,
		runtime_hr=24,
		temp_dir = os.path.join(RESULTS_DIR, "tmp")
	threads: 16
	conda:
		"../envs/sc_e2g.yml"
	shell:
		"""
		LC_ALL=C
		export BUFFER_SIZE=$(awk -v mem_mb={resources.mem_mb} -v threads={threads} 'BEGIN {{ result = mem_mb/threads/2; print int(result) }}')

		# Generate sorted bedgraph once (unnormalized coverage)
		zcat {input.frag_file} | \
			bedtools genomecov -bg -i stdin -g {params.chrSizes} | \
			sort -k1,1 -k2,2n --parallel={threads} -S $BUFFER_SIZE -T {resources.temp_dir} > {output.bedGraph_file}

		# Create unnormalized bigWig
		bedGraphToBigWig {output.bedGraph_file} {params.chrSizes} {output.bigWig_unnorm}

		# Create CPM-normalized bigWig (scale bedgraph values by 1M/fragment_count)
		frag_count=$(<{input.fragment_count})
		scale_factor=$(awk "BEGIN {{print 1000000 / $frag_count}}")
		norm_bg=$(mktemp --suffix=.bg)
		awk -v sf=$scale_factor 'BEGIN{{OFS="\t"}} {{print $1,$2,$3,$4*sf}}' {output.bedGraph_file} > "$norm_bg"
		bedGraphToBigWig "$norm_bg" {params.chrSizes} {output.bigWig_norm}
		rm -f "$norm_bg"
		"""
