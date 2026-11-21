# save true cell type barcodes
rule save_reference_barcodes:
	input:
		L3_fragments = lambda wildcards: CT_CONFIG.loc[wildcards.l3_id, "L3_fragment_file"]
	output: 
		cell_barcodes = os.path.join(RESULTS_DIR_OAK, "L3_reference_barcodes", "{l3_id}_barcodes.txt.gz")
	resources:
		mem_mb=determine_mem_mb
	shell:
		"""
		zcat {input.L3_fragments} | cut -f 4 | sort | uniq | gzip > {output.cell_barcodes}
		"""
		

# calculate number of cells per sample
rule calculate_cell_type_counts:
	input:
		ref_barcodes = expand(os.path.join(RESULTS_DIR_OAK, "L3_reference_barcodes", "{l3_id}_barcodes.txt.gz"), l3_id = CT_CONFIG["L3"])
	params:
		ct_config = config["cell_type_config"],
		n_sample = config["n_sample"],
		n_L1_sample = config["n_L1_sample"]
	output:
		l1_counts = os.path.join(RESULTS_DIR_OAK, "reference", "L3_counts_for_L1.pkl"), # cells per ct for L1 sample (from the main sample)
		l2_counts = os.path.join(RESULTS_DIR_OAK, "reference", "L2_counts.pkl"), # cells per sg
		l3_counts = os.path.join(RESULTS_DIR_OAK, "reference", "L3_counts.pkl") # cells per ct
	resources:
		mem_mb=determine_mem_mb
	run:
		import gzip
		ct_config = pd.read_csv(params.ct_config, sep = "\t").set_index("L3", drop = False)

		# total counts per cell type
		total_counts = {l3_id: sum(1 for _ in gzip.open(barcode_file, "rt")) for l3_id, barcode_file in zip(ct_config["L3"], input.ref_barcodes)}
		total_cells = sum(total_counts.values())
		sample_frac = params.n_sample / total_cells
		l1_frac = params.n_L1_sample / total_cells

		# save l1 counts and l3 counts
		l3_counts = {l3_id: round(counts * sample_frac) for l3_id, counts in total_counts.items()}
		l1_counts = {l3_id: round(counts * l1_frac) for l3_id, counts in total_counts.items()}
		with open(output.l3_counts, "wb") as file: pickle.dump(l3_counts, file)
		with open(output.l1_counts, "wb") as file: pickle.dump(l1_counts, file)

		# save l2 counts
		l2_dict = ct_config.groupby("L2")["L3"].apply(list).to_dict()
		l2_counts = {l2_id: sum(l3_counts[l3_id] for l3_id in l2_dict[l2_id]) for l2_id in l2_dict.keys()}
		with open(output.l2_counts, "wb") as file: pickle.dump(l2_counts, file)


rule get_barcodes_for_iteration:
	input: 
		ref_barcodes = expand(os.path.join(RESULTS_DIR_OAK, "L3_reference_barcodes", "{l3_id}_barcodes.txt.gz"), l3_id = CT_CONFIG["L3"]),
		l1_counts = os.path.join(RESULTS_DIR_OAK, "reference", "L3_counts_for_L1.pkl"), # cells per ct for L1 sample (from the main sample)
		l2_counts = os.path.join(RESULTS_DIR_OAK, "reference", "L2_counts.pkl"), # cells per sg
		l3_counts = os.path.join(RESULTS_DIR_OAK, "reference", "L3_counts.pkl") # cells per ct
	params:
		scratch_dir = SCRATCH_DIR,
		ct_config = config["cell_type_config"],
		seed_multiplier = config["seed_multiplier"]
	output:
		#l1_merge = os.path.join(SCRATCH_DIR, "iter{iter}", "real", "L1", "all", "sample_barcodes.txt.gz"),
		l3_from_l2_out = expand(os.path.join(RESULTS_DIR, "iter{{iter}}.{sample_type}.L3_from_L2.{l3_id}", "sample_barcodes.txt.gz"),
			sample_type = ["permute"], l3_id = L3_L2_DICT.keys()),
		l3_from_l1_out = expand(os.path.join(RESULTS_DIR, "iter{{iter}}.{sample_type}.L3_from_L1.{l3_id}", "sample_barcodes.txt.gz"),
			sample_type = ["permute"], l3_id = L3_L2_DICT.keys()),
		l3_real_out = expand(os.path.join(RESULTS_DIR, "iter{{iter}}.{sample_type}.L3.{l3_id}", "sample_barcodes.txt.gz"),
			sample_type = ["real"], l3_id = L3_L2_DICT.keys()),
	resources:
		mem_mb=determine_mem_mb
	run:
		import random
		import gzip
		random.seed(wildcards.iter * params.seed_multiplier)

		# read in count dicts, configs
		with open(input.l1_counts, "rb") as file: l1_counts = pickle.load(file)
		with open(input.l2_counts, "rb") as file: l2_counts = pickle.load(file)
		with open(input.l3_counts, "rb") as file: l3_counts = pickle.load(file)
		ct_config = pd.read_csv(params.ct_config, sep = "\t").set_index("L3", drop = False)
		l2_dict = ct_config.groupby("L2")["L3"].apply(list).to_dict()

		# get full sample & real subsets
		all_barcodes = {l3_id: gzip.open(barcode_file, "rt").read().splitlines() for l3_id, barcode_file in zip(ct_config["L3"], input.ref_barcodes)}
		sample_l3 = {l3_id: random.sample(all_barcodes[l3_id], n) for l3_id, n in l3_counts.items()}
		#sample_merge = list(itertools.chain(*sample_l3.values())) 
		sample_l2 = {l2_id: list(itertools.chain(*(sample_l3[l3_id] 
			for l3_id in l2_l3s)))
			for l2_id, l2_l3s in l2_dict.items()}
		# sample_l1 = [random.sample(sample_l3[l3_id], l1_counts[l3_id]) for l3_id in sample_l3.keys()]
		# sample_l1_merge = list(itertools.chain(*sample_l1))

		# # get permuted samples for l2
		# sample_merge = set(itertools.chain(*sample_l3.values())) 
		# permuted_l2 = {}
		# for l2_id in l2_dict.keys():
		# 	permuted_l2[l2_id] = random.sample(sorted(sample_merge), l2_counts[l2_id])
		# 	this_set = set(permuted_l2[l2_id]) # convert to set for faster filtering
		# 	sample_merge.difference_update(permuted_l2[l2_id])
		
		# get permuted samples for l3 from l1
		sample_merge = set(itertools.chain(*sample_l3.values())) 
		permuted_l3_from_l1 = {}
		for l3_id in l3_counts.keys():
			permuted_l3_from_l1[l3_id] = random.sample(sorted(sample_merge), l3_counts[l3_id])
			this_set = set(permuted_l3_from_l1[l3_id]) # convert to set for faster filtering
			sample_merge.difference_update(permuted_l3_from_l1[l3_id])

		# get permuted samples for l3 from l2
		sample_merge = set(itertools.chain(*sample_l3.values())) 
		permuted_l3_from_l2 = {}
		for l2_id, l2_barcodes in sample_l2.items():
			this_l2_set = set(l2_barcodes)
			for l3_id in l2_dict[l2_id]:
				permuted_l3_from_l2[l3_id] = random.sample(sorted(this_l2_set), l3_counts[l3_id])
				this_set = set(permuted_l3_from_l2[l3_id]) # convert to set for faster filtering
				this_l2_set.difference_update(permuted_l3_from_l2[l3_id])

		# # get permuted sample for l3 from l2
		# l3_dict = {l3_id: l2_id for l2_id, l3_ids in l2_dict.items() for l3_id in l3_ids}
		# permuted_l3_from_l2  = {l3_id: random.sample(sample_l2[l3_dict[l3_id]], l3_counts[l3_id]) for l3_id in sample_l3.keys()}

		# function to save files
		def save_barcodes(out_dir, iter, barcodes, sample_type, level, id):
			cluster_id = f"iter{iter}.{sample_type}.{level}.{id}"
			out_file = os.path.join(out_dir, cluster_id, "sample_barcodes.txt.gz")
			with gzip.open(out_file, "wt") as f:
				for bc in barcodes:
					f.write(f"{bc}\n")

		# save files
		[save_barcodes(params.scratch_dir, wildcards.iter, barcodes, "real", "L3", l3_id) for l3_id, barcodes in sample_l3.items()]
		#[save_barcodes(params.scratch_dir, wildcards.iter, barcodes, "real", "L2", l2_id) for l2_id, barcodes in sample_l2.items()]
		[save_barcodes(params.scratch_dir, wildcards.iter, barcodes, "permute", "L3_from_L2", l3_id) for l3_id, barcodes in permuted_l3_from_l2.items()]
		[save_barcodes(params.scratch_dir, wildcards.iter, barcodes, "permute", "L3_from_L1", l3_id) for l3_id, barcodes in permuted_l3_from_l1.items()]
		#save_barcodes(out_dir, sample_l1_merge, "real", "L1", "all")


def get_input_file(file_type, sample_type, level, id):
	if level == "L1":
		res = CT_CONFIG["L1_" + file_type].dropna().iloc[0]
	elif sample_type == "real":
		res = CT_CONFIG.loc[CT_CONFIG[level] == id][level + "_" + file_type].dropna().iloc[0]
	elif sample_type == "permute":
		if level == "L2":
			res = CT_CONFIG["L1_" + file_type].dropna().iloc[0]
		elif level == "L3_from_L2":
			res = CT_CONFIG.loc[CT_CONFIG["L2"] == L3_L2_DICT[id]]["L2_" + file_type].dropna().iloc[0] # get the corresponding L2 file
		elif level == "L3_from_L1":
			res = CT_CONFIG["L1_" + file_type].dropna().iloc[0]
	return res
		

rule create_sample_fragment_file:
	input:
		fragments = lambda wildcards: get_input_file("fragment_file", wildcards.sample_type, wildcards.level, wildcards.sample_id),
		sample_barcodes = os.path.join(RESULTS_DIR, "iter{iter}.{sample_type}.{level}.{sample_id}", "sample_barcodes.txt.gz")
	output: 
		sample_fragments = os.path.join(RESULTS_DIR, "iter{iter}.{sample_type}.{level}.{sample_id}", "atac_fragments.tsv.gz"),
		sample_index = os.path.join(RESULTS_DIR, "iter{iter}.{sample_type}.{level}.{sample_id}", "atac_fragments.tsv.gz.tbi")
	resources:
		mem_mb=determine_mem_mb
	threads: 8
	conda:
		"../envs/sc_e2g.yml"
	shell:
		"""
			zgrep -F -f <(zcat {input.sample_barcodes}) {input.fragments} | \
				sort -k1,1 -k2,2n --parallel {threads} | \
				bgzip > {output.sample_fragments}
			
			tabix -p bed {output.sample_fragments}
		"""


rule create_sample_rna_matrix:
	input:
		matrix = lambda wildcards: get_input_file("rna_matrix", wildcards.sample_type, wildcards.level, wildcards.sample_id),
		sample_barcodes = os.path.join(RESULTS_DIR, "iter{iter}.{sample_type}.{level}.{sample_id}", "sample_barcodes.txt.gz")
	output: 
		sample_rna = os.path.join(RESULTS_DIR, "iter{iter}.{sample_type}.{level}.{sample_id}", "rna_matrix.csv.gz"),
	resources:
		mem_mb=determine_mem_mb
	conda:
		"../envs/sc_e2g.yml"
	shell:
		"""
			cols=1,$(grep -Fx -n -f <(zcat  {input.sample_barcodes}) <(zcat {input.matrix}  | head -1 | tr ',' '\n') | cut -d: -f1 | paste -sd,)
			zcat {input.matrix}  | cut -d',' -f"$cols" | gzip > {output.sample_rna}
		"""