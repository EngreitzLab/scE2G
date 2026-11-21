# comparisons to make...
	# real 15 to permute 15
	# real 5 to permute 5
	# real 15 to permute 15 per supergroup??


rule make_pred_intersectable:
	input: 
		pred_thresholded = os.path.join(RESULTS_DIR, "{cluster}", "{model_name}", "encode_e2g_predictions_thresholded.tsv.gz")
	output:
		pred_intersectable = os.path.join(RESULTS_DIR, "{cluster}", "{model_name}", "gene_start_end_thresholded.tsv.gz")
	conda:
		"../envs/sc_e2g.yml"
	resources:
		mem_mb=determine_mem_mb
	shell:
		"""
			zcat {input.pred_thresholded} | csvtk cut -t -f TargetGene,start,end | \
				sed 1d | sort -k 1,1 -k2,2n | uniq | \
				gzip > {output.pred_intersectable}
		"""

rule compute_L1_totals:
	input:
		# real_L2 = expand(os.path.join(RESULTS_DIR, "iter{{iter}}.real.L2.{l2_id}", os.path.basename(config["model_dir_use"]), "gene_start_end_thresholded.tsv.gz"),
		# 	l2_id = L2_L3_DICT.keys()),
		# permute_L2 = expand(os.path.join(RESULTS_DIR, "iter{{iter}}.permute.L2.{l2_id}", os.path.basename(config["model_dir_use"]), "gene_start_end_thresholded.tsv.gz"),
		# 	l2_id = L2_L3_DICT.keys()),
		real_L3 = expand(os.path.join(RESULTS_DIR, "iter{{iter}}.real.L3.{l3_id}", os.path.basename(config["model_dir_use"]), "gene_start_end_thresholded.tsv.gz"),
			l3_id = L3_L2_DICT.keys()),
		permute_L3_from_L2 = expand(os.path.join(RESULTS_DIR, "iter{{iter}}.permute.L3_from_L2.{l3_id}", os.path.basename(config["model_dir_use"]), "gene_start_end_thresholded.tsv.gz"),
			l3_id = L3_L2_DICT.keys()),
		permute_L3_from_L1 = expand(os.path.join(RESULTS_DIR, "iter{{iter}}.permute.L3_from_L1.{l3_id}", os.path.basename(config["model_dir_use"]), "gene_start_end_thresholded.tsv.gz"),
			l3_id = L3_L2_DICT.keys()),
	output:
		# real_L2 = temp(os.path.join(RESULTS_DIR, "intersections", "iter{iter}", "L2_real.txt")),
		# permute_L2 = temp(os.path.join(RESULTS_DIR, "intersections", "iter{iter}", "L2_permute.txt")),
		real_L3 = (os.path.join(RESULTS_DIR, "intersections", "iter{iter}", "L3_real.txt")),
		permute_L3_from_L1 = (os.path.join(RESULTS_DIR, "intersections", "iter{iter}", "L3_from_L1_permute.txt")),
		permute_L3_from_L2 = (os.path.join(RESULTS_DIR, "intersections", "iter{iter}", "L3_from_L2_permute.txt")),
		totals = os.path.join(RESULTS_DIR, "intersections", "iter{iter}", "L1_totals.tsv")
	conda:
		"../envs/sc_e2g.yml"
	resources:
		mem_mb=determine_mem_mb
	shell:
		"""
			echo "iter L3_real_total L3_from_L2_total L3_from_L1_total" | \
				tr " " "\t" > {output.totals}
		
		## L3
			real_l3_array=($(echo {input.real_L3} | tr ' ' '\n' | shuf))
			cat ${{real_l3_array[0]}} > {output.real_L3}
			for pred in "${{real_l3_array[@]:1}}"; do
				bedtools intersect -a $pred -b {output.real_L3} -f 0.5 -F 0.5 -wa -v >> {output.real_L3}
			done

			permute_l3_from_l1_array=($(echo {input.permute_L3_from_L1} | tr ' ' '\n' | shuf))
			cat ${{permute_l3_from_l1_array[0]}} > {output.permute_L3_from_L1}
			for pred in "${{permute_l3_from_l1_array[@]:1}}"; do
				bedtools intersect -a $pred -b {output.permute_L3_from_L1} -f 0.5 -F 0.5 -wa -v >> {output.permute_L3_from_L1}
			done

			permute_l3_from_l2_array=($(echo {input.permute_L3_from_L2} | tr ' ' '\n' | shuf))
			cat ${{permute_l3_from_l2_array[0]}} > {output.permute_L3_from_L2}
			for pred in "${{permute_l3_from_l2_array[@]:1}}"; do
				bedtools intersect -a $pred -b {output.permute_L3_from_L2} -f 0.5 -F 0.5 -wa -v >> {output.permute_L3_from_L2}
			done

			real_L3_total=$(cat {output.real_L3} | wc -l)
			permute_L3_from_L1_total=$(cat {output.permute_L3_from_L1} | wc -l)
			permute_L3_from_L2_total=$(cat {output.permute_L3_from_L2} | wc -l)

			echo "{wildcards.iter} $real_L3_total $permute_L3_from_L2_total $permute_L3_from_L1_total" | \
				tr " " "\t" >> {output.totals}
		"""


rule combine_L1_totals:
	input:
		totals_list = expand(os.path.join(RESULTS_DIR, "intersections", "iter{i}", "L1_totals.tsv"),
			i = list(range(1, config["n_iterations"] + 1)))
	output:
		all_totals = os.path.join(RESULTS_DIR_OAK, "L1_totals.tsv")
	conda:
		"../envs/sc_e2g.yml"
	resources:
		mem_mb=determine_mem_mb
	shell:
		"""
			all_files=({input.totals_list})

			cat "${{all_files[0]}}" > {output.all_totals}
			tail -n +2 -q "${{all_files[@]:1}}" >> {output.all_totals}

		"""
		
