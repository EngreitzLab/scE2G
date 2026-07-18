## CI subsample generation rules
## Determines adaptive subsample fractions and generates subsampled fragment files.
## Uses FULL_RUN_CELL_CLUSTER_DF to access original (pre-subsampling) fragment files.

def get_frag_file(cluster):
    return FULL_RUN_CELL_CLUSTER_DF.loc[cluster, "atac_frag_file"]


## Determine adaptive subsample fraction based on cluster size
## Ensures each subsample meets minimum cell and fragment thresholds
rule determine_subsample_params:
    input:
        frag_file = lambda wildcards: get_frag_file(wildcards.cluster)
    params:
        target_fraction = config.get("target_fraction", 0.70),
        min_cells = config.get("min_cells", 100),
        min_fragments = config.get("min_fragments", 2000000)
    output:
        params_json = os.path.join(SCRATCH_DIR, "{cluster}", "subsample_params.json")
    conda:
        "../envs/sc_e2g.yml"
    shell:
        """
        python {CI_SCRIPTS_DIR}/determine_subsample_params.py \
            --frag_file {input.frag_file} \
            --target_fraction {params.target_fraction} \
            --min_cells {params.min_cells} \
            --min_fragments {params.min_fragments} \
            --output {output.params_json}
        """


## Generate a subsampled fragment file for one subsample
## Randomly selects a fraction of cell barcodes (adaptive fraction from step above)
rule subsample_cells:
    input:
        frag_file = lambda wildcards: get_frag_file(wildcards.cluster),
        params_json = os.path.join(SCRATCH_DIR, "{cluster}", "subsample_params.json")
    params:
        subsample_seed = config.get("subsample_seed", 42)
    output:
        frag_out = os.path.join(SCRATCH_DIR, "{cluster}", "subsample_{i}", "atac_fragments.tsv.gz"),
        frag_index = os.path.join(SCRATCH_DIR, "{cluster}", "subsample_{i}", "atac_fragments.tsv.gz.tbi")
    conda:
        "../envs/sc_e2g.yml"
    resources:
        mem_mb=16000,
        runtime=240
    shell:
        """
        unset SLURM_CPUS_PER_TASK
        python {CI_SCRIPTS_DIR}/subsample_cells.py \
            --frag_file {input.frag_file} \
            --params_json {input.params_json} \
            --subsample_index {wildcards.i} \
            --base_seed {params.subsample_seed} \
            --output {output.frag_out}
        """
