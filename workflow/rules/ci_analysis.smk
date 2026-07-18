## CI analysis rules
## Aggregates subsample predictions and computes CI statistics and QC reports

def get_subsample_predictions(wildcards):
    """Return list of prediction files for all subsamples of a given cluster."""
    return [
        os.path.join(
            SCRATCH_DIR, f"{wildcards.cluster}__sub{i}", MODEL,
            "encode_e2g_predictions.tsv.gz"
        )
        for i in SUBSAMPLE_INDICES
    ]

def get_subsample_params_json(wildcards):
    return os.path.join(SCRATCH_DIR, wildcards.cluster, "subsample_params.json")


## Collect per-subsample QC stats (cells, fragments, UMIs) from the pipeline QC output
rule collect_subsample_stats:
    input:
        qc_stats = os.path.join(SCRATCH_DIR, "qc_plots", "all_qc_stats.tsv")
    output:
        subsample_stats = os.path.join(CI_RESULTS_DIR, "{cluster}", MODEL, "subsample_qc_stats.tsv")
    run:
        import pandas as pd
        df = pd.read_table(input.qc_stats)
        pattern = wildcards.cluster + "__sub"
        df = df[df["cluster"].str.startswith(pattern)].copy()
        df["subsample_index"] = df["cluster"].str.extract(r"__sub(\d+)$")[0].astype(int)
        df = df.sort_values("subsample_index")
        keep_cols = [c for c in ["subsample_index", "cluster", "model_name", "cell_count", "fragments_total", "umi_count"] if c in df.columns]
        df[keep_cols].to_csv(output.subsample_stats, sep="\t", index=False)


## Compute CI annotations on the full run's predictions using subsample scores
rule compute_ci_annotations:
    input:
        full_predictions = lambda wildcards: os.path.join(
            FULL_RUN_DIR, wildcards.cluster, MODEL, "encode_e2g_predictions.tsv.gz"
        ),
        subsample_predictions = get_subsample_predictions,
        params_json = get_subsample_params_json
    params:
        model_threshold = lambda wildcards: CI_BIOSAMPLE_DF.loc[wildcards.cluster, "model_threshold"],
        ci_level = config.get("ci_level", 0.95),
        overlap_fraction = config.get("overlap_fraction", 0.50),
        cluster = lambda wildcards: wildcards.cluster,
        model = MODEL,
        mode = MODE,
        n_subsamples = N_SUBSAMPLES,
        scripts_dir = CI_SCRIPTS_DIR
    output:
        annotated = os.path.join(CI_RESULTS_DIR, "{cluster}", MODEL, "encode_e2g_predictions_ci_annotated.tsv.gz"),
        summary = os.path.join(CI_RESULTS_DIR, "{cluster}", MODEL, "ci_summary_stats.tsv"),
        gene_coverage = os.path.join(CI_RESULTS_DIR, "{cluster}", MODEL, "gene_coverage_stats.tsv")
    conda:
        "../envs/sc_e2g.yml"
    resources:
        mem_mb=64000,
        runtime=480
    script:
        "../scripts/ci/compute_ci_annotations.R"


## Generate consolidated CI QC report covering all clusters
rule plot_ci:
    input:
        annotated = expand(
            os.path.join(CI_RESULTS_DIR, "{cluster}", MODEL, "encode_e2g_predictions_ci_annotated.tsv.gz"),
            cluster=CLUSTERS
        ),
        summary = expand(
            os.path.join(CI_RESULTS_DIR, "{cluster}", MODEL, "ci_summary_stats.tsv"),
            cluster=CLUSTERS
        ),
        gene_coverage = expand(
            os.path.join(CI_RESULTS_DIR, "{cluster}", MODEL, "gene_coverage_stats.tsv"),
            cluster=CLUSTERS
        ),
        params_json = expand(
            os.path.join(SCRATCH_DIR, "{cluster}", "subsample_params.json"),
            cluster=CLUSTERS
        ),
        subsample_stats = expand(
            os.path.join(CI_RESULTS_DIR, "{cluster}", MODEL, "subsample_qc_stats.tsv"),
            cluster=CLUSTERS
        ),
        full_run_qc_stats = os.path.join(FULL_RUN_DIR, "qc_plots", "all_qc_stats.tsv")
    params:
        clusters = CLUSTERS,
        model = MODEL,
        mode = MODE,
        n_subsamples = N_SUBSAMPLES,
        model_thresholds = {c: CI_BIOSAMPLE_DF.loc[c, "model_threshold"] for c in CLUSTERS},
        scripts_dir = CI_SCRIPTS_DIR
    output:
        report = os.path.join(CI_RESULTS_DIR, "ci_qc_report.html")
    conda:
        "../envs/sc_e2g.yml"
    resources:
        mem_mb=64000,
        runtime=120
    script:
        "../scripts/ci/plot_ci.R"
