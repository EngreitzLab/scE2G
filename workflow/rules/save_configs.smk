rule save_reference_configs:
    """
    Save resolved config files used for scE2G, ENCODE-rE2G, and ABC
    Save biosample / cell cluster table, augmented with per-sample output files
    """
    params:
        sce2g_config = config,
        e2g_config = encode_e2g_config,
        abc_config = encode_e2g.get_abc_config(encode_e2g_config),
        biosample_config = BIOSAMPLE_DF,
        out_dir = os.path.join(RESULTS_DIR, "config"), 
        results_dir = RESULTS_DIR,
        igv_dir = IGV_DIR
    output:
        sce2g_out = os.path.join(RESULTS_DIR, "config", "scE2G_config.yml"),
        e2g_out = os.path.join(RESULTS_DIR, "config", "ENCODE_rE2G_config.yml"),
        abc_out = os.path.join(RESULTS_DIR, "config", "ABC_config.yml"),
        res_out = os.path.join(RESULTS_DIR, "config", "expanded_biosample_config.tsv")
    run:
        import os, yaml, pandas as pd

        # make sure the output dir exists
        os.makedirs(params.out_dir, exist_ok=True)

        # 1) save sce2g config
        with open(output.sce2g_out, 'w') as fh:
            yaml.safe_dump(params.sce2g_config, fh)

        # 2) save encode_re2g config
        with open(output.e2g_out, 'w') as fh:
            yaml.safe_dump(params.e2g_config, fh)

        # 3) save abc config
        with open(output.abc_out, 'w') as fh:
            yaml.safe_dump(params.abc_config, fh)

        # 4) augment biosample table
        df = params.biosample_config

        # helper to construct file paths for each sample
        def mkpaths(row):
            biosample = row["biosample"]  # adjust to your column name
            model_name = row["model_dir_base"]
            threshold = row["model_threshold"]
            out = {}

            # predictions
            out["predictions_full"] = os.path.join(params.results_dir,
                biosample, model_name, "encode_e2g_predictions.tsv.gz")
            out["predictions_thresholded"] = os.path.join(params.results_dir,
                biosample, model_name, f"encode_e2g_predictions_threshold{threshold}.tsv.gz")
            
            # optional genome‐browser track
            if params.sce2g_config["make_IGV_tracks"]:
                out["predictions_bedpe"] = os.path.join(params.IGV_dir,
                    biosample, model_name, f"encode_e2g_predictions_threshold{threshold}.tsv.gz")
                out["ATAC_bw"] = os.path.join(params.IGV_dir,
                    biosample, "ATAC_norm.bw")

            return pd.Series(out)

        df = pd.concat([df, df.apply(mkpaths, axis=1)], axis=1)
        df.to_csv(output.res_out, sep='\t', index=False)

