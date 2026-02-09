import pickle

import click
import numpy as np
import pandas as pd

SCORE_COLUMN_BASE = "E2G.Score"


def make_e2g_predictions(df_enhancers, feature_list, trained_model, tpm_threshold, epsilon):
    # transform the features
    X = df_enhancers.loc[:, feature_list]
    X = np.log(np.abs(X) + epsilon)

    with open(trained_model, "rb") as f:
        model = pickle.load(f)
    probs = model.predict_proba(X)
    
    if ("RNA_pseudobulkTPM" in df_enhancers.columns) and (tpm_threshold > 0):
        df_enhancers[SCORE_COLUMN_BASE + ".ignoreTPM"] = probs[:, 1]
        df_enhancers[SCORE_COLUMN_BASE] = [score if tpm >= tpm_threshold else 0 
            for (score, tpm) in zip(df_enhancers[SCORE_COLUMN_BASE + ".ignoreTPM"] , df_enhancers["RNA_pseudobulkTPM"])]
    else:
         df_enhancers[SCORE_COLUMN_BASE] = probs[:, 1]

    return df_enhancers


@click.command()
@click.option("--predictions", required=True)
@click.option("--feature_table_file", required=True)
@click.option("--trained_model", required=True)
@click.option("--tpm_threshold", default=0)
@click.option("--epsilon", type=float, default=0.01)
@click.option("--output_file", required=True)
def main(predictions, feature_table_file, trained_model, tpm_threshold, epsilon, output_file):
    feature_table = pd.read_csv(feature_table_file, sep="\t")
    feature_list = feature_table["feature"]

    df_enhancers = pd.read_csv(predictions, sep="\t")
    df_enhancers = df_enhancers.replace([np.inf, -np.inf], np.nan)
    df_enhancers = df_enhancers.fillna(0)

    df_enhancers = make_e2g_predictions(
        df_enhancers, feature_list, trained_model, tpm_threshold, epsilon
    )
    df_enhancers.to_csv(output_file, compression="gzip", sep="\t", index=False)


if __name__ == "__main__":
    main()
