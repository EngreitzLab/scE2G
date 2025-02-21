import pickle
import os
import click
import numpy as np
import pandas as pd
from scipy import interpolate

SCORE_COL_BASE = "E2G.Score"


def make_e2g_predictions(df_enhancers, feature_list, trained_model, epsilon):
    # transform the features
    X = df_enhancers.loc[:, feature_list]
    X = np.log(np.abs(X) + epsilon)

    with open(trained_model, "rb") as f:
        model = pickle.load(f)
    probs = model.predict_proba(X)

    df_enhancers[SCORE_COL_BASE] = probs[:, 1]

    return df_enhancers


def make_e2g_predictions_cv(df_enhancers, feature_list, cv_models, epsilon):
    score_col = SCORE_COL_BASE + ".cv"

    # get scores per chrom
    X = df_enhancers.loc[:, feature_list]
    X = np.log(np.abs(X) + epsilon)
    chr_list = np.unique(df_enhancers["chr"])

    for chr in chr_list:
        idx_test = df_enhancers[df_enhancers["chr"] == chr].index.values
        if len(idx_test) > 0:
            X_test = X.loc[idx_test, :]
            print(f"Length of X_test: {(X_test.shape)}")
            with open(os.path.join(cv_models, f"model_test_{chr}.pkl"), "rb") as f:
                model = pickle.load(f)
            probs = model.predict_proba(X_test)
            df_enhancers.loc[idx_test, score_col] = probs[:, 1]

    return df_enhancers


def calculate_quantiles(score_column):
    score_column_zero_replaced = score_column.replace(0, np.nan)
    ranks = score_column_zero_replaced.rank(
        method="average", na_option="top"
    )  # Calculate the rank, treating NaNs as the lowest possible ranks
    quantiles = (ranks - 1) / (len(score_column) - 1)  # Calculate the quantiles
    quantiles[score_column == 0] = (
        0  # Replace NaN ranks (original 0 values) with 0th quantile
    )
    return quantiles


def qnorm_scores(df_enhancers, qnorm_ref, crispr_benchmarking):
    # generate reference
    # ref_quantiles = np.linspace(0.01, 0.99, 99, endpoint=True) # [0.01, 0.02, ... , 0.99]
    # n_unique_scores = len(np.unique(qnorm_ref))
    ref_quantiles = qnorm_ref["quantile"].to_numpy()
    ref_scores = qnorm_ref["reference_score"].to_numpy()
    interpfunc = interpolate.interp1d(
        ref_quantiles, ref_scores, kind="linear", fill_value="extrapolate"
    )

    # qnorm regular score
    # score_quantile = df_enhancers[SCORE_COL_BASE + '.Score'].rank() / float(len(df_enhancers))
    score_quantile = calculate_quantiles(df_enhancers[SCORE_COL_BASE])
    df_enhancers[SCORE_COL_BASE + ".qnorm"] = interpfunc(score_quantile).clip(0)
    print(df_enhancers[SCORE_COL_BASE].describe())
    print(df_enhancers[SCORE_COL_BASE + ".qnorm"].describe())

    # qnorm cv score
    if str(crispr_benchmarking) == "True":
        cv_score_quantile = calculate_quantiles(df_enhancers[SCORE_COL_BASE + ".cv"])
        df_enhancers[SCORE_COL_BASE + ".cv.qnorm"] = interpfunc(cv_score_quantile).clip(
            0
        )

    return df_enhancers


def filter_by_tpm(df_enhancers, tpm_threshold, crispr_benchmarking):
    if ("RNA_pseudobulkTPM" not in df_enhancers.columns) or (tpm_threshold == 0):
        return df_enhancers  # don't filter

    df_enhancers[SCORE_COL_BASE + ".qnorm.ignoreTPM"] = df_enhancers[
        SCORE_COL_BASE + ".qnorm"
    ]
    df_enhancers[SCORE_COL_BASE + ".qnorm"] = [
        score if tpm >= tpm_threshold else 0
        for (score, tpm) in zip(
            df_enhancers[SCORE_COL_BASE + ".qnorm.ignoreTPM"],
            df_enhancers["RNA_pseudobulkTPM"],
        )
    ]

    if str(crispr_benchmarking) == "True":
        df_enhancers[SCORE_COL_BASE + ".cv.qnorm.ignoreTPM"] = df_enhancers[
            SCORE_COL_BASE + ".cv.qnorm"
        ]
        df_enhancers[SCORE_COL_BASE + ".cv.qnorm"] = [
            score if tpm >= tpm_threshold else 0
            for (score, tpm) in zip(
                df_enhancers[SCORE_COL_BASE + ".cv.qnorm.ignoreTPM"],
                df_enhancers["RNA_pseudobulkTPM"],
            )
        ]

    return df_enhancers


@click.command()
@click.option("--predictions", required=True)
@click.option("--feature_table_file", required=True)
@click.option("--trained_model", required=True)
@click.option("--model_dir", required=True)
@click.option("--crispr_benchmarking", required=True)
@click.option("--epsilon", type=float, default=0.01)
@click.option("--tpm_threshold", type=float, default=0)
@click.option("--output_file", required=True)
def main(
    predictions,
    feature_table_file,
    trained_model,
    model_dir,
    crispr_benchmarking,
    epsilon,
    tpm_threshold,
    output_file,
):
    feature_table = pd.read_csv(feature_table_file, sep="\t")
    feature_list = feature_table["feature"]

    # genomewide features
    df_enhancers = pd.read_csv(predictions, sep="\t")
    df_enhancers = df_enhancers.replace([np.inf, -np.inf], np.nan)
    df_enhancers = df_enhancers.fillna(0)

    # read in qnorm ref
    qnorm_file = os.path.join(
        model_dir, "qnorm_reference.tsv.gz"
    )  # col name = "E2G.Score.ignoreTPM"
    qnorm_ref = pd.read_csv(qnorm_file, sep="\t")

    df_enhancers = make_e2g_predictions(
        df_enhancers, feature_list, trained_model, epsilon
    )  # add "E2G.Score"

    if str(crispr_benchmarking) == "True":
        cv_models = os.path.join(model_dir, "cv_models")
        df_enhancers = make_e2g_predictions_cv(
            df_enhancers, feature_list, cv_models, epsilon
        )  # add "E2G.Score.cv"

    df_enhancers = qnorm_scores(
        df_enhancers, qnorm_ref, crispr_benchmarking
    )  # add "E2G.Score.qnorm" &  if crispr_benchmarking, "E2G.Score.cv.qnorm"
    df_enhancers = filter_by_tpm(
        df_enhancers, tpm_threshold, crispr_benchmarking
    )  # if required, add "E2G.Score.qnorm.ignoreTPM" &  if crispr_benchmarking, "E2G.Score.cv.qnorm.ignoreTPM" and filter the other columns

    df_enhancers.to_csv(output_file, compression="gzip", sep="\t", index=False)


if __name__ == "__main__":
    main()
