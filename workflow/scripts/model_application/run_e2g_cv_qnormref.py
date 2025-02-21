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
	score_col =  SCORE_COL_BASE + ".cv"

	# get scores per chrom
	X = df_enhancers.loc[:, feature_list]
	X = np.log(np.abs(X) + epsilon)
	chr_list = np.unique(df_enhancers["chr"])

	for chr in chr_list:
		idx_test = df_enhancers[df_enhancers["chr"] == chr].index.values
		if len(idx_test) > 0:
			X_test = X.loc[idx_test, :]
			with open(os.path.join(cv_models, f"model_test_{chr}.pkl"), "rb") as f:
				model = pickle.load(f)
			probs = model.predict_proba(X_test)
			df_enhancers.loc[idx_test, score_col] = probs[:, 1]

	return df_enhancers

def filter_by_tpm(df_enhancers, tpm_threshold, crispr_benchmarking):
	if ("RNA_pseudobulkTPM" not in df_enhancers.columns) or (tpm_threshold == 0):
		return df_enhancers  # don't filter

	df_enhancers[SCORE_COL_BASE + ".ignoreTPM"] = df_enhancers[
		SCORE_COL_BASE 
	]
	df_enhancers[SCORE_COL_BASE ] = [
		score if tpm >= tpm_threshold else 0
		for (score, tpm) in zip(
			df_enhancers[SCORE_COL_BASE + ".ignoreTPM"],
			df_enhancers["RNA_pseudobulkTPM"],
		)
	]

	if str(crispr_benchmarking) == "True":
		df_enhancers[SCORE_COL_BASE + ".cv.ignoreTPM"] = df_enhancers[
			SCORE_COL_BASE + ".cv"
		]
		df_enhancers[SCORE_COL_BASE + ".cv"] = [
			score if tpm >= tpm_threshold else 0
			for (score, tpm) in zip(
				df_enhancers[SCORE_COL_BASE + ".cv.ignoreTPM"],
				df_enhancers["RNA_pseudobulkTPM"],
			)
		]

	return df_enhancers

def get_qnorm_quantiles(input_df, n_scores):
	all_scores = input_df["E2G.Score"]
	n_unique_scores = int(n_scores)

	ref_quantiles = np.linspace(
		0, 1, n_unique_scores, endpoint=True
	)  # [0, 0.01, 0.02, ... , 0.99, 1]
	ref_scores = np.quantile(all_scores, ref_quantiles)

	out_df = pd.DataFrame({"quantile": ref_quantiles, "reference_score": ref_scores})

	return out_df


@click.command()
@click.option("--predictions", required=True)
@click.option("--feature_table_file", required=True)
@click.option("--trained_model", required=True)
@click.option("--model_dir", required=True)
@click.option("--epsilon", type=float, default=0.01)
@click.option("--tpm_threshold", type=float, default=0)
@click.option("--n_scores", type=float, default=0)
@click.option("--output_file_pred", required=True)
@click.option("--output_file_qnorm", required=True)

## modified functions from run_e2g_cv.py to create qnorm reference and full-depth predictions to benchmark to find threshold; do not use for other purposes
def main(
	predictions,
	feature_table_file,
	trained_model,
	model_dir,
	epsilon,
	tpm_threshold,
	n_scores,
	output_file_pred,
	output_file_qnorm
):
	feature_table = pd.read_csv(feature_table_file, sep="\t")
	feature_list = feature_table["feature"]

	# genomewide features
	df_enhancers = pd.read_csv(predictions, sep="\t")
	df_enhancers = df_enhancers.replace([np.inf, -np.inf], np.nan)
	df_enhancers = df_enhancers.fillna(0)

	# make non-cv predictions
	df_enhancers = make_e2g_predictions(
		df_enhancers, feature_list, trained_model, epsilon
	)  # add "E2G.Score"

	# make cv predictions
	cv_models = os.path.join(model_dir, "cv_models")
	df_enhancers = make_e2g_predictions_cv(
			df_enhancers, feature_list, cv_models, epsilon
		)  # add "E2G.Score.cv"
	
	# save qnorm reference before TPM filter
	out_df = get_qnorm_quantiles(df_enhancers, n_scores)
	out_df.to_csv(output_file_qnorm, compression="gzip", sep="\t", index=False)

	# TPM filter (for benchmarking)
	df_enhancers = filter_by_tpm(
		df_enhancers, tpm_threshold, "True"
	) 

	df_enhancers.to_csv(output_file_pred, compression="gzip", sep="\t", index=False)

if __name__ == "__main__":
	main()

# python workflow/scripts/model_application/run_e2g_cv_qnormref.py \
# 	--predictions results/2025_0206_new_gene_universe/K562_Xu_pl/genomewide_features.tsv.gz \
# 	--feature_table_file models/multiome_powerlaw_v3/feature_table.tsv \
# 	--trained_model models/multiome_powerlaw_v3/model.pkl \
# 	--model_dir models/multiome_powerlaw_v3 \
# 	--epsilon 0.01 \
# 	--tpm_threshold 1 \
# 	--n_scores 10000 \
# 	--output_file_pred results/2025_0206_new_gene_universe/K562_Xu_pl/multiome_powerlaw_v3/genomewide_predictions.tsv.gz \
# 	--output_file_qnorm models/multiome_powerlaw_v3/qnorm_reference.tsv.gz

# python workflow/scripts/model_application/run_e2g_cv_qnormref.py \
# 	--predictions results/2025_0206_new_gene_universe/K562_Xu_pl/genomewide_features.tsv.gz \
# 	--feature_table_file models/scATAC_powerlaw_v3/feature_table.tsv \
# 	--trained_model models/scATAC_powerlaw_v3/model.pkl \
# 	--model_dir models/scATAC_powerlaw_v3 \
# 	--epsilon 0.01 \
# 	--tpm_threshold 0 \
# 	--n_scores 10000 \
# 	--output_file_pred results/2025_0206_new_gene_universe/K562_Xu_pl/scATAC_powerlaw_v3/genomewide_predictions.tsv.gz \
# 	--output_file_qnorm models/scATAC_powerlaw_v3/qnorm_reference.tsv.gz