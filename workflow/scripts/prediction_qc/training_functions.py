import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve, auc
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LogisticRegression


## statistic functions for delta auPR/precision to be used for scipy.stats.bootstrap


# modify auc to handle single positive
def auc_mod(recall, precision):
    if len(recall) < 2:
        return 0
    else:
        return auc(recall, precision)


# the first precision & recall values are precision=class_balance and recall=1, not corresponding to a threshold. we do not want this point to affect our performance.
def precision_recall_curve_modified(y_true, y_pred):
    precision_full, recall_full, thresholds_full = precision_recall_curve(
        y_true, y_pred
    )
    return precision_full[1:], recall_full[1:], thresholds_full


def statistic_delta_aupr(
    y_true, y_pred_full, y_pred_ablated
):  # return aupr_ablated-aupr_full
    precision_full, recall_full, thresholds_full = precision_recall_curve_modified(
        y_true, y_pred_full
    )
    aupr_full = auc_mod(recall_full, precision_full)

    precision_ablated, recall_ablated, thresholds_ablated = (
        precision_recall_curve_modified(y_true, y_pred_ablated)
    )
    aupr_ablated = auc_mod(recall_ablated, precision_ablated)
    delta_aupr = aupr_ablated - aupr_full

    return delta_aupr


def statistic_aupr(y_true, y_pred_full):
    precision_full, recall_full, thresholds_full = precision_recall_curve_modified(
        y_true, y_pred_full
    )
    aupr_full = auc_mod(recall_full, precision_full)

    return aupr_full


# note: precision at 70% recall (vs precision at constant threshold chosen for 70% recall)
def statistic_delta_precision(
    y_true, y_pred_full, y_pred_ablated
):  # return precision_ablated-precision_full
    precision_full, recall_full, thresholds_full = precision_recall_curve_modified(
        y_true, y_pred_full
    )
    idx_recall_full_70_pct = np.argsort(np.abs(recall_full - 0.7))[0]
    precision_full_at_70_pct_recall = precision_full[idx_recall_full_70_pct]

    precision_ablated, recall_ablated, thresholds_ablated = (
        precision_recall_curve_modified(y_true, y_pred_ablated)
    )
    idx_recall_ablated_70_pct = np.argsort(np.abs(recall_ablated - 0.7))[0]
    precision_ablated_at_70_pct_recall = precision_ablated[idx_recall_ablated_70_pct]
    delta_precision = (
        precision_ablated_at_70_pct_recall - precision_full_at_70_pct_recall
    )

    return delta_precision


def statistic_delta_precision_at_threshold(
    y_true, y_pred_full, y_pred_ablated, thresh_full, thresh_ablated
):  # return precision_ablated-precision_full
    precision_full, recall_full, thresholds_full = precision_recall_curve_modified(
        y_true, y_pred_full
    )
    if thresh_full != None:
        idx_threshold = np.argsort(np.abs(thresholds_full - thresh_full))[0]
        precision_at_threshold = precision_full[idx_threshold - 1]
    else:
        precision_at_threshold = 0

    precision_ablated, recall_ablated, thresholds_ablated = (
        precision_recall_curve_modified(y_true, y_pred_ablated)
    )
    if thresh_ablated != None:
        idx_threshold_ablated = np.argsort(np.abs(thresholds_ablated - thresh_ablated))[
            0
        ]
        precision_ablated_at_threshold = precision_ablated[idx_threshold_ablated - 1]
    else:
        precision_ablated_at_threshold = 0

    delta_precision = precision_ablated_at_threshold - precision_at_threshold

    return delta_precision


# note: precision at 70% recall (vs precision at constant threshold chosen for 70% recall)
def statistic_precision(y_true, y_pred_full):
    precision_full, recall_full, thresholds_full = precision_recall_curve_modified(
        y_true, y_pred_full
    )
    idx_recall_full_70_pct = np.argsort(np.abs(recall_full - 0.7))[0]
    precision_full_at_70_pct_recall = precision_full[idx_recall_full_70_pct]
    return precision_full_at_70_pct_recall


# return threshold for 70% recall
def threshold_70_pct_recall(y_true, y_pred_full):
    precision_full, recall_full, thresholds_full = precision_recall_curve_modified(
        y_true, y_pred_full
    )
    if np.max(recall_full) > 0.7:
        idx_recall_full_70_pct = np.argsort(np.abs(recall_full - 0.7))[
            0
        ]  # find index of recall closest to 0.7
        thresh = thresholds_full[
            idx_recall_full_70_pct + 1
        ]  # adjust for removing first row of rec + prec
    else:
        thresh = None
    return thresh


def threshold_at_target_recall(y_true, y_pred_full, recall_target):
    precision_full, recall_full, thresholds_full = precision_recall_curve_modified(
        y_true, y_pred_full
    )
    if np.max(recall_full) > recall_target:
        idx_recall_at_target = np.argsort(np.abs(recall_full - recall_target))[
            0
        ]  # find index of recall closest to 0.7
        thresh = thresholds_full[
            idx_recall_at_target + 1
        ]  # adjust for removing first row of rec + prec
    else:
        thresh = None
    return thresh


# note: precision at constant threshold chosen for 70% recall (vs precision at 70% recall)
def statistic_precision_at_threshold(y_true, y_pred_full, threshold):
    if threshold != None:
        precision_full, recall_full, thresholds_full = precision_recall_curve_modified(
            y_true, y_pred_full
        )
        idx_threshold = np.argsort(np.abs(thresholds_full - threshold))[
            0
        ]  # index of thresholds_full for value closest threshold
        precision_at_threshold = precision_full[
            idx_threshold - 1
        ]  # precision at corresponding index
    else:
        precision_at_threshold = 0

    return precision_at_threshold


def statistic_recall_at_threshold(y_true, y_pred_full, threshold):
    precision_full, recall_full, thresholds_full = precision_recall_curve_modified(
        y_true, y_pred_full
    )
    if threshold < np.max(thresholds_full):
        idx_threshold = np.argsort(np.abs(thresholds_full - threshold))[
            0
        ]  # index of thresholds_full for value closest threshold
        recall_at_threshold = recall_full[idx_threshold - 1]
    else:
        recall_at_threshold = 0
    return recall_at_threshold


# bootstrap p-values for delta (aupr/precision)
def bootstrap_pvalue(delta, res_delta):
    # Generate boostrap distribution of delta under null hypothesis (important centering step to get sampling distribution under the null)
    delta_boot_distribution = (
        res_delta.bootstrap_distribution - res_delta.bootstrap_distribution.mean()
    )

    # Calculate proportion of bootstrap samples with at least as strong evidence against null
    pval = np.mean(np.abs(delta_boot_distribution) >= np.abs(delta))

    return pval


# assumes necessary features are present in X and Y, features are already transformed
def train_and_predict_once(df_dataset, X, Y, feature_list, model_name, params):
    X = X.loc[:, feature_list]
    idx = np.arange(len(Y))  # number of elements

    chr_list = np.unique(df_dataset["chr"])
    if len(chr_list) > 1:
        for chr in chr_list:
            idx_test = df_dataset[df_dataset["chr"] == chr].index.values
            if len(idx_test) > 0:
                idx_train = np.delete(idx, idx_test)
                X_test = X.loc[idx_test, :]
                X_train = X.loc[idx_train, :]
                Y_train = Y[idx_train]
                model = LogisticRegression(**params).fit(X_train, Y_train)
                probs = model.predict_proba(X_test)  # calculate scores
                df_dataset.loc[idx_test, model_name + ".Score"] = probs[:, 1]
    else:
        model = LogisticRegression(**params).fit(X, Y)
        probs = model.predict_proba(X)  # calculate scores
        df_dataset.loc[:, model_name + ".Score"] = probs[:, 1]

    return df_dataset
