import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(
    0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "workflow", "scripts", "model_application")
)
from run_e2g_cv import calculate_quantiles, filter_by_tpm, qnorm_scores


def test_calculate_quantiles():
    scores = pd.Series([0, 0, 0.2, 0.5, 0.9])
    quantiles = calculate_quantiles(scores)
    assert quantiles.tolist() == pytest.approx([0.0, 0.0, 0.5, 0.75, 1.0])


def test_qnorm_scores_interpolates_against_reference():
    df = pd.DataFrame({"E2G.Score": [0, 0, 0.2, 0.5, 0.9]})
    qnorm_ref = pd.DataFrame(
        {"quantile": [0.0, 0.25, 0.5, 0.75, 1.0], "reference_score": [0.0, 0.1, 0.3, 0.6, 1.0]}
    )
    out = qnorm_scores(df.copy(), qnorm_ref, "False")
    assert out["E2G.Score.qnorm"].tolist() == pytest.approx([0.0, 0.0, 0.3, 0.6, 1.0])
    assert "E2G.Score.cv.qnorm" not in out.columns


def test_qnorm_scores_also_qnorms_cv_column_when_crispr_benchmarking():
    df = pd.DataFrame({"E2G.Score": [0, 0, 0.2, 0.5, 0.9], "E2G.Score.cv": [0, 0.1, 0.2, 0.5, 0.9]})
    qnorm_ref = pd.DataFrame(
        {"quantile": [0.0, 0.25, 0.5, 0.75, 1.0], "reference_score": [0.0, 0.1, 0.3, 0.6, 1.0]}
    )
    out = qnorm_scores(df.copy(), qnorm_ref, "True")
    assert out["E2G.Score.cv.qnorm"].tolist() == pytest.approx([0.0, 0.1, 0.3, 0.6, 1.0])


def test_filter_by_tpm_noop_when_rna_column_missing():
    df = pd.DataFrame({"E2G.Score.qnorm": [0.1, 0.5, 0.9]})
    out = filter_by_tpm(df.copy(), 1.0, "False")
    pd.testing.assert_frame_equal(out, df)


def test_filter_by_tpm_noop_when_threshold_zero():
    df = pd.DataFrame({"E2G.Score.qnorm": [0.1, 0.5, 0.9], "RNA_pseudobulkTPM": [0.5, 2.0, 5.0]})
    out = filter_by_tpm(df.copy(), 0, "False")
    pd.testing.assert_frame_equal(out, df)


def test_filter_by_tpm_zeros_scores_below_threshold():
    df = pd.DataFrame({"E2G.Score.qnorm": [0.1, 0.5, 0.9], "RNA_pseudobulkTPM": [0.5, 2.0, 5.0]})
    out = filter_by_tpm(df.copy(), 1.0, "False")
    assert out["E2G.Score.qnorm"].tolist() == pytest.approx([0.0, 0.5, 0.9])
    assert out["E2G.Score.qnorm.ignoreTPM"].tolist() == pytest.approx([0.1, 0.5, 0.9])
