import os
import sys

import numpy as np
import pytest

sys.path.insert(
    0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "workflow", "scripts", "prediction_qc")
)
from training_functions import (
    auc_mod,
    statistic_aupr,
    statistic_precision_at_threshold,
    statistic_recall_at_threshold,
    threshold_at_target_recall,
)


def test_auc_mod_single_point_returns_zero():
    assert auc_mod(np.array([1.0]), np.array([0.5])) == 0


def test_auc_mod_two_points():
    assert auc_mod(np.array([0.0, 1.0]), np.array([1.0, 1.0])) == pytest.approx(1.0)


def test_statistic_aupr_perfect_separation():
    y_true = np.array([0, 0, 1, 1])
    y_pred = np.array([0.1, 0.2, 0.8, 0.9])
    assert statistic_aupr(y_true, y_pred) == pytest.approx(1.0)


def test_statistic_aupr_imperfect_separation():
    y_true = np.array([0, 1, 0, 1, 1])
    y_pred = np.array([0.2, 0.4, 0.5, 0.6, 0.9])
    assert statistic_aupr(y_true, y_pred) == pytest.approx(0.9027777777777777)


def test_threshold_at_target_recall_achievable():
    y_true = np.array([0, 1, 0, 1, 1])
    y_pred = np.array([0.2, 0.4, 0.5, 0.6, 0.9])
    assert threshold_at_target_recall(y_true, y_pred, 0.5) == pytest.approx(0.5)


def test_threshold_at_target_recall_unachievable_returns_none():
    y_true = np.array([0, 1, 0, 1, 1])
    y_pred = np.array([0.2, 0.4, 0.5, 0.6, 0.9])
    assert threshold_at_target_recall(y_true, y_pred, 1.5) is None


def test_statistic_precision_at_threshold_none_returns_zero():
    y_true = np.array([0, 1, 0, 1, 1])
    y_pred = np.array([0.2, 0.4, 0.5, 0.6, 0.9])
    assert statistic_precision_at_threshold(y_true, y_pred, None) == 0


def test_statistic_precision_at_threshold():
    y_true = np.array([0, 1, 0, 1, 1])
    y_pred = np.array([0.2, 0.4, 0.5, 0.6, 0.9])
    assert statistic_precision_at_threshold(y_true, y_pred, 0.5) == pytest.approx(0.6666666666666666)


def test_statistic_recall_at_threshold_beyond_max_returns_zero():
    y_true = np.array([0, 1, 0, 1, 1])
    y_pred = np.array([0.2, 0.4, 0.5, 0.6, 0.9])
    assert statistic_recall_at_threshold(y_true, y_pred, 10.0) == 0
