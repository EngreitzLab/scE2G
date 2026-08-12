import glob
import logging
import os
import time
import unittest
from typing import Dict

import numpy as np
import pandas as pd
import yaml
from utils import get_biosample_names, get_filtered_dataframe, hash_large_file, run_cmd

logging.basicConfig(level=logging.INFO)

CONFIG_FILE = "tests/config/test_config.yml"
with open(CONFIG_FILE, "r") as file:
    CONFIG = yaml.safe_load(file)
COLUMNS_TO_COMPARE: Dict[str, type] = {
    "chr": str,
    "start": np.int64,
    "end": np.int64,
    "name": str,
    "class": str,
    "TargetGene": str,
    "ABC.Score": np.float64,
    "E2G.Score.qnorm": np.float64,
}
TEST_OUTPUT_DIR = CONFIG["results_dir"]
EXPECTED_OUTPUT_DIR = f"tests/expected_output/{CONFIG['TEST_CONFIG_NAME']}"
# ALL_PUTATIVE_PRED_FILE = "multiome_powerlaw_v2/scE2G_predictions.tsv.gz" # file too large to upload; just compare thresholded
# THRESHOLDED_PRED_FILE_PATTERN = (
#     "K562_chr22_cluster1/multiome_powerlaw_v2/scE2G_predictions_threshold*[0-9].tsv.gz"
# )
THRESHOLDED_PRED_FILE = "multiome_powerlaw_v3/scE2G_predictions_threshold0.177.tsv.gz"

# Kendall correlation output is small enough (~500KB) to compare as a full fixture.
KENDALL_FILE = "Kendall/Pairs.Kendall.tsv.gz"
KENDALL_COLUMNS_TO_COMPARE: Dict[str, type] = {
    "chr": str,
    "start": np.int64,
    "end": np.int64,
    "TargetGene": str,
    "Kendall": np.float64,
}

# ARC integration and pre-threshold predictions carry one row per candidate
# enhancer-gene pair genome-wide (~50-100MB) -- too large to git-commit as
# full fixtures, so these are checked via row count + a hash of identity
# columns and rounded score columns instead (see hash_large_file()).
IDENTITY_COLS = ["chr", "start", "end", "name", "class", "TargetGene"]
LARGE_FILE_CHECKS = {
    "ARC/EnhancerPredictionsAllPutative_ARC.tsv.gz": ["ABC.Score", "Kendall", "ARC.E2G.Score"],
    "multiome_powerlaw_v3/scE2G_predictions.tsv.gz": [
        "ABC.Score",
        "Kendall",
        "ARC.E2G.Score",
        "E2G.Score",
        "E2G.Score.qnorm",
    ],
}
PRE_THRESHOLD_PRED_FILE = "multiome_powerlaw_v3/scE2G_predictions.tsv.gz"

# On this test dataset, no distal (non-promoter) element has ever cleared the
# final threshold (confirmed pre-existing, not a regression), so the
# thresholded-file comparison never exercises distal scoring at all. Guard
# against distal scoring silently breaking (e.g. going all-zero/NaN) by
# checking the pre-threshold file directly instead.
DISTAL_CLASSES = ["genic", "intergenic"]

class scE2GTest(unittest.TestCase):
    def compare_prediction_file(self, biosample: str, pred_file, cols_to_compare=COLUMNS_TO_COMPARE) -> None:
        test_file = os.path.join(TEST_OUTPUT_DIR, biosample, pred_file)
        expected_file = os.path.join(EXPECTED_OUTPUT_DIR, biosample, pred_file)
        print(f"Comparing biosample: {biosample} for pred_file: {pred_file}")
        pd.testing.assert_frame_equal(
            get_filtered_dataframe(test_file, cols_to_compare),
            get_filtered_dataframe(expected_file, cols_to_compare),
        )

    def compare_large_file_checksum(self, biosample: str, pred_file: str, score_cols) -> None:
        test_file = os.path.join(TEST_OUTPUT_DIR, biosample, pred_file)
        expected_file = os.path.join(EXPECTED_OUTPUT_DIR, biosample, pred_file + ".rowcount_and_hash.txt")
        print(f"Comparing checksum for biosample: {biosample}, pred_file: {pred_file}")

        row_count, digest = hash_large_file(test_file, IDENTITY_COLS, score_cols)
        with open(expected_file) as f:
            expected_row_count, expected_digest = f.read().strip().split("\n")

        self.assertEqual(
            row_count,
            int(expected_row_count),
            msg=f"{pred_file}: row count changed ({row_count} vs expected {expected_row_count})",
        )
        self.assertEqual(
            digest,
            expected_digest,
            msg=f"{pred_file}: identity/rounded-score hash changed even though row count matches",
        )

    def check_distal_elements_have_nonzero_scores(self, biosample: str, pred_file: str) -> None:
        test_file = os.path.join(TEST_OUTPUT_DIR, biosample, pred_file)
        print(f"Checking distal element scores for biosample: {biosample}, pred_file: {pred_file}")
        df = pd.read_csv(test_file, sep="\t", usecols=["class", "E2G.Score.qnorm"])
        distal = df[df["class"].isin(DISTAL_CLASSES)]

        self.assertGreater(
            len(distal), 0, msg=f"{pred_file}: no distal ({'/'.join(DISTAL_CLASSES)}) candidate elements found at all"
        )
        self.assertGreater(
            distal["E2G.Score.qnorm"].max(),
            0,
            msg=f"{pred_file}: all distal element scores are exactly zero -- distal scoring may be broken",
        )
        self.assertGreater(
            distal["E2G.Score.qnorm"].std(),
            0,
            msg=f"{pred_file}: distal element scores have zero variance -- distal scoring may be degenerate",
        )

    # def compare_thresholded_prediction_file(self, biosample: str) -> None:
    #     test_files = glob.glob(
    #         os.path.join(TEST_OUTPUT_DIR, biosample, THRESHOLDED_PRED_FILE_PATTERN)
    #     )
    #     expected_files = glob.glob(
    #         os.path.join(EXPECTED_OUTPUT_DIR, biosample, THRESHOLDED_PRED_FILE_PATTERN)
    #     )
    #     if len(test_files) != 1:
    #         raise Exception(
    #             f"Multiple or no test thresholded files found. Please clean up. {test_files}"
    #         )
    #     if len(expected_files) != 1:
    #         raise Exception(
    #             f"Multiple or no expected thresholded files found. Please clean up. {expected_files}"
    #         )
    #     test_file = test_files[0]
    #     expected_file = expected_files[0]
    #     print(
    #         f"Comparing biosample: {biosample} for pred_file: {os.path.basename(test_file)}"
    #     )
    #     pd.testing.assert_frame_equal(
    #         get_filtered_dataframe(test_file, COLUMNS_TO_COMPARE),
    #         get_filtered_dataframe(expected_file, COLUMNS_TO_COMPARE),
    #     )

    def run_test(self, config_file: str) -> None:
        start = time.time()
        cmd = f"snakemake -j1 -F --configfile {config_file} --use-conda"
        run_cmd(cmd)
        time_taken = time.time() - start

        #biosample_names = get_biosample_names(CONFIG["cell_clusters"])
        biosample_names = ["K562_cluster1_chr22p"]
        for biosample in biosample_names:
            self.compare_prediction_file(biosample, KENDALL_FILE, KENDALL_COLUMNS_TO_COMPARE)
            for large_file, score_cols in LARGE_FILE_CHECKS.items():
                self.compare_large_file_checksum(biosample, large_file, score_cols)
            self.check_distal_elements_have_nonzero_scores(biosample, PRE_THRESHOLD_PRED_FILE)
            self.compare_prediction_file(biosample, THRESHOLDED_PRED_FILE)
            #self.compare_thresholded_prediction_file(biosample)

        # Make sure the test doesn't take too long
        # May need to adjust as more biosamples are added
        max_time = 60 * 50  # 50 min
        self.assertLessEqual(
            time_taken,
            max_time,
            msg=f"Running scE2G took too long: {int(time_taken/60)} minutes",
        )

    def test_sce2g_apply(self) -> None:
        self.run_test(CONFIG_FILE)


if __name__ == "__main__":
    unittest.main()
