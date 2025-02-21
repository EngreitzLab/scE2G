import glob
import logging
import os
import time
import unittest
from typing import Dict

import numpy as np
import pandas as pd
import yaml
from utils import get_biosample_names, get_filtered_dataframe, run_cmd

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
# ALL_PUTATIVE_PRED_FILE = "multiome_powerlaw_v2/encode_e2g_predictions.tsv.gz" # file too large to upload; just compare thresholded
# THRESHOLDED_PRED_FILE_PATTERN = (
#     "K562_chr22_cluster1/multiome_powerlaw_v2/encode_e2g_predictions_threshold*[0-9].tsv.gz"
# )
THRESHOLDED_PRED_FILE = "multiome_powerlaw_v3/encode_e2g_predictions_threshold0.149.tsv.gz"

class scE2GTest(unittest.TestCase):
    def compare_prediction_file(self, biosample: str, pred_file) -> None:
        test_file = os.path.join(TEST_OUTPUT_DIR, biosample, pred_file)
        expected_file = os.path.join(EXPECTED_OUTPUT_DIR, biosample, pred_file)
        print(f"Comparing biosample: {biosample} for pred_file: {pred_file}")
        pd.testing.assert_frame_equal(
            get_filtered_dataframe(test_file, COLUMNS_TO_COMPARE),
            get_filtered_dataframe(expected_file, COLUMNS_TO_COMPARE),
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
            self.compare_prediction_file(biosample, THRESHOLDED_PRED_FILE)
            #self.compare_thresholded_prediction_file(biosample)

        # Make sure the test doesn't take too long
        # May need to adjust as more biosamples are added, but we should keep
        # tests quick, so don't run rE2G on all chromosomes
        max_time = 60 * 30  # 30 min
        self.assertLessEqual(
            time_taken,
            max_time,
            msg=f"Running scE2G took too long: {int(time_taken/60)} minutes",
        )

    def test_sce2g_apply(self) -> None:
        self.run_test(CONFIG_FILE)


if __name__ == "__main__":
    unittest.main()
