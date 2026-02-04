#!/bin/bash
# Manual test script to validate Singularity container workflow
# Run this on a system with Singularity installed (e.g., Stanford cluster)
#
# Usage:
#   conda activate run_snakemake
#   bash tests/test_singularity.sh

set -e

echo "=== Testing Singularity container workflow ==="

# Clean previous test output
rm -rf tests/test_output/powerlaw_v3_models/

# Run snakemake with container
echo "Running snakemake with --use-singularity..."
snakemake --configfile tests/config/test_config.yml \
    -j1 -F --use-conda --use-singularity

# Compare outputs
echo "Comparing outputs to expected..."
python -c "
import pandas as pd
import numpy as np

COLUMNS_TO_COMPARE = {
    'chr': str,
    'start': np.int64,
    'end': np.int64,
    'name': str,
    'class': str,
    'TargetGene': str,
    'ABC.Score': np.float64,
    'E2G.Score.qnorm': np.float64,
}

def get_filtered_dataframe(file_path, columns):
    df = pd.read_csv(file_path, sep='\t')
    df = df[list(columns.keys())]
    for col, dtype in columns.items():
        df[col] = df[col].astype(dtype)
    return df.sort_values(['chr', 'start', 'end', 'TargetGene']).reset_index(drop=True)

test_file = 'tests/test_output/powerlaw_v3_models/K562_cluster1_chr22p/multiome_powerlaw_v3/scE2G_predictions_threshold0.177.tsv.gz'
expected_file = 'tests/expected_output/powerlaw_v3_models/K562_cluster1_chr22p/multiome_powerlaw_v3/scE2G_predictions_threshold0.177.tsv.gz'

test_df = get_filtered_dataframe(test_file, COLUMNS_TO_COMPARE)
expected_df = get_filtered_dataframe(expected_file, COLUMNS_TO_COMPARE)

pd.testing.assert_frame_equal(test_df, expected_df)
print('Singularity container test PASSED!')
"

echo "=== Test complete ==="
