#!/bin/bash
# Used when we want to change the expected output to match the latest test output run
set -e

if [ $# -ne 1 ]; then
    echo "Usage: $0 <results_dir>"
    exit 1
fi
results_dir=$1
biosample=K562_cluster1_chr22p

# remove old expected output
rm -rf expected_output/powerlaw_v3_models/*

# small files: copy in full
mkdir -p expected_output/powerlaw_v3_models/${biosample}/multiome_powerlaw_v3/
mkdir -p expected_output/powerlaw_v3_models/${biosample}/Kendall/
cp ${results_dir}/${biosample}/multiome_powerlaw_v3/scE2G_predictions_threshold0.177_stats.tsv expected_output/powerlaw_v3_models/${biosample}/multiome_powerlaw_v3/
cp ${results_dir}/${biosample}/multiome_powerlaw_v3/scE2G_predictions_threshold0.177.tsv.gz expected_output/powerlaw_v3_models/${biosample}/multiome_powerlaw_v3/
cp ${results_dir}/${biosample}/Kendall/Pairs.Kendall.tsv.gz expected_output/powerlaw_v3_models/${biosample}/Kendall/

# large files (ARC output, pre-threshold predictions): too big to commit in
# full (one row per candidate enhancer-gene pair genome-wide), so store just
# a row count + identity/rounded-score hash instead. See
# tests/utils.py:hash_large_file() and tests/test_sce2g_apply.py.
mkdir -p expected_output/powerlaw_v3_models/${biosample}/ARC/
python3 - "$results_dir" "$biosample" <<'EOF'
import sys
sys.path.insert(0, ".")
from utils import hash_large_file

results_dir, biosample = sys.argv[1], sys.argv[2]
checks = {
    "ARC/EnhancerPredictionsAllPutative_ARC.tsv.gz": ["ABC.Score", "Kendall", "ARC.E2G.Score"],
    "multiome_powerlaw_v3/scE2G_predictions.tsv.gz": [
        "ABC.Score", "Kendall", "ARC.E2G.Score", "E2G.Score", "E2G.Score.qnorm",
    ],
}
identity_cols = ["chr", "start", "end", "name", "class", "TargetGene"]

for rel_path, score_cols in checks.items():
    row_count, digest = hash_large_file(
        f"{results_dir}/{biosample}/{rel_path}", identity_cols, score_cols
    )
    out_path = f"expected_output/powerlaw_v3_models/{biosample}/{rel_path}.rowcount_and_hash.txt"
    with open(out_path, "w") as f:
        f.write(f"{row_count}\n{digest}\n")
    print(f"wrote {out_path}: {row_count} rows, {digest}")
EOF
