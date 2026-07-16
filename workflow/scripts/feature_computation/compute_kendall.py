"""Compute Kendall's tau-b for candidate enhancer-gene pairs via fast_kendall_sc.

Standalone CLI invoked as a subprocess from compute_kendall.R. R still handles
all I/O (RDS/h5ad/csv.gz/10x reading), Seurat normalization, and GTF/TSS500
gene-name mapping -- none of that is the bottleneck. This script only does the
actual per-pair correlation computation, which previously took ~7 hours in
pure R at ~10M-pair scale.

Inputs (all written by R): the (already gene-name-mapped, normalized) RNA
matrix and the (already binarized) ATAC matrix in MatrixMarket format
(genes/peaks x cells), a plain-text list of row names for each, and a
tab-separated (TargetGene, PeakName) pairs file with no header, one row per
candidate pair, in the exact output order desired.

Output: one Kendall tau-b value per line (or "NA" where undefined), in the
same order as the input pairs file, for R to read back positionally.
"""

import argparse
import os


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rna-mtx", required=True, help="RNA matrix, MatrixMarket format, genes x cells")
    parser.add_argument("--rna-genes", required=True, help="One gene name per line, matching --rna-mtx row order")
    parser.add_argument("--atac-mtx", required=True, help="ATAC matrix, MatrixMarket format, peaks x cells, binarized")
    parser.add_argument("--atac-peaks", required=True, help="One peak name per line, matching --atac-mtx row order")
    parser.add_argument("--pairs", required=True, help="TSV, no header, columns: TargetGene, PeakName")
    parser.add_argument("--out", required=True, help="Output path: one Kendall value per line, same order as --pairs")
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()

    # rayon's global thread pool is sized once, on first use -- must be set
    # before fast_kendall_sc (and its Rust extension) is imported.
    os.environ.setdefault("RAYON_NUM_THREADS", str(max(1, args.threads)))

    import numpy as np
    import scipy.io
    from fast_kendall_sc import batch_kendall_tau

    rna_matrix = scipy.io.mmread(args.rna_mtx).T.tocsc()
    atac_matrix = scipy.io.mmread(args.atac_mtx).T.tocsc()

    with open(args.rna_genes) as f:
        gene_names = [line.rstrip("\n") for line in f]
    with open(args.atac_peaks) as f:
        peak_names = [line.rstrip("\n") for line in f]

    gene_to_index = {name: i for i, name in enumerate(gene_names)}
    peak_to_index = {name: i for i, name in enumerate(peak_names)}

    pair_gene_indices = []
    pair_peak_indices = []
    with open(args.pairs) as f:
        for line in f:
            gene_name, peak_name = line.rstrip("\n").split("\t")
            pair_gene_indices.append(gene_to_index[gene_name])
            pair_peak_indices.append(peak_to_index[peak_name])

    gene_indices = np.array(pair_gene_indices, dtype=np.int64)
    peak_indices = np.array(pair_peak_indices, dtype=np.int64)

    results = batch_kendall_tau(rna_matrix, atac_matrix, gene_indices, peak_indices)

    with open(args.out, "w") as f:
        for value in results:
            f.write("NA\n" if np.isnan(value) else f"{value}\n")


if __name__ == "__main__":
    main()
