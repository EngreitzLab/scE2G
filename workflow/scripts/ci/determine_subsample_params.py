"""
Determine adaptive subsampling parameters for a cell cluster.

Counts unique barcodes and total fragments in a fragment file, then computes
the actual subsample fraction that ensures minimum cell and fragment thresholds
are met while targeting a configurable fraction.
"""

import gzip
import json
import os
import sys
import click


@click.command()
@click.option("--frag_file", required=True, help="Path to bgzipped fragment file (5-column TSV)")
@click.option("--target_fraction", type=float, default=0.70, help="Target fraction of cells per subsample")
@click.option("--min_cells", type=int, default=100, help="Minimum cells per subsample")
@click.option("--min_fragments", type=int, default=2000000, help="Minimum fragments per subsample")
@click.option("--output", required=True, help="Output JSON file path")
def main(frag_file, target_fraction, min_cells, min_fragments, output):
    barcodes = set()
    n_fragments = 0

    open_fn = gzip.open if frag_file.endswith(".gz") else open
    with open_fn(frag_file, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            barcodes.add(fields[3])
            n_fragments += 1

    n_cells = len(barcodes)

    # Compute minimum fraction needed to satisfy thresholds
    min_fraction_for_cells = min_cells / n_cells if n_cells > 0 else 1.0
    min_fraction_for_frags = min_fragments / n_fragments if n_fragments > 0 else 1.0
    min_required = max(min_fraction_for_cells, min_fraction_for_frags)

    actual_fraction = min(1.0, max(target_fraction, min_required))

    n_cells_per_subsample = round(n_cells * actual_fraction)
    n_fragments_estimated = round(n_fragments * actual_fraction)

    if actual_fraction >= 0.95 and n_cells < min_cells:
        print(
            f"WARNING: cluster has only {n_cells} cells total; even 100% subsampling "
            f"does not reach min_cells={min_cells}. CIs may be uninformative.",
            file=sys.stderr,
        )
    elif actual_fraction >= 0.95:
        print(
            f"WARNING: actual_fraction={actual_fraction:.2f} (cluster is small); "
            "subsamples will be nearly identical and CIs will be tight.",
            file=sys.stderr,
        )

    params = {
        "n_cells": n_cells,
        "n_fragments": n_fragments,
        "target_fraction": target_fraction,
        "actual_fraction": actual_fraction,
        "n_cells_per_subsample": n_cells_per_subsample,
        "n_fragments_estimated": n_fragments_estimated,
        "min_cells": min_cells,
        "min_fragments": min_fragments,
    }

    out_dir = os.path.dirname(os.path.abspath(output))
    os.makedirs(out_dir, exist_ok=True)
    with open(output, "w") as fh:
        json.dump(params, fh, indent=2)

    print(
        f"Cluster: {n_cells} cells, {n_fragments} fragments. "
        f"Subsample fraction: {actual_fraction:.3f} "
        f"({n_cells_per_subsample} cells, ~{n_fragments_estimated} fragments per subsample)"
    )


if __name__ == "__main__":
    main()
