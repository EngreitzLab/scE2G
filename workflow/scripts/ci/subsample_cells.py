"""
Generate a subsampled fragment file by randomly selecting a fraction of cell barcodes.

Reads subsample parameters from the JSON produced by determine_subsample_params.py,
selects barcodes using a deterministic seed (base_seed + subsample_index), filters
the fragment file to those barcodes, and outputs a bgzipped + tabix-indexed file.
"""

import gzip
import json
import os
import random
import subprocess
import sys
import tempfile
import click


@click.command()
@click.option("--frag_file", required=True, help="Path to bgzipped fragment file")
@click.option("--params_json", required=True, help="Path to subsample_params.json")
@click.option("--subsample_index", type=int, required=True, help="Index of this subsample (0-based)")
@click.option("--base_seed", type=int, default=42, help="Base random seed; actual seed = base_seed + subsample_index")
@click.option("--output", required=True, help="Output bgzipped fragment file path (.tsv.gz)")
def main(frag_file, params_json, subsample_index, base_seed, output):
    with open(params_json) as fh:
        params = json.load(fh)

    n_cells_to_select = params["n_cells_per_subsample"]
    seed = base_seed + subsample_index

    # First pass: collect all unique barcodes
    barcodes = set()
    open_fn = gzip.open if frag_file.endswith(".gz") else open
    with open_fn(frag_file, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                barcodes.add(line.rstrip("\n").split("\t")[3])

    unique_barcodes = list(barcodes)
    n_total = len(unique_barcodes)

    if n_cells_to_select > n_total:
        print(
            f"WARNING: requested {n_cells_to_select} cells but only {n_total} available; using all.",
            file=sys.stderr,
        )
        n_cells_to_select = n_total

    rng = random.Random(seed)
    selected = set(rng.sample(unique_barcodes, n_cells_to_select))

    print(
        f"Subsample {subsample_index}: selected {len(selected)} of {n_total} cells (seed={seed})",
        file=sys.stderr,
    )

    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)

    # Second pass: filter fragments to selected barcodes, pipe to bgzip
    # Write to a temp uncompressed file first, then bgzip + tabix
    tmp_fd, tmp_tsv = tempfile.mkstemp(suffix=".tsv", dir=os.path.dirname(os.path.abspath(output)))
    os.close(tmp_fd)
    try:
        with open_fn(frag_file, "rt") as fh_in, open(tmp_tsv, "w") as fh_out:
            for line in fh_in:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if fields[3] in selected:
                    fh_out.write(line)

        # bgzip compress
        subprocess.run(["bgzip", "-f", tmp_tsv], check=True)
        bgzipped = tmp_tsv + ".gz"
        os.rename(bgzipped, output)

        # tabix index
        subprocess.run(["tabix", "-p", "bed", output], check=True)

    finally:
        if os.path.exists(tmp_tsv):
            os.remove(tmp_tsv)

    # Count output fragments for logging
    n_out = 0
    with gzip.open(output, "rt") as fh:
        for _ in fh:
            n_out += 1
    print(f"Subsample {subsample_index}: wrote {n_out} fragments to {output}", file=sys.stderr)


if __name__ == "__main__":
    main()
