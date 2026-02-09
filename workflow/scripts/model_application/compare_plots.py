import glob
import os

import click
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from generate_plots import VALUE_KEY
from matplotlib.backends.backend_pdf import PdfPages

"""
py workflow/scripts/model_application/compare_plots.py --v1_dir /oak/stanford/groups/engreitz/Users/atan5133/ENCODE_rE2G/results/ENCODE_v1 --v2_dhs_dir /oak/stanford/groups/engreitz/Users/atan5133/encode_dataset_processing/results/dhs_only --v2_dhs_hic_dir /oak/stanford/groups/engreitz/Users/atan5133/encode_dataset_processing/results/dhs_hic --output_file comparison_plots.pdf
"""


def get_v1_dhs_encode_ids(dir):
    files = os.listdir(dir)
    encode_ids = []
    for file in files:
        exp_id = file.split("_")[-1]
        encode_ids.append(exp_id)
    return set(encode_ids)


def get_matching_experiments(v1_encode_ids, dir):
    """
    Returns list of matching (exp_id, encode_id)
    """
    matching = set()
    for cluster in os.listdir(dir):
        if "ENC" not in cluster:
            continue
        exp_id = [
            encode_id for encode_id in cluster.split("_") if encode_id.startswith("ENC")
        ][0]
        files_with_encode_id = glob.glob(os.path.join(dir, cluster, "*", "*ENC*"))
        for file in files_with_encode_id:
            file = os.path.basename(file)
            encode_id = file[file.index("ENC") :].split(".")[0].split("_")[0]
            if encode_id in v1_encode_ids:
                matching.add((exp_id, encode_id))
    return matching


def get_num_eg_links(id, dir):
    file = glob.glob(os.path.join(dir, f"*{id}*", "Metrics", "*_stats.tsv"))[0]
    df = pd.read_csv(file, sep="\t").set_index("Metric")
    return df.loc["num_enh_gene_links", VALUE_KEY]


def plot_enh_gene_links(v1_dir, v2_dir, matching_exp_ids, title):
    plt.clf()
    v1 = []
    v2 = []
    for exp_id, encode_id in matching_exp_ids:
        v1.append(get_num_eg_links(encode_id, v1_dir))
        v2.append(get_num_eg_links(exp_id, v2_dir))

    ax = sns.scatterplot(x=v1, y=v2, label=f"n={len(v1)}")
    plt.axis("equal")
    ax.set_xlabel("V1")
    ax.set_ylabel("V2")
    ax.set_title(title)
    plt.legend()
    return ax.get_figure()


@click.command()
@click.option("--v1_dir")
@click.option("--v2_dhs_dir")
@click.option("--v2_dhs_hic_dir")
@click.option("--output_file")
def main(v1_dir, v2_dhs_dir, v2_dhs_hic_dir, output_file):
    v1_encode_ids = get_v1_dhs_encode_ids(v1_dir)
    v2_dhs_matching_ids = get_matching_experiments(v1_encode_ids, v2_dhs_dir)
    v2_dhs_hic_matching_ids = get_matching_experiments(v1_encode_ids, v2_dhs_hic_dir)

    with PdfPages(output_file) as pdf_writer:
        pdf_writer.savefig(
            plot_enh_gene_links(
                v1_dir, v2_dhs_dir, v2_dhs_matching_ids, "E-G Links DHS Only"
            )
        )
        pdf_writer.savefig(
            plot_enh_gene_links(
                v1_dir,
                v2_dhs_hic_dir,
                v2_dhs_hic_matching_ids,
                "E-G Links V1 DHS vs V2 DHS+HIC",
            )
        )


if __name__ == "__main__":
    main()
