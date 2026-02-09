import os
from pathlib import Path
from typing import List

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

VALUE_KEY = "Value"
SEQ_DEPTH_METRIC = "num_sequencing_reads"
METRICS = [
    "num_enh",
    "num_genes",
    "num_enh_gene_links",
    "num_genes_with_1_enh_min",
    "mean_num_genes_per_enh",
    "mean_num_enh_per_gene",
    "mean_num_enh_per_gene_no_prom",
    "mean_log10_dist_to_tss",
    "mean_enh_region_size",
]


def load_and_save_stat_files(stat_files, output_table) -> List[pd.DataFrame]:
    results = {}
    results_no_ind = {}
    for stat_file in stat_files:
        cell_cluster = Path(stat_file).parts[-3]
        model_name = Path(stat_file).parts[-2]
        pred_id = f"{cell_cluster}_{model_name}"
        results[pred_id] = pd.read_csv(stat_file, sep="\t").set_index("Metric").assign(
            cell_cluster=cell_cluster,
            model_name=model_name
        )
        results_no_ind[pred_id] = pd.read_csv(stat_file, sep="\t").assign(
            cell_cluster=cell_cluster,
            model_name=model_name
        )

    # save merged dataframe
    pd.concat(results_no_ind, ignore_index=True).to_csv(
        output_table, sep="\t", index=False
    )

    return results


def snake_to_title(snake_case):
    words = snake_case.split("_")
    # Capitalize the first letter of each word except the first one
    title = " ".join(word.capitalize() for word in words)
    return title


def plot_cdf(stats, metric, biosample_type=None):
    plt.clf()
    points = [df.loc[metric, VALUE_KEY] for df in stats.values()]
    mean, median = np.mean(points) / 1e6, np.median(points) / 1e6
    ax = sns.kdeplot(points, fill=True, bw_adjust=0.5, color="#006eae", label="Density")
    title = snake_to_title(metric)
    if biosample_type:
        title = f"{biosample_type}: {title}"

    ax.set_title(title)
    plt.xlabel("Data")
    plt.ylabel("Density")
    plt.scatter(
        [], [], label=f"n={len(points)}\nMean={mean:.2f}M\nMedian={median:.2f}M"
    )
    plt.legend()
    return ax.get_figure()


def plot_violin(stats, metric, biosample_type=None):
    plt.clf()
    points = [df.loc[metric, VALUE_KEY] for df in stats.values()]
    mean, median = np.mean(points), np.median(points)
    ax = sns.violinplot(y=points)
    title = snake_to_title(metric)
    if biosample_type:
        title = f"{biosample_type}: {title}"

    ax.set_title(title)
    plt.scatter([], [], label=f"n={len(points)}\nMean={mean:.2f}\nMedian={median:.2f}")
    plt.legend()
    return ax.get_figure()


def plot_swarm_box(stats, metric, biosample_type=None):
    plt.clf()
    points = [df.loc[metric, VALUE_KEY] for df in stats.values()]
    print(points)
    mean, median = np.mean(points), np.median(points)
    ax = sns.swarmplot(y=points)
    ax = sns.boxplot(y=points)
    title = snake_to_title(metric)
    if biosample_type:
        title = f"{biosample_type}: {title}"

    ax.set_title(title)
    plt.scatter([], [], label=f"n={len(points)}\nMean={mean:.2f}\nMedian={median:.2f}")
    plt.legend()
    return ax.get_figure()


def plot_metric_vs_seq_depth(stats, metric, biosample_type=None):
    plt.clf()
    points = [df.loc[metric, VALUE_KEY] for df in stats.values()]
    seq_depths = [df.loc[SEQ_DEPTH_METRIC, VALUE_KEY] / 1e6 for df in stats.values()]
    ax = sns.scatterplot(x=seq_depths, y=points, label=f"n={len(points)}")
    title = f"{snake_to_title(metric)} vs Sequencing Depth"
    if biosample_type:
        title = f"{biosample_type}: {title}"
    ax.set_title(title)
    plt.xlabel("Sequencing depth (millions of reads or fragments)")
    plt.legend()
    return ax.get_figure()


def plot_scatter(stat_dfs, metadata_df, metric):
    plt.clf()
    X, Y = [], []
    K562_X, K562_Y = [], []
    for cell_cluster in stat_dfs:
        row = metadata_df[metadata_df["CellClusterID"] == cell_cluster].iloc[0]
        num_fragments = row["nCells"] * row["MeanATACFragmentsPerCell"]
        metric_val = stat_dfs[cell_cluster].loc[metric, VALUE_KEY]
        if "K562" in row["ManualAnnotationLabel"]:
            K562_X.append(int(num_fragments))
            K562_Y.append(metric_val)
        else:
            X.append(int(num_fragments))
            Y.append(metric_val)

    ax = sns.scatterplot(x=X, y=Y, label=f"n={len(X)}")
    ax = sns.scatterplot(x=K562_X, y=K562_Y, color="red", label=f"K562 n={len(K562_X)}")
    plt.xscale("log")
    plt.axvline(x=2e6, color="#c5373d", linestyle="-", label=f"2 million fragments")
    ax.set_xlabel("Number unique fragments")
    ax.set_ylabel(snake_to_title(metric))
    ax.set_title(snake_to_title(metric))
    plt.legend()

    return ax.get_figure()


def add_intro_page(pdf_writer):
    plt.clf()
    plt.axis("off")
    plt.text(
        0.5,
        1,
        "QC Plots for ENCODE_rE2G\n",
        ha="center",
        va="center",
        fontsize="large",
        wrap=True,
    )

    plt.text(
        0.5,
        0.5,
        "N represents number of datasets\n\n"
        "The metrics represent E-G pairs after applying a threshold\n\n"
        "Metric definitions are found here: https://github.com/EngreitzLab/ENCODE_rE2G/blob/main/workflow/scripts/get_stats.py",
        ha="center",
        va="center",
        fontsize="medium",
        wrap=True,
    )
    pdf_writer.savefig()


def save_outlier_stats(stat_dfs, metric, pdf_writer):
    pred_values = {}
    for this_df in stat_dfs.values():
        pred_id = f"{this_df['cell_cluster'].iloc[0]}_{this_df['model_name'].iloc[0]}"
        pred_values[pred_id] =  this_df.loc[metric, VALUE_KEY]

    sorted_items = sorted(pred_values.items(), key=lambda item: item[1])

    top_5 = [f"{item[0]}: {item[1]}" for item in list(reversed(sorted_items))[:5]]
    bottom_5 = [f"{item[0]}: {item[1]}" for item in sorted_items[:5]]

    plt.clf()
    plt.axis("off")
    plt.text(
        0.5,
        1,
        f"{metric}",
        ha="center",
        va="center",
        fontsize="large",
        wrap=True,
    )
    nl = "\n"
    plt.text(
        0.5,
        0.5,
        f"Top 5 datasets:{nl}"
        f"{nl.join(top_5)}"
        f"{nl*3}"
        f"Bottom 5 datasets:{nl}"
        f"{nl.join(bottom_5)}"
        f"{nl*3}",
        ha="center",
        va="center",
        fontsize="medium",
        wrap=True,
    )
    pdf_writer.savefig()


def split_stats(stats, metadata_df):
    cell_stats, tissue_stats = {}, {}
    for cluster, stat in stats.items():
        try:
            experiment = [
                encode_id
                for encode_id in cluster.split("_")
                if encode_id.startswith("ENC")
            ][0]
            matching_row = metadata_df[
                metadata_df["DNase Experiment accession"] == experiment
            ].iloc[0]
        except Exception as e:
            import pdb

            pdb.set_trace()
            print()
        if matching_row["Biosample type"] == "tissue":
            tissue_stats[cluster] = stat
        else:
            cell_stats[cluster] = stat

    return cell_stats, tissue_stats


@click.command()
@click.option("--output_file", type=str, default="qc_plots.pdf")
@click.option("--output_table", type=str, default="all_qc_stats.tsv")
@click.argument("stat_files", nargs=-1, type=click.Path(exists=True))
@click.option("--y2ave_metadata", type=str)
@click.option("--encode_metadata", type=str)
def main(output_file, output_table, stat_files, y2ave_metadata, encode_metadata):
    if y2ave_metadata:
        metadata_df = pd.read_csv(y2ave_metadata, sep="\t")
    stats = load_and_save_stat_files(stat_files, output_table)
    if encode_metadata:
        metadata_df = pd.read_csv(encode_metadata, sep="\t")
        cell_stats, tissue_stats = split_stats(stats, metadata_df)
    with PdfPages(output_file) as pdf_writer:
        add_intro_page(pdf_writer)
        pdf_writer.savefig(plot_cdf(stats, SEQ_DEPTH_METRIC))
        for metric in METRICS:
            if encode_metadata:
                for biosample_type, subtype_stats in [
                    ("All", stats),
                    ("Cells", cell_stats),
                    ("Tissues", tissue_stats),
                ]:
                    pdf_writer.savefig(
                        plot_violin(subtype_stats, metric, biosample_type)
                    )
                    pdf_writer.savefig(
                        plot_swarm_box(subtype_stats, metric, biosample_type)
                    )
                    pdf_writer.savefig(
                        plot_metric_vs_seq_depth(subtype_stats, metric, biosample_type)
                    )
                    save_outlier_stats(subtype_stats, metric, pdf_writer)
            else:
                pdf_writer.savefig(plot_violin(stats, metric))
                pdf_writer.savefig(plot_swarm_box(stats, metric))
                pdf_writer.savefig(plot_metric_vs_seq_depth(stats, metric))
                save_outlier_stats(stats, metric, pdf_writer)
            if y2ave_metadata:
                pdf_writer.savefig(plot_scatter(stats, metadata_df, metric))


if __name__ == "__main__":
    main()
