import os
import click
import pandas as pd
import numpy as np


def generate_num_sum_enhancers(
    pred_file,
    enh_df,
    distance_threshold,
    chr_sizes,
    enh_midpoint,
    enh_expanded,
    pred_slim,
    enh_pred_int,
    out_num,
    out_sum,
):

    # write enhancer midpoint file
    enh_df["midpoint"] = (enh_df["start"] + enh_df["end"]) / 2
    enh_df["midpoint"] = enh_df["midpoint"].apply(np.int64)
    enh_df[["chr", "midpoint", "midpoint", "name"]].to_csv(
        enh_midpoint,
        sep="\t",
        index=False,
        header=False,
    )

    # expand enhancer to distance bounds
    os.system(
        "bedtools slop -b {} -i {} -g {} > {}".format(
            distance_threshold, enh_midpoint, chr_sizes, enh_expanded
        )
    )

    # select columns from EnhancerPredictionsAllPutative
    os.system(
        "zcat {} | csvtk cut -t -f chr,start,end,name,activity_base | sed '1d' > {}".format(
            pred_file, pred_slim
        )
    )

    # intersect with expanded enhancer regions and count n enhancers
    os.system(
        "bedtools intersect -a {} -b {} -wa -wb | sort -u > {}".format(
            enh_expanded,
            pred_slim,
            enh_pred_int,
        )
    )

    # read in interesct file
    int_df = pd.read_csv(
        enh_pred_int,
        sep="\t",
        header=None,
    )

    #
    int_df_nonself = int_df[
        int_df[3] != int_df[7]
    ]  # filter out target enhancer from the count
    num_enh = (
        int_df_nonself.groupby([3]).size().reset_index(name="count")
    )  # count number of intersecting elements per enhancer
    sum_enh = (
        int_df_nonself.groupby([3])[8].sum().reset_index(name="sum")
    )  # sum activity of intersecting elements per enhancer

    num_enh.to_csv(
        out_num,
        sep="\t",
        header=False,
        index=False,
    )

    sum_enh.to_csv(
        out_sum,
        sep="\t",
        header=False,
        index=False,
    )


@click.command()
@click.option("--abc_predictions")
@click.option("--enhancer_list")
@click.option("--distance_threshold_kb")
@click.option("--chr_sizes")
@click.option("--enh_midpoint")
@click.option("--enh_expanded")
@click.option("--pred_slim")
@click.option("--enh_pred_int")
@click.option("--out_num")
@click.option("--out_sum")
def main(
    abc_predictions,
    enhancer_list,
    distance_threshold_kb,
    chr_sizes,
    enh_midpoint,
    enh_expanded,
    pred_slim,
    enh_pred_int,
    out_num,
    out_sum,
):

    enh_df = pd.read_csv(
        enhancer_list, sep="\t", usecols=["chr", "start", "end", "name"]
    )

    distance_threshold = int(distance_threshold_kb) * 1000  # convert from kb to bp

    generate_num_sum_enhancers(
        abc_predictions,
        enh_df,
        distance_threshold,
        chr_sizes,
        enh_midpoint,
        enh_expanded,
        pred_slim,
        enh_pred_int,
        out_num,
        out_sum,
    )


if __name__ == "__main__":
    main()
