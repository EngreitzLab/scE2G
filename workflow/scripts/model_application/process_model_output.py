import click
import pandas as pd
import os


def write_connections_bedpe_format(pred, outfile, score_column):
    # Output a 2d annotation file with EP connections in bedpe format for loading into IGV
    pred = pred.drop_duplicates()

    towrite = pd.DataFrame()

    towrite["chr1"] = pred["chr"]
    towrite["x1"] = pred["start"]
    towrite["x2"] = pred["end"]
    towrite["chr2"] = pred["chr"]
    towrite["y1"] = pred["TargetGeneTSS"]
    towrite["y2"] = pred["TargetGeneTSS"]
    towrite["name"] = pred["TargetGene"] + "_" + pred["name"]
    towrite["score"] = pred[score_column]
    towrite["strand1"] = "."
    towrite["strand2"] = "."

    towrite.to_csv(outfile, header=False, index=False, sep="\t")


@click.command()
@click.option("--predictions_file", required=True)
@click.option("--score_column", required=True)
@click.option("--bedpe_output", required=True)
def main(predictions_file, score_column, bedpe_output):
    pred_thresh = pd.read_csv(predictions_file, sep="\t")
    pred_thresh = pred_thresh.loc[pred_thresh.TargetGeneIsExpressed]
    write_connections_bedpe_format(pred_thresh, bedpe_output, score_column)


if __name__ == "__main__":
    main()
