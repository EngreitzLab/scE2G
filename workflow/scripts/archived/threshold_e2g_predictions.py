import click
import pandas as pd


def threshold_predictions(all_putative, threshold, include_self_promoter):
    # filter by threshold and promoter category
    filtered_predictions = all_putative[
        all_putative["ENCODE-rE2G.Score.qnorm"] >= threshold
    ]

    # process promoters & self-promoters
    if include_self_promoter:
        return filtered_predictions[
            (filtered_predictions["class"] != "promoter")
            | filtered_predictions["isSelfPromoter"]
        ]
    else:
        return filtered_predictions[filtered_predictions["class"] != "promoter"]


@click.command()
@click.option("--all_predictions_file", required=True)
@click.option("--threshold", type=float, required=True)
@click.option("--include_self_promoter", type=bool, default=True)
@click.option("--output_file", required=True)
def main(all_predictions_file, threshold, include_self_promoter, output_file):
    all_predictions = pd.read_csv(all_predictions_file, sep="\t")
    filtered_predictions = threshold_predictions(
        all_predictions, threshold, include_self_promoter
    )
    filtered_predictions.to_csv(output_file, compression="gzip", sep="\t", index=False)


if __name__ == "__main__":
    main()
