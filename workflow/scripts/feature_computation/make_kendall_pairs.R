## Make enhancer-gene pairs for computing Kendall correlation 
## using non-extended peaks

# Load required packages
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(data.table)
})

# Import parameters from Snakemake
narrowPeak_path = snakemake@input[["narrowPeak"]]
allPutative_path = snakemake@input[["allPutative"]]
kendallPairs_path = snakemake@output[["kendallPairs"]]


# Read the narrowPeak (only chr/start/end are used downstream)
bed.narrowPeak = makeGRangesFromDataFrame(
  fread(narrowPeak_path, select = 1:3, col.names = c("chr","start","end")),
  starts.in.df.are.0based = TRUE)

# Read the Enhancer-Gene pairs in ABC prediction.
# Only chr/start/end (for findOverlaps) and TargetGene (assigned to
# output pairs) are used; skip the ~20 other ABC columns.
bed.allPutative = makeGRangesFromDataFrame(
  fread(allPutative_path, select = c("chr","start","end","TargetGene")),
  keep.extra.columns = TRUE,
  starts.in.df.are.0based = TRUE)

# Find overlaps between narrowPeak and enhancers in ABC prediction
overlaps.res = findOverlaps(bed.narrowPeak,
                            bed.allPutative)

# Extract overlapping peaks and assign target genes to them
pairs.e2g = bed.narrowPeak[overlaps.res@from]
mcols(pairs.e2g)[,"TargetGene"] = mcols(bed.allPutative)[overlaps.res@to,"TargetGene"]

# Generate peak name. GRanges starts are 1-based inclusive; subtract 1
# so PeakName uses the original 0-based BED start (matches ATAC matrix
# rownames written by generate_atac_matrix.R).
mcols(pairs.e2g)[,"PeakName"] =
  paste(seqnames(pairs.e2g),
        start(pairs.e2g) - 1L,
        end(pairs.e2g),
        sep = "-")

# Generate pair name for each pair
mcols(pairs.e2g)[,"PairName"] = 
  paste(mcols(pairs.e2g)[,"PeakName"],
        mcols(pairs.e2g)[,"TargetGene"],
        sep = "_")

# Sort pairs by PairName and remove duplicates
pairs.e2g = pairs.e2g[order(mcols(pairs.e2g)[,"PairName"])]
pairs.e2g = pairs.e2g[!duplicated(mcols(pairs.e2g)[,"PairName"])]

# Convert the Enhancer-Gene pairs to a data frame
df.pairs.e2g = 
  as.data.frame(pairs.e2g)[,c("seqnames",
                              "start",
                              "end",
                              "TargetGene",
                              "PeakName",
                              "PairName")]

# Rename the columns for clarity
colnames(df.pairs.e2g) =
  c("chr",
    "start",
    "end",
    "TargetGene",
    "PeakName",
    "PairName")

# GRanges starts are 1-based inclusive; subtract 1 to restore the
# 0-based BED convention of the input narrowPeak/AllPutative files.
df.pairs.e2g$start = df.pairs.e2g$start - 1L

# Write the enhancer-gene pairs to a file
fwrite(df.pairs.e2g,
       file = kendallPairs_path,
       row.names = F,
       col.names = T,
       quote = F,
       sep = "\t")

