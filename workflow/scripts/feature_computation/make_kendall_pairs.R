## Make enhancer-gene pairs for computing Kendall correlation 
## using non-extended peaks

# Load required packages
suppressPackageStartupMessages({
  library(genomation)
  library(GenomicRanges)
  library(data.table)
})

# Import parameters from Snakemake
narrowPeak_path = snakemake@input[["narrowPeak"]]
allPutative_path = snakemake@input[["allPutative"]]
kendallPairs_path = snakemake@output[["kendallPairs"]]


# Read the narrowPeak
bed.narrowPeak = readGeneric(narrowPeak_path)

# Read the Enhancer-Gene pairs in ABC prediction
bed.allPutative = readGeneric(allPutative_path,
                              keep.all.metadata = T,
                              header = T)

# Find overlaps between narrowPeak and enhancers in ABC prediction
overlaps.res = findOverlaps(bed.narrowPeak,
                            bed.allPutative)

# Extract overlapping peaks and assign target genes to them
pairs.e2g = bed.narrowPeak[overlaps.res@from]
mcols(pairs.e2g)[,"TargetGene"] = mcols(bed.allPutative)[overlaps.res@to,"TargetGene"]

# Generate peak name
mcols(pairs.e2g)[,"PeakName"] = 
  paste(seqnames(pairs.e2g),
        start(pairs.e2g),
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

# Write the enhancer-gene pairs to a file
fwrite(df.pairs.e2g,
       file = kendallPairs_path,
       row.names = F,
       col.names = T,
       quote = F,
       sep = "\t")

