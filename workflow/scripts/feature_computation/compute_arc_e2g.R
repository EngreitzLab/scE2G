# Load required packages
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(genomation)
  library(data.table)
})

## Define functions --------------------------------------------------------------------------------

# Add the maximum Kendall correlation from bed.E2G.B as a new column in bed.E2G.A
AddMax = function(bed.E2G.A,
                  bed.E2G.B,
                  colname.TargetGene.A = "TargetGene",
                  colname.TargetGene.B = "TargetGene",
                  colname.Kendall.A = "Kendall",
                  colname.Kendall.B = "Kendall"){
  
  # Rename seqnames of df.E2G.A as chr_TargetGene
  df.E2G.A.chr_rename = as.data.frame(bed.E2G.A)[,1:3]
  df.E2G.A.chr_rename[,"seqnames"] = paste(seqnames(bed.E2G.A),
                                           mcols(bed.E2G.A)[,colname.TargetGene.A],
                                           sep = "_")
  bed.E2G.A.chr_rename = GRanges(df.E2G.A.chr_rename)
  rm(df.E2G.A.chr_rename)
  
  # Rename seqnames of df.E2G.B as chr_TargetGene
  df.E2G.B.chr_rename = as.data.frame(bed.E2G.B)[,1:3]
  df.E2G.B.chr_rename[,"seqnames"] = paste(seqnames(bed.E2G.B),
                                           mcols(bed.E2G.B)[,colname.TargetGene.B],
                                           sep = "_")
  bed.E2G.B.chr_rename = GRanges(df.E2G.B.chr_rename)
  rm(df.E2G.B.chr_rename)
  
  # Find overlaps between E-G pairs in bed.E2G.A and bed.E2G.B
  overlaps.res = findOverlaps(bed.E2G.A.chr_rename,
                              bed.E2G.B.chr_rename)
  rm(bed.E2G.A.chr_rename)
  rm(bed.E2G.B.chr_rename)
  df.overlaps = as.data.frame(overlaps.res)
  rm(overlaps.res)
  
  # Retain maximum Kendall correlation values
  df.overlaps[,colname.Kendall.B] = mcols(bed.E2G.B)[df.overlaps$subjectHits,colname.Kendall.B]
  rm(bed.E2G.B)
  df.overlaps.sort = df.overlaps[order(df.overlaps[,colname.Kendall.B],decreasing = T),]
  df.overlaps.sort2 = df.overlaps.sort[order(df.overlaps.sort$queryHits),]
  df.overlaps.sort2.unique = df.overlaps.sort2[!duplicated(df.overlaps.sort2$queryHits),]
  mcols(bed.E2G.A)[df.overlaps.sort2.unique$queryHits,
                   colname.Kendall.A] = 
    df.overlaps.sort2.unique[,colname.Kendall.B]
  
  return(bed.E2G.A)
}

# Integrate ABC and Kendall into ARC-E2G
IntegrateABC = function(bed.E2G,
                         colname.ABC,
                         colname.feature,
                         colname.output){
  
  index.exclude_na = !is.na(mcols(bed.E2G)[,colname.ABC]) &
    mcols(bed.E2G)[,colname.ABC] > 0 &
    !is.na(mcols(bed.E2G)[,colname.feature])
  
  bed.E2G.filter = bed.E2G[index.exclude_na]
  
  sd.lnABC = sd(log(mcols(bed.E2G.filter)[,colname.ABC]))
  sd.feature = sd(mcols(bed.E2G.filter)[,colname.feature])
  
  lm.res = 
    lm(scale(log(mcols(bed.E2G.filter)[,colname.ABC])) ~
         scale(mcols(bed.E2G.filter)[,colname.feature]))
  
  beta = 1/coef(lm.res)[2]
  
  mcols(bed.E2G)[,colname.output] = mcols(bed.E2G)[,colname.ABC]
  mcols(bed.E2G)[index.exclude_na,colname.output] =
    exp(log(mcols(bed.E2G)[index.exclude_na,colname.ABC]) +
          mcols(bed.E2G)[index.exclude_na,colname.feature] * sd.lnABC / sd.feature / beta)
  
  return(bed.E2G)
}

# Import parameters from Snakemake
abc_predictions_path = snakemake@input$abc_predictions
kendall_predictions_path = snakemake@input$kendall_predictions
arc_predictions_path = snakemake@output$arc_predictions

# Load ABC data
pairs.E2G.ABC = readGeneric(abc_predictions_path,
                            keep.all.metadata = T,
                            header = T)

# Load Kendall data
pairs.E2G.Kendall = GRanges(fread(kendall_predictions_path))
pairs.E2G.Kendall = 
  pairs.E2G.Kendall[!is.na(mcols(pairs.E2G.Kendall)[,"Kendall"])]

# Add Kendall correlation as a new column of ABC prediction
pairs.E2G.ABC = AddMax(pairs.E2G.ABC,
                       pairs.E2G.Kendall)

# Compute ARC-E2G
message(snakemake@params$abc_score_col)
pairs.E2G.ABC = IntegrateABC(pairs.E2G.ABC,
                              snakemake@params$abc_score_col,
                              "Kendall",
                              "ARC.E2G.Score")

# Copy mean log normalized gene expression
df.gene_exp = mcols(pairs.E2G.Kendall)[,c("TargetGene",
                                          "mean_log_normalized_rna",
                                          "RnaDetectedPercent",
                                          "RnaPseudobulkTPM")]
df.gene_exp = as.data.frame(df.gene_exp)
df.gene_exp = df.gene_exp[!duplicated(df.gene_exp$TargetGene),]
rownames(df.gene_exp) = df.gene_exp$TargetGene
mcols(pairs.E2G.ABC)[,c("mean_log_normalized_rna",
                        "RnaDetectedPercent",
                        "RnaPseudobulkTPM")] = 
  df.gene_exp[pairs.E2G.ABC$TargetGene,c("mean_log_normalized_rna",
                                         "RnaDetectedPercent",
                                         "RnaPseudobulkTPM")]

# Write output to file
df.output = as.data.frame(pairs.E2G.ABC)
colnames(df.output)[1] = "chr"
fwrite(df.output,
       file = arc_predictions_path,
       row.names = F,
       quote = F,
       sep = "\t")
