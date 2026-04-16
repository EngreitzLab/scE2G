## Generate single-cell ATAC-seq matrix for a specific cluster

# Load required packages
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(data.table)
  library(Signac)
  library(anndata)
  library(tools)
  library(Seurat)
})

# Import parameters from Snakemake
kendall_pairs_path = snakemake@input[["kendall_pairs_path"]]
atac_frag_path = snakemake@input[["atac_frag_path"]]
rna_matrix_path = snakemake@input[["rna_matrix_path"]]
cell_bc_path = snakemake@input[["cell_barcodes_path"]]
atac_matrix_path = snakemake@output[["atac_matrix_path"]]

max_cell_count = snakemake@params[["max_cell_count"]]

# Read the enhancer-gene pairs and extract unique peaks.
# Only chr/start/end (for FeatureMatrix) and PeakName (for dedup and
# row naming) are used; skip TargetGene and PairName.
pairs.e2g = makeGRangesFromDataFrame(
  fread(kendall_pairs_path, select = c("chr","start","end","PeakName")),
  keep.extra.columns = TRUE,
  starts.in.df.are.0based = TRUE)
bed.peaks = pairs.e2g[!duplicated(mcols(pairs.e2g)[,"PeakName"])]
mcols(bed.peaks) = NULL

# Read the rna matrix to extract cell name
if (file_ext(rna_matrix_path) == "h5ad") {
rna_matrix <- t(read_h5ad(rna_matrix_path)$X)
} else if (file_ext(rna_matrix_path) == "gz") {
rna_matrix = read.csv(rna_matrix_path,
						row.names = 1,
						check.names = F)
} else if (file.info(rna_matrix_path)$isdir) { # assume sparse matrix format
rna_matrix = Read10X(rna_matrix_path, gene.column=1)
} else {
message("Please provide a supported RNA matrix format.")
}

rna.cells = colnames(rna_matrix)
message("Number of cells in RNA matrix: ", length(rna.cells))
message("Example cell: ", rna.cells[1])

# if we need to subset based on atac cells, intersect with fragment file cells 
if (file_test("-f", cell_bc_path)) {
  atac.cells <- readLines(cell_bc_path)
  message("Number of cells from fragment file: ", length(atac.cells))
  message("Example cell: ", atac.cells[1])

  cells.use <- intersect(rna.cells, atac.cells)
} else {
	cells.use <- rna.cells
}

# If the cell count exceeds the max_cell_count, 
# randomly extract max_cell_count cells
if(length(cells.use) > max_cell_count){
  set.seed(123)
  cells.use <- sample(cells.use, max_cell_count)
}

names(cells.use) <- cells.use
message("Number of cells to use: ", length(cells.use))

# Create a list to store Signac Fragment object
list.fragments = list()
list.fragments[[1]] =
  CreateFragmentObject(path = atac_frag_path,
                       cells = cells.use)

# Construct the ATAC-seq matrix for the specified cluster
atac.matrix <- FeatureMatrix(
  fragments = list.fragments,
  features = bed.peaks,
  cells = cells.use
)
# Signac::FeatureMatrix names rows via GRangesToString (1-based start).
# Rewrite rownames using the 0-based BED start so they match the
# PeakName column in Pairs.Kendall.tsv.gz, which compute_kendall.R
# uses to subset the matrix. FeatureMatrix preserves the order of
# the input `features`, so this is a positional rename.
rownames(atac.matrix) <- paste(seqnames(bed.peaks),
                               start(bed.peaks) - 1L,
                               end(bed.peaks),
                               sep = "-")
message("Example cell in ATAC matrix: ", colnames(atac.matrix)[1])


# Save ATAC-seq matrix
saveRDS(atac.matrix,
	atac_matrix_path)
