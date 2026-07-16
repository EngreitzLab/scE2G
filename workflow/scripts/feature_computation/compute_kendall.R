## Compute Kendall correlation

# Load required packages
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(genomation)
  library(Signac)
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(anndata)
  library(tools)
  library(dplyr)
  library(tibble)
})

options(scipen = 999)

## Define functions --------------------------------------------------------------------------------

# Writes the (filtered) RNA/ATAC data to disk for the Python correlation step
# and returns the paths, without holding on to anything the caller couldn't
# already free. Split out from run_kendall_python() specifically so the
# caller can rm() its own references to data.RNA/data.ATAC (potentially very
# large, atlas-scale matrices) in between the two calls -- R passes arguments
# by reference until modified, so rm()ing them *inside* a function would only
# drop that function's local binding and free nothing, since the caller's
# variables would still keep the underlying object alive. Without freeing
# them here, R and the Python subprocess would both hold a full copy of the
# same data at once during system2()'s call in run_kendall_python().
prepare_kendall_inputs = function(bed.E2G,
                                  data.RNA,
                                  data.ATAC,
                                  colname.gene_name = "gene_name",
                                  colname.enhancer_name = "peak_name") {

  # Filter E2G pairs based on presence in RNA and ATAC data
  bed.E2G.filter =
    bed.E2G[mcols(bed.E2G)[,colname.gene_name] %in% rownames(data.RNA) &
              mcols(bed.E2G)[,colname.enhancer_name] %in% rownames(data.ATAC)]

  tmp_dir = tempfile("kendall_")
  dir.create(tmp_dir)

  rna_mtx_path = file.path(tmp_dir, "rna_matrix.mtx")
  atac_mtx_path = file.path(tmp_dir, "atac_matrix.mtx")
  rna_genes_path = file.path(tmp_dir, "rna_genes.txt")
  atac_peaks_path = file.path(tmp_dir, "atac_peaks.txt")
  pairs_path = file.path(tmp_dir, "pairs.tsv")

  writeMM(data.RNA, rna_mtx_path)
  writeMM(data.ATAC, atac_mtx_path)
  writeLines(rownames(data.RNA), rna_genes_path)
  writeLines(rownames(data.ATAC), atac_peaks_path)
  fwrite(
    data.frame(
      TargetGene = mcols(bed.E2G.filter)[, colname.gene_name],
      PeakName = mcols(bed.E2G.filter)[, colname.enhancer_name]
    ),
    file = pairs_path, sep = "\t", col.names = FALSE, quote = FALSE
  )

  list(
    bed.E2G.filter = bed.E2G.filter,
    tmp_dir = tmp_dir,
    rna_mtx_path = rna_mtx_path,
    atac_mtx_path = atac_mtx_path,
    rna_genes_path = rna_genes_path,
    atac_peaks_path = atac_peaks_path,
    pairs_path = pairs_path
  )
}

# Invokes the Python correlation script (backed by the fast_kendall_sc Rust
# package) against inputs already written by prepare_kendall_inputs(), and
# reads the results back onto the filtered pairs. This is the ~7-hour
# bottleneck in the old pure-R implementation, now delegated to Python.
run_kendall_python = function(prepared, cores = 1, colname.output = "Kendall", python_script_path) {
  on.exit(unlink(prepared$tmp_dir, recursive = TRUE), add = TRUE)

  results_path = file.path(prepared$tmp_dir, "kendall_results.txt")

  exit_code = system2(
    "python",
    args = c(
      shQuote(python_script_path),
      "--rna-mtx", shQuote(prepared$rna_mtx_path),
      "--rna-genes", shQuote(prepared$rna_genes_path),
      "--atac-mtx", shQuote(prepared$atac_mtx_path),
      "--atac-peaks", shQuote(prepared$atac_peaks_path),
      "--pairs", shQuote(prepared$pairs_path),
      "--out", shQuote(results_path),
      "--threads", cores
    )
  )
  if (exit_code != 0) {
    stop("compute_kendall.py failed with exit code ", exit_code)
  }

  bed.E2G.filter = prepared$bed.E2G.filter
  kendall_values = as.numeric(readLines(results_path))
  if (length(kendall_values) != length(bed.E2G.filter)) {
    stop(
      "Expected ", length(bed.E2G.filter), " Kendall values from compute_kendall.py, got ",
      length(kendall_values)
    )
  }
  mcols(bed.E2G.filter)[, colname.output] = kendall_values

  return(bed.E2G.filter)
}

# helper function to parse GTF
extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- unlist(strsplit(gtf_attributes, " "))
  if(att_of_interest %in% att){
    return(gsub("\"|;","", att[which(att %in% att_of_interest)+1]))
  } else {
    return(NA)}
}

# map gene names from RNA matrix and expression metrics to gene reference used by scE2G via Ensembl ID
map_gene_names <- function(rna_matrix, df_exp, gene_gtf_path, abc_genes_path){
	gene_ref <- fread(gene_gtf_path, header = FALSE, sep = "\t") %>%
		setNames(c("chr","source","type","start","end","score","strand","phase","attributes")) %>%
		dplyr::filter(type == "gene")
	gene_ref$gene_ref_name <- unlist(lapply(gene_ref$attributes, extract_attributes, "gene_name"))
	gene_ref$Ensembl_ID <- unlist(lapply(gene_ref$attributes, extract_attributes, "gene_id"))
	gene_ref <- dplyr::select(gene_ref, gene_ref_name, Ensembl_ID) %>%
		mutate(Ensembl_ID = sub("\\.\\d+$", "", Ensembl_ID)) %>% # remove decimal digits 
		distinct()
	
	abc_genes <- fread(abc_genes_path, col.names = c("chr", "start", "end", "name", "score", "strand", "Ensembl_ID", "gene_type")) %>%
		dplyr::select(name, Ensembl_ID) %>%
		rename(abc_name = name) %>%
		left_join(gene_ref, by = "Ensembl_ID") %>%
    filter(!is.na(gene_ref_name)) %>%
		group_by(Ensembl_ID) %>% # remove cases where multiple genes map to one ensembl ID
		filter(n() == 1) %>%
		ungroup()

	gene_key <- abc_genes$abc_name
	names(gene_key) <- abc_genes$gene_ref_name
  message("Number of genes in 1:1 key: ", length(gene_key))

	# remove genes not in our gene universe	
	row_sub <- intersect(rownames(rna_matrix), names(gene_key)) # gene ref names
	rna_matrix_filt <- rna_matrix[row_sub,] # still gene ref names
	rownames(rna_matrix_filt) <- gene_key[row_sub] # converted to abc names
  message("Number of genes in RNA matrix before matching: ", length(rownames(rna_matrix)))
  message("Number of genes in RNA matrix after matching: ", length(rownames(rna_matrix_filt)))


	# do the same for expression df
	df_exp_filt <- df_exp[row_sub,]
	rownames(df_exp_filt) <- gene_key[row_sub]

	return(list(rna_matrix_filt, df_exp_filt))
}
## -------------------------------------------------------------------------------------------------

# Import parameters from Snakemake
kendall_pairs_path = snakemake@input$kendall_pairs_path
atac_matrix_path = snakemake@input$atac_matrix
rna_matrix_path = snakemake@input$rna_matrix
gene_gtf_path = snakemake@params$gene_gtf
abc_genes_path = snakemake@params$abc_genes
python_script_path = snakemake@params$python_script
kendall_predictions_path = snakemake@output$kendall_predictions
umi_count_path = snakemake@output$umi_count
cell_count_path = snakemake@output$cell_count
gex_out_path = snakemake@output$all_gex

cores = as.integer(snakemake@threads)

# Load candidate E-G pairs
pairs.E2G = readGeneric(kendall_pairs_path,
                        keep.all.metadata = T,
                        header = T)

# Load scATAC matrix
matrix.atac_count = readRDS(atac_matrix_path)
matrix.atac = BinarizeCounts(matrix.atac_count)
rm(matrix.atac_count)

# Load scRNA matrix
if (file_ext(rna_matrix_path) %in% c("h5ad", "h5")) {
  matrix.rna_count <- t(read_h5ad(rna_matrix_path)$X)
} else if (file_ext(rna_matrix_path) == "gz") {
  matrix.rna_count = read.csv(rna_matrix_path,
                              row.names = 1,
                              check.names = F)
  matrix.rna_count = Matrix(as.matrix(matrix.rna_count), sparse = TRUE)
} else if (file.info(rna_matrix_path)$isdir) { # assume sparse matrix format
	matrix.rna_count = Read10X(rna_matrix_path, gene.column=1)
} else {
	message("Please provide a supported RNA matrix format.")
}

matrix.rna_count = matrix.rna_count[,colnames(matrix.atac)]

# save number of UMIs (pre-filtering) and cells
n_umi = sum(matrix.rna_count)
write(n_umi, file = umi_count_path)

n_cells = ncol(matrix.rna_count)
write(n_cells, file = cell_count_path)

# Normalize scRNA matrix 
matrix.rna = NormalizeData(matrix.rna_count)

# Compute gene expression measurements
df.exp_inf = data.frame(mean_log_normalized_rna = rowMeans(matrix.rna),
                        RnaDetectedPercent = rowSums(matrix.rna_count > 0) / ncol(matrix.rna_count),
                        RnaPseudobulkTPM =  rowSums(matrix.rna_count) / sum(matrix.rna_count)*10^6,
                        row.names = rownames(matrix.rna_count))

# subset (normalized) RNA matrix and map names to ABC gene reference; also subset the gene expression measurements
gene_filtered_out = map_gene_names(matrix.rna, df.exp_inf, gene_gtf_path, abc_genes_path)
matrix.rna_filt <- gene_filtered_out[[1]]
df.exp_filt <-  gene_filtered_out[[2]]

df.exp_filt.to_save <- df.exp_filt %>% 
  rownames_to_column(var = "TargetGene") %>% 
  select(TargetGene,
    RNA_meanLogNorm = mean_log_normalized_rna,
    RNA_pseudobulkTPM = RnaPseudobulkTPM,
    RNA_percentCellsDetected = RnaDetectedPercent)

fwrite(df.exp_filt.to_save, 
       file = gex_out_path,
       row.names = F,
       quote = F,
       sep = "\t")

# Compute Kendall correlation: write inputs, then free R's (potentially very
# large, atlas-scale) copies before starting the Python subprocess, so the R
# and Python processes don't both need the full matrices in memory at once.
kendall_inputs = prepare_kendall_inputs(pairs.E2G,
                                        matrix.rna_filt,
                                        matrix.atac,
                                        colname.gene_name = "TargetGene",
                                        colname.enhancer_name = "PeakName")

rm(matrix.rna, matrix.rna_count, matrix.rna_filt, matrix.atac)
gc()

pairs.E2G = run_kendall_python(kendall_inputs,
                               cores = cores,
                               colname.output = "Kendall",
                               python_script_path = python_script_path)

# add gene expression metrics to E2G pairs
mcols(pairs.E2G)[,c("mean_log_normalized_rna",
                    "RnaDetectedPercent",
                    "RnaPseudobulkTPM")] = 
  df.exp_filt[pairs.E2G$TargetGene,]

# Write output to file
df.pairs.E2G = 
  as.data.frame(pairs.E2G)[,c("seqnames",
                              "start",
                              "end",
                              "TargetGene",
                              "PeakName",
                              "PairName",
                              "mean_log_normalized_rna",
                              "RnaDetectedPercent",
                              "RnaPseudobulkTPM",
                              "Kendall")]
colnames(df.pairs.E2G) = 
  c("chr",
    "start",
    "end",
    "TargetGene",
    "PeakName",
    "PairName",
    "mean_log_normalized_rna",
    "RnaDetectedPercent",
    "RnaPseudobulkTPM",
    "Kendall")
fwrite(df.pairs.E2G,
       file = kendall_predictions_path,
       row.names = F,
       quote = F,
       sep = "\t")
