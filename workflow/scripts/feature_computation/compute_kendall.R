## Compute Kendall correlation

# Load required packages
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(genomation)
  library(foreach)
  library(Signac)
  library(Seurat)
  library(Rcpp)
  library(data.table)
  library(Matrix)
  library(anndata)
  library(tools)
  library(dplyr)
})

## Define functions --------------------------------------------------------------------------------

# Calculate the difference between concordant and disconcordant pairs from a sorted logical matrix
cppFunction('
NumericVector count_diff(LogicalMatrix y_matrix_sorted) {
    int n = y_matrix_sorted.nrow();
    int m = y_matrix_sorted.ncol();
    NumericVector result(m);
    for (int j = 0; j < m; j++) {
        long long concordant = 0;
        long long disconcordant = 0;
        long long cumsum = 0;
        for (int i = 0; i < n; i++) {
            bool tmp = y_matrix_sorted(i, j);
            cumsum += tmp;
            if (tmp) {
                disconcordant += (i + 1 - cumsum);
            } else {
                concordant += cumsum;
            }
        }
        result[j] = static_cast<double>(concordant - disconcordant);
    }
    return result;
}
')

# Compute Kendall correlation between a single gene and multiple enhancers
kendall_one_gene = function(x, y.matrix){
  
  # Sort x in decreasing order and accordingly sort y.matrix
  ord = order(x, 
              decreasing = T)
  x.sorted = x[ord]
  y.matrix.sorted = 
    y.matrix[ord, ,drop = F]
  
  # Calculate initial differences between concordant and disconcordant pairs
  n.diff = count_diff(as.matrix(y.matrix.sorted))
  
  # Adjust differences for ties in x
  x.ties = unique(x.sorted[duplicated(x.sorted)])
  for (x.tie in x.ties) {
    n.diff = 
      n.diff - 
      count_diff(as.matrix(y.matrix.sorted[x.sorted == x.tie, ,drop = F]))
  }
  
  # Calculate Kendall's tau-b coefficient
  l = length(x)
  s = colSums(y.matrix)
  tx = table(x)
  
  n0 = choose(l, 2)
  n1 = sum(choose(tx, 2))
  n2 = (s*(s-1) + (l-s)*(l-s-1))/2
  
  tau_b = n.diff / sqrt((n0 - n1) * (n0 - n2))
  
  return(tau_b)
}


# Compute Kendall correlation between a mutliple genes and multiple enhancers
kendall_mutliple_genes = function(bed.E2G,
                                  data.RNA,
                                  data.ATAC,
                                  colname.gene_name = "gene_name",
                                  colname.enhancer_name = "peak_name",
                                  colname.output = "Kendall") {
  
  # Filter E2G pairs based on presence in RNA and ATAC data
  bed.E2G.filter = 
    bed.E2G[mcols(bed.E2G)[,colname.gene_name] %in% rownames(data.RNA) &
              mcols(bed.E2G)[,colname.enhancer_name] %in% rownames(data.ATAC)] 
  

  
  # Compute Kendall correlation for each gene
  bed.E2G.output <- foreach(gene.name = unique(mcols(bed.E2G.filter)[,colname.gene_name]),
                            .combine = 'c') %do% {
                              
                              bed.E2G.tmp <- bed.E2G.filter[mcols(bed.E2G.filter)[,colname.gene_name] == gene.name]
                              
                              mcols(bed.E2G.tmp)[, colname.output] = 
                                kendall_one_gene(as.numeric(data.RNA[gene.name, ]),
                                                 t(data.ATAC[mcols(bed.E2G.tmp)[,colname.enhancer_name], , drop = F]))
                              bed.E2G.tmp
                            }
  return(bed.E2G.output)
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

	# remove genes not in our gene universe	
	row_sub <- intersect(rownames(rna_matrix), names(gene_key)) # gene ref names
	rna_matrix_filt <- rna_matrix[row_sub,] # still gene ref names
	rownames(rna_matrix_filt) <- gene_key[row_sub] # converted to abc names

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
kendall_predictions_path = snakemake@output$kendall_predictions
umi_count_path = snakemake@output$umi_count

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

# Compute Kendall correlation
pairs.E2G = kendall_mutliple_genes(pairs.E2G,
                                   matrix.rna_filt,
                                   matrix.atac,
                                   colname.gene_name = "TargetGene",
                                   colname.enhancer_name = "PeakName",
                                   colname.output = "Kendall")

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
