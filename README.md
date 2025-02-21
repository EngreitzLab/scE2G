# scE2G
This project is licensed under the terms of the MIT license.

The preprint describing scE2G is available [here](https://www.biorxiv.org/content/10.1101/2024.11.23.624931v1).

Input: Single-cell ATAC-seq or paired ATAC and RNA-seq (mulitome) data per cell cluster

Output: Genomewide enhancer-gene regulatory link predictions per cell cluster

The pipeline consists of the following components:
1. Compute ABC model predictions for each cell cluster
2. Generate E2G features from ABC predictions
3. If running scE2G (Multiome): compute Kendall correlation and/or ARC-E2G score for each cell cluster
4. Combine 2 & 3 to construct a feature file to be used as input to predictive model
5. (optional) Train predictive model using CRISPR-validated E-G pairs from K562 dataset
6. Apply trained model to make predictions by assigning a score to each E-G pair

## Set up
Clone the repo and set it up for submodule usage
```
git clone --recurse-submodules https://github.com/EngreitzLab/scE2G.git
git config --global submodule.recurse true
```

When running for the first time, the conda environments have to be setup. 
We highly recommend using the environment specified in `workflow/envs/run_snakemake.yml`, which specifies the exact package versions compatible with the pipeline.
For speed, we include mamba in this recommended environment.

```
conda config --set channel_priority flexible  # Make sure your conda config uses flexible channel packaging to prevent unsatisfiable errors
conda env create -f workflow/envs/run_snakemake.yml
conda activate run_snakemake

```

## Apply model
Before running this workflow, users should perform clustering to define cell clusters through regular single-cell analysis (such as Seurat & Signac).
Required input data includes (refer to the example data in the `resources/example_chr22_multiome_cluster` folder):
1. Pseudobulk fragment files and their corresponding *.tbi index files in the same directory for each cell cluster
	- Must be sorted by coordinates and gzipped. If it is not sorted, you can use the sortBed tool from bedtools: `sortBed atac_fragments.unsorted.tsv > atac_fragments.tsv`.
	- Must have 5 columns (no header) corresponding to chr, start, end, cell_name, read_count (usually just 1)
	- The cell_name column should correspond to the cell names in the RNA count matrix. All cells in the RNA matrix must be represented in the fragment file and vice versa.
	- To create a .tbi index, use `bgzip` from [HTSlib](https://github.com/samtools/htslib) instead of gzip to compress the fragment file: `bgzip atac_fragments.tsv`, then generate the corresponding .tbi index file using `tabix -p bed atac_fragments.tsv.gz`.
2. For scE2G (Multiome): RNA count matrix (gene x cell) for each cell cluster
	- Use unnormalized (raw) counts
	- Ensure there are not duplicated gene names
	- Accepted formats for RNA matrix:
		1. .csv.gz format
		2. .h5ad or .h5 format
		3. Directory with sparse matrix files from 10X (matrix.mtx, genes.tsv or features.tsv, and barcodes.tsv). Note that files are read with the argument `gene.column = 1`. 
	- By default, gene names in the RNA matrix are mapped to our gene reference file via Ensembl ID using gene annotations from the GENCODE v43 release (included in `resources/gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz`). You can modify the annotations to another GENCODE version in in the `gene_annotation` field of `config/config.yaml`. For example, CellRanger uses the GENCODE v32 gene annotations by default; this annotation file can be downloaded [here](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz).

To configure the pipeline:
- Modify `config/config.yaml` to specify paths for results_dir.
- Modify `config/config_cell_clusters.tsv` to specify the cluster name, ATAC fragment file path, and RNA matrix path. Leave the remaining columns empty.
- If running scE2G (ATAC), also leave the RNA matrix path empty, but still include the column.
- Specify model directory from those in `models/` to be applied


Running the pipeline:
```
snakemake -j1 --use-conda
```
This command make take a while the first time you run it, as it needs to build the conda environments. 
But if it takes more than 1 hour, that's usually a bad sign, meaning that you're not using mamba and/or need more memory to build the environment.

Output will show up in the `results/` directory by default, with the structure `results/cell_cluster/model_name/encode_e2g_predictions.tsv.gz`. The score column to use is `E2G.Score.qnorm`. 

## Train model

**Important: Only train models for biosamples matching the corresponding CRISPR data (in this case, K562)**

Modify `config/config_training.yaml` with your model and cell_cluster configs
- `model_config` has columns: model, dataset, ABC_directory, feature_table, polynomial (do you want to use polynomial features?) 
Note that trained models generated using polynomial features cannot directly be used in the **Apply model** workflow
- `cell_cluster_config` has rows representing each "dataset"  in `model_config`, where each "dataset" must correspond to a "cluster" in `cell_cluster_config`
    - If an ABC_directory is not specified for a dataset, its entry in `cell_cluster_config` must also contain the required ABC biosample parameters
    - TO DO: specify how to generate and formats for Kendall parameters
- To apply a trained model, it must contain the following files: 1) `model.pkl`, 2) `feature_table.tsv`, 3) `score_threshold_.XX` , 5) `tpm_threshold_YY` (YY=0 if ATAC-only model), 4) `qnorm_reference.tsv.gz` (single column with header `E2G.Score` that contains raw scores for genomewide predictions)

Running the pipeline:
```
snakemake -s workflow/Snakefile_training -j1 --use-conda
```
Output will show up in the `results_training/` directory by default
