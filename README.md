# scE2G

A computational pipeline for predicting genome-wide enhancer-gene regulatory links from single-cell ATAC-seq or paired ATAC and RNA-seq (multiome) data.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-2024.11.23.624931v1-red.svg)](https://www.biorxiv.org/content/10.1101/2024.11.23.624931v1)
[![CircleCI](https://dl.circleci.com/status-badge/img/gh/EngreitzLab/scE2G/tree/main.svg?style=shield)](https://dl.circleci.com/status-badge/redirect/gh/EngreitzLab/scE2G/tree/main)

<hr>

## Overview

**Input:** Single-cell ATAC-seq or paired ATAC and RNA-seq (multiome) data per cell cluster

**Output:** Genome-wide enhancer-gene regulatory link predictions per cell cluster

### Pipeline components

1. **ABC model predictions** - Compute ABC model predictions for each cell cluster
2. **E2G feature generation** - Generate element-gene features from ABC predictions  
3. **Correlation analysis** *(multiome only)* - Compute Kendall correlation and/or ARC-E2G score for each cell cluster
4. **Feature integration** - Combine components 2 & 3 to construct feature file for predictive model
5. **Model training** *(optional)* - Train predictive model using CRISPR-validated E-G pairs from K562 dataset
6. **Prediction** - Apply trained model to assign scores to each element-gene pair

<hr>

## Installation

### Clone repository

```bash
# Clone the repository with submodules
git clone --recurse-submodules https://github.com/EngreitzLab/scE2G.git

# Initialize and update nested submodules
cd scE2G
git submodule update --init --recursive
```

### Set up environment

We highly recommend using the provided conda environment for compatibility:

```bash
# Configure conda for flexible channel packaging
conda config --set channel_priority flexible

# Create and activate environment
conda env create -f workflow/envs/run_snakemake.yml
conda activate run_snakemake
```

<hr>

## Usage

### Prerequisites

Before running scE2G, perform clustering to define cell clusters through standard single-cell analysis (e.g., Seurat & Signac).

### Input data requirements

See example data in `resources/example_chr22_multiome_cluster/` folder.

#### 1. Pseudobulk fragment files

- **Format:** One file per cell cluster with corresponding `*.tbi` index files
- **Requirements:**
  - Sorted by coordinates, compressed using `bgzip`, and indexed with `tabix`
  - 5 columns (no header) corresponding to: `chr`, `start`, `end`, `cell_name`, `read_count`
  - Cell names must match RNA count matrix column names
  - Must contain exact set of cells represented in RNA matrix
 
**Preparation:**
```bash
# Sort if needed
sort -k1,1 -k2,2n atac_fragments.unsorted.tsv > atac_fragments.tsv

# Compress and index
bgzip atac_fragments.tsv
tabix -p bed atac_fragments.tsv.gz
```

#### 2. RNA count matrix *(multiome only)*

- **Format:** Gene × cell matrix for each cell cluster
- **Requirements:**
  - Use unnormalized (raw) counts
  - No duplicated gene names
  - Supported formats:
    - `.csv.gz`
    - `.h5ad` or `.h5` (may require matching `anndata` version)
    - Sparse matrix directory (`matrix.mtx.gz`, `genes.tsv.gz`/`features.tsv.gz`, `barcodes.tsv.gz`)

**Gene mapping:** By default, genes are mapped via Ensembl ID using GENCODE v43 annotations. Modify `gene_annotation` in `config/config.yaml` for different versions (e.g., GENCODE v32 for CellRanger data).

<details>
<summary><h4>Advanced file format options</h4></summary>

**Pre-processed fragment files**

If your fragment files are already properly sorted and filtered to main chromosomes, you can skip preprocessing steps by setting `fragments_preprocessed: True` in your config file. Use this option only if you are certain that your files:
1. Are sorted with `sort -k1,1 -k2,2n` (chromosome then numerical position)
2. Only contain fragments on chromosomes present in the chromosome sizes reference file

This option allows you to skip the sorting and filtering steps of the pipeline, which can be very resource intensive for large fragment files.

**Cell filtering configurations**

The default pipeline settings assume the each cell cluster has a corresponding fragment file and RNA matrix that contain the exact same cells. If you instead have an RNA matrix containing cells from many clusters, you can avoid making cluster-specific matrices by setting `RNA_matrix_filtered: False` in your config file, and using the same RNA matrix for all clusters in your `cell_cluster_config`. The pipeline will then use the *intersection of cells contained in the cluster-specific fragment file and combined RNA matrix* to compute the Kendall correlation.
Please note:
1. You still must provide cluser-specific fragment files
2. The RNA matrix must meet the formatting requirements indicated above
3. The memory requirements to load a very large RNA matrix may exceed the default estimations in the pipeline.

</details>

### Configuration

1. **Main config** - Edit `config/config.yaml`:
   - Set `results_dir` path
   
2. **Cell clusters** - Edit `config/config_cell_clusters.tsv`:
   - Specify cluster name, ATAC fragment file path, RNA matrix path
   - For ATAC-only analysis: leave RNA matrix path empty but include the column
   
3. **Model selection** - Specify path to model directory, or a comma-separated list of multiple models. Current supported models estimate contact using a power law function of genomic distance:
   - For multiome predictions: `models/multiome_powerlaw_v3`
   - For ATAC-only predictions: `models/scATAC_powerlaw_v3`

### Running the pipeline

```bash
snakemake -j1 --use-conda --configfile config/config.yaml
```

> **Note:** First run may take time to build conda environments. If it exceeds 1 hour, ensure you're using mamba and have sufficient memory.

### Output

#### Key outputs
- **Tabular predictions**
  - `{results_dir}/{cell_cluster}/{model_name}/scE2G_predictions.tsv.gz`: All putative enhancer-gene predictions for a cell cluster
  - `{results_dir}/{cell_cluster}/{model_name}/scE2G_predictions_threshold{model_threshold}.tsv.gz`: Thresholded predictions containing enhancer-gene pairs that pass score threshold and other filtering steps
  - Key score column in these files is `E2G.Score.qnorm` (quantile-normalized scE2G score)
- **Genome-browser files** (produced if `make_IGV_tracks: True` in your config file)
  - `{IGV_dir}/{cell_cluster}/ATAC_norm.bw`: bigWig file with read-depth normalized pseudobulk ATAC signal 
  - `{IGV_dir}/{cell_cluster}/scE2G_predictions_threshold{model_threshold}.bedpe`: bedpe file with filtered enhnancer-gene predictions
- **QC report**
  - `{results_dir}/qc_plots/predictions_qc_report.html`: Report summarizing properties of predictions in comparison to reference values

<details>
<summary><h4>Full output structure</h4></summary>
  
```
{results_dir}/                                         # Main results directory
├── {cell_cluster}/                                      # Outputs for each cell cluster
│   ├── ActivityOnly_features.tsv.gz                       # All element-gene pairs with activity-based features
│   ├── ActivityOnly_plus_external_features.tsv.gz         # All element-gene pairs with activity-based and other features
│   ├── ARC/                                               # ARC-E2G results and intermediate files (multiome only)
│   ├── external_features_config.tsv                       # Configuration for external features
│   ├── feature_table.tsv                                  # Feature table reflecting all models designated for the cell cluster
│   ├── genomewide_features.tsv.gz                         # Genome-wide element-gene pairs with features, formatted for model application
│   ├── Kendall/                                           # Kendall correlation results and intermediate results (multiome only)
│   ├── {model_name}/                                      # scE2G model-specific results (e.g., multiome_powerlaw_v3)
│   │   ├── scE2G_predictions.tsv.gz                         # All predictions with scores
│   │   ├── scE2G_predictions_threshold{threshold}.tsv.gz    # Filtered predictions
│   │   ├── scE2G_predictions_threshold{threshold}_stats.tsv # Properties of filtered predictions
│   │   ├── scE2G_element_list.tsv.gz                        # List of candidate elements and associated features
│   │   └── scE2G_gene_list.tsv.gz                           # List of genes in reference file and associated features
│   ├── Neighborhoods/                                    # Results from "Neighborhoods" step of ABC
│   ├── new_features/                                     # Additional computed features
│   ├── Peaks/                                            # Results from "Peaks" step of ABC
│   ├── Predictions/                                      # Results from "Predictions" step of ABC
│   ├── processed_genes_file.bed                          # Processed gene annotations
│   ├── tagAlign/                                         # Results from converting fragment to tagAlign file
│   └── to_generate.txt                                   # Pipeline generation tracking file
└── qc_plots/                                           # Quality control outputs across clusters
  └── predictions_qc_report.html                          # Report of prediction properties compared to reference
  └── [PDFs of prediction property plots]

{IGV_dir (= results_dir if not defined}/             # Genome-browser results directory (only if make_IGV_tracks is True)
├── {cell_cluster}/                                    # Outputs for each cell cluster
│   ├── ATAC.bw                                          # bigWig with unnormalized pseudobulk ATAC signal
│   ├── ATAC_norm.bw                                     # bigWig with read-depth normalized pseudobulk ATAC signal
└── └── {model_name}/                                    # scE2G model results (e.g., multiome_powerlaw_v3)
       └── scE2G_predictions_threshold{threshold}.bedpe   # bedpe file corresponding to filtered predictions
```
  
</details>

<hr>

## Model training

> **⚠️ Important:** Only train models for biosamples matching the corresponding CRISPR data (currently K562).

### Configuration

Edit `config/config_training.yaml`:

- **`model_config`** columns:
  - `model`: Model name
  - `dataset`: Dataset identifier  
  - `ABC_directory`: ABC results directory
  - `feature_table`: Feature table path
  - `polynomial`: Use polynomial features? (Note: models with polynomial features cannot be directly used in Apply model workflow)

- **`cell_cluster_config`** rows:
  - Each "dataset" in `model_config` must correspond to a "cluster" here
  - If `ABC_directory` not specified, must contain required ABC biosample parameters

### Assembling model directory to use in application workflow

Each model directory must contain:
1. `model.pkl`
2. `feature_table.tsv` 
3. `score_threshold_{score_threshold}`, where `score_threshold` is a value from 0–1 (e.g., `0.177`)
4. `tpm_threshold_{tpm_threshold}`, where `tpm_threshold` is any non-negative value  (use 0 for ATAC-only models)
5. `qnorm_reference.tsv.gz` (single column with header `E2G.Score` containing raw scores)

### Running training workflow

```bash
snakemake -s workflow/Snakefile_training -j1 --use-conda
```

Output appears at the path to `results_dir` specified in `config_training.yaml`.

<hr>

### License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### Citation

If you use scE2G in your research, please cite our preprint:
[bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2024.11.23.624931v1)

### Support

For questions and issues, please use the [GitHub Issues](https://github.com/EngreitzLab/scE2G/issues) page.
