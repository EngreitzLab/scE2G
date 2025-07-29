# scE2G

A computational pipeline for predicting genome-wide enhancer-gene regulatory links from single-cell ATAC-seq or paired ATAC and RNA-seq (multiome) data.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-2024.11.23.624931v1-red.svg)](https://www.biorxiv.org/content/10.1101/2024.11.23.624931v1)

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

#### 2. RNA count matrix *(Multiome only)*

- **Format:** Gene × cell matrix for each cell cluster
- **Requirements:**
  - Use unnormalized (raw) counts
  - No duplicated gene names
  - Supported formats:
    - `.csv.gz`
    - `.h5ad` or `.h5` (may require matching `anndata` version)
    - Sparse matrix directory (`matrix.mtx.gz`, `genes.tsv.gz`/`features.tsv.gz`, `barcodes.tsv.gz`)

**Gene mapping:** By default, genes are mapped via Ensembl ID using GENCODE v43 annotations. Modify `gene_annotation` in `config/config.yaml` for different versions (e.g., GENCODE v32 for CellRanger data).

### Configuration

1. **Main config** - Edit `config/config.yaml`:
   - Set `results_dir` path
   
2. **Cell clusters** - Edit `config/config_cell_clusters.tsv`:
   - Specify cluster name, ATAC fragment file path, RNA matrix path
   - For ATAC-only analysis: leave RNA matrix path empty but include the column
   
3. **Model selection** - Specify path to model directory in `models/` directory (e.g., `models/multiome_powerlaw_v3`)

### Running the pipeline

```bash
snakemake -j1 --use-conda --configfile config/config.yaml
```

> **Note:** First run may take time to build conda environments. If it exceeds 1 hour, ensure you're using mamba and have sufficient memory.

### Output

Results appear at the `results_dir` path with structure:
```
{results_dir}/{cell_cluster}/{model_name}/scE2G_predictions.tsv.gz
{results_dir}/{cell_cluster}/{model_name}/scE2G_predictions_threshold{model_threshold}.tsv.gz
```

**Key score column:** `E2G.Score.qnorm`

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

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use scE2G in your research, please cite our preprint:
[bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2024.11.23.624931v1)

## Support

For questions and issues, please use the [GitHub Issues](https://github.com/EngreitzLab/scE2G/issues) page.
