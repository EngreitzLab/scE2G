#!/bin/bash
#SBATCH -p engreitz,owners,normal
#SBATCH -t 24:00:00  
#SBATCH --mem=64G
#SBATCH -o run_snakemake.%j.txt
#SBATCH -e run_snakemake.%j.txt
#SBATCH --job-name=run_snakemake


# load Conda (if necessary)
source ~/.bashrc  # Ensure Conda is available in non-interactive shells
conda activate run_snakemake

this_dir=$OAK/Users/sheth/scE2G_permute_cell_types/scE2G

# run Snakemake
snakemake --configfile config/config.yaml --use-conda --profile $HOME/.config/snakemake/slurm --unlock

snakemake --configfile config/config.yaml --use-conda --profile $HOME/.config/snakemake/slurm \
    --rerun-incomplete \
    --keep-going \
    --max-jobs-per-second 0.5 \
    --max-status-checks-per-second 5

# Exit 
exit $?
