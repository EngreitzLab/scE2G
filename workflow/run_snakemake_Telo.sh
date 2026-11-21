#!/bin/bash
#SBATCH -p engreitz,owners,normal
#SBATCH -t 24:00:00  
#SBATCH --mem=8G
#SBATCH -o telo.%j.txt
#SBATCH -e telo.%j.txt
#SBATCH --job-name=telo


# load conda (if necessary)
source ~/.bashrc  # ensure Conda is available in non-interactive shells
conda activate run_snakemake

this_dir=$OAK/Users/sheth/scE2G_temp/scE2G

# run Snakemake
Snakefile=$this_dir/workflow/Snakefile
config=$this_dir/config/mine/config.yaml

snakemake -s $Snakefile --configfile $config --use-conda --profile $HOME/.config/snakemake/slurm --unlock

snakemake -s $Snakefile --configfile $config --use-conda --profile $HOME/.config/snakemake/slurm \
    --rerun-incomplete \
    --keep-going \
    --max-jobs-per-second 0.5 \
    --max-status-checks-per-second 5 \
    --rerun-triggers software-env

# Exit 
exit $?