#!/bin/bash
#SBATCH -p engreitz,owners,normal
#SBATCH -t 24:00:00  
#SBATCH --mem=8G
#SBATCH -o ds2_train.%j.txt
#SBATCH -e ds2_train.%j.txt
#SBATCH --job-name=ds2_train


# load conda (if necessary)
source ~/.bashrc  # ensure Conda is available in non-interactive shells
conda activate run_snakemake

this_dir=$OAK/Users/sheth/scE2G_temp/scE2G

# run Snakemake
snakemake -s $this_dir/workflow/Snakefile_training --configfile $this_dir/config/mine/config_training_ds2.yaml --use-conda --profile $HOME/.config/snakemake/slurm --unlock

snakemake -s $this_dir/workflow/Snakefile_training --configfile $this_dir/config/mine/config_training_ds2.yaml --use-conda --profile $HOME/.config/snakemake/slurm \
    --rerun-incomplete \
    --keep-going \
    --max-jobs-per-second 0.5 \
    --max-status-checks-per-second 5

# Exit 
exit $?
