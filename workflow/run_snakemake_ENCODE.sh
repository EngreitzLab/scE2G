#!/bin/bash
#SBATCH -p engreitz,owners,normal
#SBATCH -t 24:00:00  
#SBATCH --mem=8G
#SBATCH -o ENCODE.%j.txt
#SBATCH -e ENCODE.%j.txt
#SBATCH --job-name=ENCODE


# load conda (if necessary)
source ~/.bashrc  # ensure Conda is available in non-interactive shells
conda activate run_snakemake

this_dir=$OAK/Users/sheth/scE2G_temp/scE2G

# run Snakemake
snakemake -s $this_dir/workflow/Snakefile --configfile $this_dir/config/mine/config_ENCODE.yaml --use-conda --profile $HOME/.config/snakemake/slurm --unlock

snakemake -s $this_dir/workflow/Snakefile --configfile $this_dir/config/mine/config_ENCODE.yaml --use-conda --profile $HOME/.config/snakemake/slurm \
    --rerun-incomplete \
    --keep-going \
    --max-jobs-per-second 0.5 \
    --max-status-checks-per-second 5 \
    --set-resources compute_kendall:mem_mb=240000 compute_kendall:runtime=1440
# Exit 
exit $?
