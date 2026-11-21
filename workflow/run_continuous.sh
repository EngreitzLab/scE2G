#!/bin/bash
#SBATCH -p engreitz,owners,normal
#SBATCH -t 24:00:00
#SBATCH --mem 8G
#SBATCH -o run_continuous.%j.txt
#SBATCH -e run_continuous.%j.txt
#SBATCH --job-name=run_continuous


res_dir=2025_0220_iter10_seed13
this_dir=$OAK/Users/sheth/scE2G_permute_cell_types/scE2G
final_out=$this_dir/results/$res_dir/L1_totals.tsv

# Run the Snakemake pipeline
sbatch $this_dir/workflow/run_snakemake.sh --parsable

while true; do
    if squeue -u shethm | grep -q "run_snak"; then
        echo "run_snakemake is running. Checking again in 10 minutes..."
    else
        echo "run_snakemake is NOT running. Resubmitting..."
        sbatch $this_dir/workflow/run_snakemake.sh --parsable
    fi
    sleep 600  # Wait 10 minutes before checking again
done


# # Wait for Snakemake to finish before deciding whether to resubmit
# wait

# # If Snakemake completes successfully, mark it as done
# if [ -f "$final_out" ]; then
#     echo "Pipeline completed successfully. No resubmission needed."
#     exit 0
# else
#     echo "Pipeline not finished. Resubmitting..."
#     sbatch --dependency=afterany:$SLURM_JOBID "$0"
# fi