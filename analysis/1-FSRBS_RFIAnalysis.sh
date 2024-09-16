#!/bin/bash -l

#SBATCH --job-name=FSRBS_RFIAnalysis
#SBATCH --time=2-00:00
#SBATCH -o ./out/slurm_%A-%a.out
#SBATCH -e ./err/slurm_%A-%a.err
#SBATCH --cpus-per-task=20
#SBATCH --array=1
enable_lmod
module load container_env R
crun R CMD BATCH --no-save --no-restore 1-FSRBS_RFIAnalysis.R ./out/1-FSRBS_RFIAnalysis_$SLURM_ARRAY_TASK_ID.Rout


