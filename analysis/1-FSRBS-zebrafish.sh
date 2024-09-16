#!/bin/bash -l

#SBATCH --job-name=zebrafish
#SBATCH --time=2-00:00
#SBATCH -o ./out/slurm.out
#SBATCH -e ./err/slurm.err
#SBATCH --cpus-per-task=20
#SBATCH --array=1
enable_lmod
module load container_env R
crun R CMD BATCH --no-save --no-restore 1-FSRBS_zebrafish.R ./out/1-FSRBS_zebrafish_$SLURM_ARRAY_TASK_ID.Rout


