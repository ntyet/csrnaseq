#!/bin/bash -l

#SBATCH --job-name=FSRBS_ModelSize_Simulation
#SBATCH --time=0-01:00
#SBATCH -o ./out/slurm_%A-%a.out
#SBATCH -e ./err/slurm_%A-%a.err
#SBATCH --cpus-per-task=1
enable_lmod
module load container_env R
crun R CMD BATCH --no-save --no-restore 2-FSRBS_ModelSize.R ./out/2-FSRBS_ModelSize.Rout


