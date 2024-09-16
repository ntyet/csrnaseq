#!/bin/bash -l

#SBATCH --job-name=FSR_Sim_zebrafish
#SBATCH --time=2-00:00
#SBATCH -o ./out/slurm.out
#SBATCH -e ./err/slurm.err
#SBATCH --cpus-per-task=20
#SBATCH --array=1-8
enable_lmod
module load container_env R
crun R CMD BATCH --no-save --no-restore 3-Simulation-zebrafish.R ./out/3-Simulation-zebrafish_$SLURM_ARRAY_TASK_ID.Rout


