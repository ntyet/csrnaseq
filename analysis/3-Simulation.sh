#!/bin/bash -l

#SBATCH --job-name=FSR_Sim
#SBATCH --time=2-00:00
#SBATCH -o ./out/slurm_%A-%a.out
#SBATCH -e ./err/slurm_%A-%a.err
#SBATCH --cpus-per-task=20
#SBATCH --array=1-24
enable_lmod
module load container_env R
crun R CMD BATCH --no-save --no-restore 3-Simulation.R ./out/3-Simulation_$SLURM_ARRAY_TASK_ID.Rout


