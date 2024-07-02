#!/bin/bash -l

#SBATCH --job-name=FSR_Sim_RFI_projection_on_Line
#SBATCH --time=2-00:00
#SBATCH -o ./out/slurm.out
#SBATCH -e ./err/slurm.err
#SBATCH --cpus-per-task=20
#SBATCH --array=1-24
enable_lmod
module load container_env R
crun R CMD BATCH --no-save --no-restore 3-Simulation-RFI_projection_on_Line.R ./out/3-Simulation-RFI_projection_on_Line_$SLURM_ARRAY_TASK_ID.Rout


