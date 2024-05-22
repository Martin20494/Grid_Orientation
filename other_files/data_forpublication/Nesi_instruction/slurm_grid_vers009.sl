#!/bin/bash -e

#SBATCH --account         niwa03440
#SBATCH --job-name        lisflood-test
#SBATCH --nodes           1
#SBATCH --ntasks          16
#SBATCH --mem             8G
#SBATCH --time            10:00:00            #Time is in the format of DD-HH:MM:SS
#SBATCH --array           1-50
#SBATCH --output          slurmout/%A_%a.out

module purge
module load netCDF/4.8.1-iimpi-2022a

export PATH=/nesi/project/niwa03440/martin/software2/LISFLOOD-FP-trunk/build:$PATH

srun -n ${SLURM_NTASKS} lisflood /nesi/nobackup/niwa03440/tmn52/GRID_floodevents_002/vers009/param_${SLURM_ARRAY_TASK_ID}/par_${SLURM_ARRAY_TASK_ID}.par
