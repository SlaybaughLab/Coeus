#!/bin/sh

#SBATCH --time=02:30:00
# Job name:
#SBATCH --job-name=mc20
# QoS:
#SBATCH --partition=savio2
# Account:
#SBATCH --ntasks=20
#SBATCH --output=arrayJob_%A_%a.out
#SBATCH --error=arrayJob_%A_%a.err
# Array:
#SBATCH --array=1,2,3,4,5,7,8,9,10

module load openmpi
mpirun mcnp6.mpi i=mcnp_pop$SLURM_ARRAY_TASK_ID.txt o=mcnp_pop$SLURM_ARRAY_TASK_ID.out
