#!/bin/bash
# Job name:
#SBATCH --job-name=coeus
#
# Partition:
#SBATCH --partition=savio
#
# QoS:
# #SBATCH --qos=nuclear_normal
#SBATCH --qos=savio_normal
# #SBATCH --qos=savio_debug
#
# Account:
# #SBATCH --account=co_nuclear
#SBATCH --account=fc_neutronics
#
# Processors:
#SBATCH --nodes=1
#
# Wall clock limit:
#SBATCH --time=72:00:00
#
# SLURM Output File
#SBATCH --output=slurm.out
#
# SLURM Error File
#SBATCH --error=slurm.err

## Run command
python Coeus.py --r=n
