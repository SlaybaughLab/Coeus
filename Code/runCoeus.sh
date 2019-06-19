#!/bin/bash
# Job name:
#SBATCH --job-name=coeus
#
# Partition:
#SBATCH --partition=savio2
#
# QoS:
# #SBATCH --qos=nuclear_savio_normal
# #SBATCH --qos=savio_normal
#SBATCH --qos=savio_debug
#
# Account:
# #SBATCH --account=co_nuclear
#SBATCH --account=fc_neutronics
#
# Processors:
#SBATCH --nodes=1
#
# Wall clock limit:
#SBATCH --time=00:30:00
#
# SLURM Output File
#SBATCH --output=slurm.out
#
# SLURM Error File
#SBATCH --error=slurm.err

. .bash_profile
## Run command
#python Coeus.py --r=n --log=INFO --qos=nuclear_savio_normal --account=co_nuclear --timeout=00:30:00 --partition=savio --schedule=slurm
python Coeus.py --r=n --log=INFO --qos=savio_debug --account=fc_neutronics --timeout=00:30:00 --partition=savio2 --scheduler=slurm
