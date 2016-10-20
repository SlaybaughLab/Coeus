To set up the environment variables to run on Savio, you need to modify the .bashrc or .bash_profile to include:

# User specific modules
module load intel
module load openmpi
module load python/2.7.8
module load numpy
module load scipy
module load matplotlib
module load setuptools
module load tables
module load nose
module load mpi4py

export DATAPATH=/global/scratch/co_nuclear/MCNP/MCNP_DATA/

export PYTHONIOENCODING=utf-8

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/global/home/groups/co_nuclear/pyne/lib

PATH=$PATH:/global/home/users/sbogetic/bin:/global/home/groups/co_nuclear/bin

PATH=$PATH:/global/home/users/sbogetic/bin:/global/home/groups/co_nuclear/LANL/MCNP5/bin/

PATH=$PATH:$HOME/bin:/global/home/groups/co_nuclear/bin:/global/home/groups/co_nuclear/ADVANTG/bin

PYTHONPATH=$PYTHONPATH:/global/home/groups/co_nuclear/python-pkgs/pp-1.6.5

PYTHONPATH=$PYTHONPATH:/global/home/groups/co_nuclear/python-pkgs/pyDOE-0.3.8-py2.7.egg


export PATH


# Aliases
alias groupdir="cd /global/home/groups/co_nuclear"
alias nuc='squeue -q nuclear_normal'