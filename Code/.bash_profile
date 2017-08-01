# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi

# Modules
module load python/2.7.8
module load numpy
module load scipy
module load matplotlib
module load setuptools
module load tables
module load nose
module load mpi4py

# Aliases

alias groupdir="cd /global/home/groups/co_nuclear"

alias stat='squeue | grep youdongz'

alias nuc='squeue -q nuclear_normal'

alias scratch="cd /global/scratch/youdongz/"

# wwall -j jobid


# User specific environment and startup programs
PATH=$PATH:$HOME/bin:/global/home/groups/co_nuclear/bin:/global/home/groups/co_nuclear/ADVANTG/bin
DATAPATH=/global/scratch/co_nuclear/MCNP/MCNP_DATA/

export PATH
export DATAPATH

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/global/home/groups/co_nuclear/pyne/lib
export PYTHONPATH=$PYTHONPATH:/global/home/groups/co_nuclear/python-pkgs

export PYTHONPATH=$PYTHONPATH:/global/home/groups/co_nuclear/pyne/lib/python2.7/site-packages
