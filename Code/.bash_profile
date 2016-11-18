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
alias nuc='squeue -q nuclear_normal'

# User specific environment and startup programs
PATH=$PATH:$HOME/bin:/global/home/groups/co_nuclear/bin:/global/home/groups/co_nuclear/ADVANTG/bin
DATAPATH=/global/scratch/co_nuclear/MCNP/MCNP_DATA/

export PATH
export DATAPATH

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/global/home/groups/co_nuclear/pyne/lib

PYTHONPATH=$PYTHONPATH:/global/home/groups/co_nuclear/python-pkgs/pp-1.6.5
PYTHONPATH=$PYTHONPATH:/global/home/groups/co_nuclear/python-pkgs/pyDOE-0.3.8-py2.7.egg
export PYTHONPATH=$PYTHONPATH:/global/home/groups/co_nuclear/pyne/lib/python2.7/site-packages
