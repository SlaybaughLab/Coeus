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




To get the latest repo on Savio:
You will need to generate an SSH key.  Follow the steps here:

https://help.github.com/articles/connecting-to-github-with-ssh/

Then clone the repo in your desired folder location.  



## Licensing Information
\license <a href='../../licensing/LICENSE'>GNU GPLv3.0+ </a>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
