#!/usr/bin/env python

#######################################################################################################
#
# Program : Coeus.py.
#
# Contains : The main program for Coeus. Coeus is a fully functional and encapsulated design tool 
#            for Energy Tuning Assemblies.  Designs are optimized to set design criteria using contraints
#            and spectrum objective functions using the Gnowee metaheuristic algorithm.  ADVANTG and MCNP
#            are used to perform the radiation transport calculations.  
#
#            Coeus is designed to run on a cluster.  Use Coeus_local to run on a desktop.  
#
# Author : James Bevins
#
# Last Modified: 24Oct16
#
#######################################################################################################

import numpy as np
import copy as cp

from ETA_Utilities import ETA_Parameters
from Gnowee_Utilities import Gnowee_Settings, Parent, Calc_Fitness, Pop_Update, Timeline
from Metaheuristics import Mat_Levy_Flights, Cell_Levy_Flights, Elite_Crossover, Partial_Inversion, Two_opt
from Metaheuristics import Crossover, Three_opt, Discard_Cells, Mutate
from ADVANTG_Utilities import ADVANTG_Settings, Print_ADVANTG_Input
from NuclearData import Build_Matlib, Calc_Moderating_Ratio
from MCNP_Utilities import MCNP_Settings, MCNP_Geometry, Print_MCNP_Input, Read_Tally_Output
from Utilities import Run_Transport, Event, Meta_Stats

import time
import shutil
import logging
import os
import sys
sys.path.insert(0,os.path.abspath(os.getcwd())+'/Sampling')

import argparse

#-------------------------------------------------------------------------------------------------------------# 
### Local Function definitions

# Print MCNP input Files
def print_MCNP_input_files(step):
    global new_pop
    
    idents=[]
    run_particles=[]
    for i in range(0,len(new_pop)):
        Print_MCNP_Input(eta_params,new_pop[i].geom,new_pop[i].rset,mat_lib,new_pop[i].ident,adv_print=True)
        idents.append(new_pop[i].ident)
        run_particles.append(new_pop[i].rset.nps)
        for m in range(9,len(new_pop[i].geom.matls)):
            if new_pop[i].geom.matls[m]==eta_params.fissile_mat:
                sys.exit()
    logger.info('Gen {} {} finished at {} sec\n'.format(history.tline[-1].g,step,time.time() - start_time))
    return idents, run_particles

# Run MCNP
def run_MCNP_on_algo(args, algo, update_gen, update_feval):
    global ids, particles, pop
    
    if len(ids)>0:
        Run_Transport(ids,args.qos, args.account, args.paritition, particles,code='mcnp6.mpi')
        logger.info('Finished running MCNP at {} sec\n'.format(time.time() - start_time))
    
        # Calculate Fitness
        Calc_Fitness(ids, new_pop, eta_params.spectrum[:,1], eta_params.min_fiss, eta_params.max_weight)
        (changes,feval)=Pop_Update(pop, new_pop, mcnp_set.nps, eta_params, mat_lib, Run_Transport, rr=False) 
        pop=history.update(pop, update_gen, update_feval)
        stats.update(algo,(changes, update_feval + feval))
#-------------------------------------------------------------------------------------------------------------# 
        
"""
@package docstring

    Entry point for the Coeus program.  

    All inputs are optional.  The program will load run inputs in the following order:

    1) User specified path

    2) Run directory defaults (Note: naming convention for files to be loaded from run directory must follow
    default naming convention shown in parameters)

    3) Preset program defaults.  

    This order will be followed for each of the settings files so that some may be ommitted if desired. 

    Parameters
    ==========
obj_path ('obj'): str
    The name and path for the objective function (spectra) file location. The format is a comma delimited spectrum
    with the first column being the upper bin boundaries and the second column the bin flux/fluence.
    [default = ../Inputs/obj_spectrum.csv]
eta_constraints_path ('eta'): str
    The name and path for the file containing the ETA design constraints. The format is a comma delimited key word input 
    file. All keywords are optional.  Non-specified keywords will default to preset program values. 
    [default = ../Inputs/eta_constraints.csv]
gnowee_settings_path ('gs'): str
    The name and path for the file containing the Gnowee search settings. The format is a comma delimited key word input file.
    All keywords are optional.  Non-specified keywords will default to preset program values. 
    [default = ../Inputs/gnowee_settings.csv]
advantg_settings_path ('adv'): str
    The name and path for the file containing the advantg settings. The format is a comma delimited key word input file.
    All keywords are optional.  Non-specified keywords will default to preset program values. 
    [default = ../Inputs/advantg_settings.csv]
mcnp_settings_path ('mcnp'): str
    The name and path for the file containing the mcnp settings. The format is a comma delimited key word input file.
    All keywords are optional.  Non-specified keywords will default to preset program values. 
    [default = ../Inputs/mcnp_settings.csv]
material_lib_path ('mat'): str
    The name and path for the file containing the materials to be included in the problem. The format is a comma delimited 
    key word input file.
    [default = ../Inputs/eta_materials_compendium.csv]
source_path ('src'): str
    The name and path for the file containing the starting neutron source distribution. The format is a comma delimited 
    key word input file. All keywords are optional.  Non-specified keywords will default to preset program values. 
    [default = ../Inputs/source.csv]
restart ('r'): boolean
        Optional input to indicate the that search process will start with a preinitialized population.  All members of 
        the population must have full initialization inputs to work. 

"""

start_time=time.time()     #Start Clock
rundir=os.path.abspath(os.path.join(os.path.abspath(os.getcwd()),os.pardir))+"/Results/Population/"
    
# Set logging options
if os.path.exists('{}/Results'.format(os.path.abspath(os.path.join(os.path.abspath(os.getcwd()), os.pardir)))):    
    if os.path.isfile('{}/Results/logfile.log'.format(os.path.abspath(os.path.join(os.path.abspath(os.getcwd()), os.pardir)))):
        os.remove('{}/Results/logfile.log'.format(os.path.abspath(os.path.join(os.path.abspath(os.getcwd()), os.pardir))))
else:
    os.mkdir('{}/Results'.format(os.path.abspath(os.path.join(os.path.abspath(os.getcwd()), os.pardir))))
logger = logging.getLogger('Coeus')
fh = logging.FileHandler('{}/Results/logfile.log'.format(os.path.abspath(os.path.join(os.path.abspath(os.getcwd()), os.pardir))))
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
logger.addHandler(fh)
logger.setLevel(logging.INFO)
logger.info('Started Coeus:\n')

# Create the output folder
if os.path.exists('{}/Results/Population'.format(os.path.abspath(os.path.join(os.path.abspath(os.getcwd()), os.pardir))))==False:    
    os.mkdir('{}/Results/Population'.format(os.path.abspath(os.path.join(os.path.abspath(os.getcwd()), os.pardir))))
    
# Set print options to print full numpy arrays
np.set_printoptions(threshold='nan')
    
# Parse the input arguments
logger.info('Reading inputs and initializing settings:')
parser = argparse.ArgumentParser()
parser.add_argument('--r', nargs='?', default='n', help='Boolean indicator for if an initial population is supplied.  This initial population must be in the form of MCNP input decks in the Coeus standard directories.  Options are y or n.  [default = n]')
parser.add_argument('--obj', nargs='?', default=os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/Inputs/obj_spectrum.csv',help='The name and path for the objective function (spectra) file location. The format is a comma delimited spectrum with the first column being the upper bin boundaries and the second column the bin flux/fluence.  [default = ../Inputs/obj_spectrum.csv]')
parser.add_argument('--eta', nargs='?', default=os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/Inputs/eta_constraints.csv',help='The name and path for the file containing the ETA design constraints. The format is a comma delimited key word input file. All keywords are optional.  Non-specified keywords will default to preset program values. [default = ../Inputs/eta_constraints.csv]')
parser.add_argument('--gs', nargs='?', default=os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/Inputs/gnowee_settings.csv',help='The name and path for the file containing the Gnowee search settings. The format is a comma delimited key word input file. All keywords are optional.  Non-specified keywords will default to preset program values.   [default = ../Inputs/gnowee_settings.csv]')
parser.add_argument('--adv', nargs='?', default=os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/Inputs/advantg_settings.csv',help='The name and path for the file containing the advantg settings. The format is a comma delimited key word input file. All keywords are optional.  Non-specified keywords will default to preset program values. [default = ../Inputs/advantg_settings.csv]')
parser.add_argument('--mcnp', nargs='?', default=os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/Inputs/mcnp_settings.csv',help='The name and path for the file containing the mcnp settings. The format is a comma delimited key word input file. All keywords are optional.  Non-specified keywords will default to preset program values. [default = ../Inputs/mcnp_settings.csv]')
parser.add_argument('--mat', nargs='?', default=os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/Inputs/eta_materials_compendium.csv',help='The name and path for the file containing the materials to be included in the problem. The format is a comma delimited key word input file. [default = ../Inputs/eta_materials_compendium.csv]')
parser.add_argument('--src', nargs='?', default=os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/Inputs/source.csv',help='The name and path for the file containing the starting neutron source distribution. The format is a comma delimited key word input file. All keywords are optional. Non-specified keywords will default to preset program values. [default = ../Inputs/source.csv]')
parser.add_argument('--qos', nargs='?', default='nuclear_normal')
parser.add_argument('--account', nargs='?', default='co_nuclear')
parser.add_argument('--partition', nargs='?', default='savio')


args = parser.parse_args()    

# Assign optional inputs to variables:
if args.obj:
    obj_path=args.obj 
if args.eta:
    eta_path=args.eta
if args.gs:
    gs_path=args.gs
if args.adv:
    advantg_path=args.adv
if args.mcnp:
    mcnp_path=args.mcnp
if args.mat:
    materials_library_path=args.mat
if args.src:
    source_path=args.src

# Create ETA_Params object
eta_params=ETA_Parameters()
   
# Test path for objective function. Call read_obj if file exists
if os.path.isfile(obj_path): 
    logger.info("\nLoading objective spectra file located at: {}".format(obj_path))
    ETA_Parameters.read_obj(eta_params,obj_path)
else:
    logger.info("\nNo user supplier objective spectra file located.  Program default values to be used instead.")
    
# Test path for constraint file. Call read_constraint if file exists
if os.path.isfile(eta_path): 
    logger.info("\nLoading ETA constraints file located at: {}".format(eta_path))
    ETA_Parameters.read_constraints(eta_params,eta_path)
else:
    logger.info("\nNo user supplier ETA constraints file located.  Program default values to be used instead.")
    
    
# Create Gnowee Settings object
g_set=Gnowee_Settings()
    
# Test path for Gnowee settings file. Call read_settings if file exists
if os.path.isfile(gs_path): 
    logger.info("\nLoading Gnowee settings file located at: {}".format(gs_path))
    Gnowee_Settings.read_settings(g_set,gs_path)
else:
    logger.info("\nNo user supplier Gnowee Search settings file located.  Program default values to be used instead.")
        
# Create ADVANTG Settings object
advantg_set=ADVANTG_Settings()
   
# Test path for ADVANTG settings file. Call read_settings if file exists
if os.path.isfile(advantg_path): 
    logger.info("\nLoading ADVANTG settings file located at: {}".format(advantg_path))
    ADVANTG_Settings.read_settings(advantg_set, advantg_path)
else:
    logger.info("\nNo user supplier ADVANTG settings file located.  Program default values to be used instead.")

        
# Create MCNP Settings object
mcnp_set=MCNP_Settings(eta_params)
 
# Test path for MCNP settings file. Call read_settings if file exists
if os.path.isfile(mcnp_path): 
    logger.info("\nLoading MCNP settings file located at: {}".format(mcnp_path))
    MCNP_Settings.read_settings(mcnp_set, mcnp_path)
else:
    logger.info("\nNo user supplier MCNP settings file located.  Program default values to be used instead.")
    
# Test path for source file. Call read_source if file exists
if os.path.isfile(source_path): 
    logger.info("\nLoading source file located at: {}\n".format(source_path))
    MCNP_Settings.read_source(mcnp_set, source_path)
else:
    logger.info("\nNo user supplier source file located.  Program default values to be used instead.\n")       
    
# Build Materials Library
mat_lib=Build_Matlib(materials_library_path)
    
# Initialize the population
pop=[]
# Create baseline ETA geometry based on ETA constraints
base_eta=MCNP_Geometry()
base_eta.init_geom(eta_params, mat_lib)
if args.r=='y':
    for i in range(0,g_set.p):
        eta=MCNP_Geometry()
        nps=eta.read_geom(rundir+str(i)+"/ETA.inp", mat_lib)
        pop.append(Parent(i, eta_params, eta, g_set, mcnp_set, mat_lib, [eta_params.fissile_mat,'Au'], i,build_geom=False))
        pop[-1].rset.nps=nps
else:
    for i in range(0,g_set.p):
        pop.append(Parent(i, eta_params, base_eta, g_set, mcnp_set, mat_lib, [eta_params.fissile_mat,'Au'], i))
        pop[-1].geom.fin_geom(eta_params, mat_lib)
logger.info('Finished reading inputs and initializing settings in {} sec'.format(time.time() - start_time)) 

# Print MCNP input Files
ids=[]
particles=[]
for i in range(0,g_set.p):
    Print_MCNP_Input(eta_params,pop[i].geom,pop[i].rset,mat_lib,i,adv_print=False)
    ids.append(i)
    if args.r=='y':
        particles.append(pop[i].rset.nps)
    else:
        particles.append(mcnp_set.nps)
    
# Print ADVANTG input Files
for i in ids:
    Print_ADVANTG_Input(eta_params,pop[i].geom,advantg_set,i,cluster=True)
logger.info('Finished printing initial input files at {} sec\n'.format(time.time() - start_time))

# Run ADVANTG
Run_Transport(ids,code='advantg', qos=args.qos, account=args.account, partition=args.paritition)
logger.info('Finished running ADVANTG at {} sec\n'.format(time.time() - start_time))

# Run MCNP
for i in ids:
    Print_MCNP_Input(eta_params,pop[i].geom,pop[i].rset,mat_lib,i,adv_print=True)
Run_Transport(ids,particles,code='mcnp6.mpi', qos=args.qos, account=args.account, partition=args.paritition)
logger.info('Finished running MCNP at {} sec\n'.format(time.time() - start_time))

# Calculate Fitness
Calc_Fitness(ids, pop, eta_params.spectrum[:,1], eta_params.min_fiss, eta_params.max_weight)

# Save the output files
for c in pop:
    shutil.copyfile(rundir+str(c.ident)+'/tmp/ETA.out', rundir+str(c.ident)+'/ETA.out')
    
# Create and store first event in timeline and Meta_Stats
stats=Meta_Stats()
stats.write(header=True)
history=Timeline()
pop=history.update(pop, 1, g_set.p)
logger.info('Calculated fitness, saved files, and added to timeline at {} sec\n'.format(time.time() - start_time))

# Calculate Moderating ratios
mod_rat=Calc_Moderating_Ratio(mat_lib)

######## Partial Inversion ########
new_pop=Partial_Inversion(pop,mod_rat,mat_lib,g_set)
(ids,particles)=print_MCNP_input_files('Partial Inversion')
run_MCNP_on_algo(args,"part_inv", 0, int(g_set.p))

        
# Iterate until termination criterion met
converge = False
while history.tline[-1].g <= g_set.gm and history.tline[-1].e <= g_set.em and converge==False:

    logger.info('Generation {} with {} function evaluations completed started at {} sec\n'.format(history.tline[-1].g,history.tline[-1].e,time.time() - start_time))    
    
    ######## Levy flight permutation of materials ########
    new_pop=Mat_Levy_Flights(pop, mat_lib, mod_rat, g_set, [eta_params.fissile_mat,'Au'])
    (ids,particles)=print_MCNP_input_files("Levy flight permutation of materials")
    run_MCNP_on_algo(args,"mat_levy", 0, int(g_set.p*g_set.fl))

        
    ######## Levy flight permutation of cells ########
    new_pop=Cell_Levy_Flights(pop,eta_params,g_set)      
    (ids,particles)=print_MCNP_input_files("Levy flight permutation of cells")
    run_MCNP_on_algo(args,"cell_levy", 0, int(g_set.p*g_set.fl))

        
    ######## Elite_Crossover ########
    new_pop=Elite_Crossover(pop,mod_rat,eta_params,mat_lib,g_set,[eta_params.fissile_mat,'Au'])
    (ids,particles)=print_MCNP_input_files("Elite Crossover")
    run_MCNP_on_algo(args,"elite_cross", 0, 1)
            
            
    ######## Mutate ########
    new_pop=Mutate(pop, eta_params, g_set)
    (ids,particles)=print_MCNP_input_files("Mutation Operator")
    run_MCNP_on_algo(args,"mutate", 0, int(g_set.p))
      

    ######## Crossover ########
    new_pop=Crossover(pop,g_set)
    (ids,particles)=print_MCNP_input_files("Crossover")
    run_MCNP_on_algo(args,"crossover", 0, int(g_set.p*g_set.fe))
        
        
    ######## 2-opt ########
    if eta_params.max_horiz >= 4:
        new_pop=Two_opt(pop,g_set)
        (ids,particles)=print_MCNP_input_files("2-opt")
        run_MCNP_on_algo(args,"two_opt", 0, int(g_set.p*g_set.fe))
        
        
    ######## 3-opt ########
    if eta_params.max_horiz >= 6:
        new_pop=Three_opt(pop,g_set)
        (ids,particles)=print_MCNP_input_files("3-opt")
        run_MCNP_on_algo(args,"three_op", 0, int(g_set.p))
        
        
    ######## Discard Cells ########
    new_pop=Discard_Cells(pop,mat_lib,g_set)
    (ids,particles)=print_MCNP_input_files("Discard Cells")
    run_MCNP_on_algo(args,"discard", 1, int(g_set.p*g_set.fd))
    stats.write()
                
    ######## Test Convergence ########
    # Test generational convergence
    if history.tline[-1].g > g_set.sl:
        if history.tline[-1].g > history.tline[-2].g+g_set.sl:
            converge=True
            logger.info("Generational Stall {}".format(str(history.tline[-1])))
        elif (history.tline[-2].f-history.tline[-1].f)/history.tline[-2].f < g_set.ct:
            if history.tline[-1].g > history.tline[-2].g+g_set.sl:
                converge=True
                logger.info("Generational Convergence")
                
    # Test fitness convergence
    if abs((history.tline[-1].f-g_set.of)/g_set.of) <= g_set.ot:
        converge=True
        logger.info("Fitness Convergence {}".format(str(history.tline[-1])))        
        
    ######## Update weight window maps ########
    # Print MCNP input Files
    ids=[]
    for i in range(0,g_set.p):
        Print_MCNP_Input(eta_params,pop[i].geom,pop[i].rset,mat_lib,i,adv_print=False)
        ids.append(i)

    # Print ADVANTG input Files
    for i in ids:
        Print_ADVANTG_Input(eta_params,pop[i].geom,advantg_set,i,cluster=True)

    # Run ADVANTG
    Run_Transport(ids,code='advantg', qos=args.qos, account=args.account, partition=args.paritition)

#Determine execution time    
logger.info('The optimization history is:{}\n'.format(history.tline))
logger.info('Total run time was {} sec'.format(time.time() - start_time))
