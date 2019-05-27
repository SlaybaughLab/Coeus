"""!
@file Coeus.py
@package Coeus

@brief The main program for Coeus. Coeus is a fully functional and encapsulated 
design tool for nuclear design problems. Designs are optimized to achieve
design criteria using user defined contraint and objective functions using the 
Gnowee metaheuristic algorithm.  ADVANTG and MCNP are used to perform the
radiation transport calculations. 

Coeus is designed to run on a cluster with a job submission script.  Local runs
on a PC are not longer supported.

All inputs are optional.  The program will load run inputs in the following
order: \n 

    1) User specified path \n 
    2) Run directory defaults (Note: naming convention for files to be loaded
       from run directory must follow default naming convention shown in
       parameters)\n 
    3) Preset program defaults.  \n 

This order will be followed for each of the settings files so that some may be 
ommitted if desired. \n 

@param input_path('inp'):str [default = ../Inputs/user_inputs.txt] 
    The name and path for the user input file location. The format is a space 
    delimited key word file as specified in the UserInputs class.
@param eta_constraints_path('eta'):str 
    [default = ../Inputs/eta_constraints.csv] 
    The name and path for the file containing the ETA design constraints. The 
    format is a comma delimited key word input file. All keywords are optional.
    Non-specified keywords will default to preset program values. 
@param gnowee_settings_path('gs'):str [default = ../Inputs/gnowee_settings.csv] 
    The name and path for the file containing the Gnowee search settings. The 
    format is a comma delimited key word input file. All keywords are optional.  
    Non-specified keywords will default to preset program values. 
@param advantg_settings_path('adv'):str 
    [default = ../Inputs/advantg_settings.csv] 
    The name and path for the file containing the advantg settings. 
    The format is a comma delimited key word input file. All keywords are 
    optional.  Non-specified keywords will default to preset program values. 
@param mcnp_settings_path('mcnp'):str [default = ../Inputs/mcnp_settings.csv] 
    The name and path for the file containing the mcnp settings. 
    The format is a comma delimited key word input file. All keywords are 
    optional.  Non-specified keywords will default to preset program values.
@param material_lib_path('mat'):str 
    [default = ../Inputs/eta_materials_compendium.csv]
    The name and path for the file containing the materials to be included in 
    the problem. The format is a comma delimited key word input file.
@param source_path('src'):str [default = ../Inputs/source.csv]
    The name and path for the file containing the starting neutron source 
    distribution. The format is a comma delimited key word input file. All 
    keywords are optional.  Non-specified keywords will default to preset 
    program values. \n 
@param restart('r'):boolean
    Optional input to indicate the that search process will start with a 
    preinitialized population.  All members of  the population must have full 
    initialization inputs to work. 

@author James Bevins

@date 27May19
"""


import numpy as np

from ETA_Utilities import ETA_Parameters

from Gnowee_Utilities import Gnowee_Settings, Parent, Calc_Fitness, Pop_Update
from Gnowee_Utilities import Timeline

from Metaheuristics import Mat_Levy_Flights, Cell_Levy_Flights 
from Metaheuristics import Elite_Crossover, Partial_Inversion, Two_opt
from Metaheuristics import Crossover, Three_opt, Discard_Cells, Mutate

from ADVANTG_Utilities import ADVANTG_Settings, Print_ADVANTG_Input

# Delete in near future.  Maybe modify Build_Matlib for mat library
#from NuclearData import Build_Matlib #, Calc_Moderating_Ratio

from MCNP_Utilities import MCNP_Settings, MCNP_Geometry, Print_MCNP_Input

from Utilities import Run_Transport, Meta_Stats

from UserInputs import UserInputs

import time
import shutil
import logging
import os
import sys

import argparse


#-----------------------------------------------------------------------------#
# Local Function definitions

"""!
Prints the radiation transport input files given a set of Gnowee generated
new parameters.

@param step: \e string \n
	The current operator name. \n
@param radCode: \e string \n
	String indicating the name of the radiation transport code to be used. \n
"""
def print_transport_input(step, radCode='MCNP'):
    global logger, history, start_time, new_pop, eta_params, mat_lib, objFunc
    idents, run_particles=[],[]
	
	# Loop over updated population and print MCNP input files
    for i in range(0,len(new_pop)):
        if radCode == "MCNP":
    	    Print_MCNP_Input(eta_params, objFunc.objective,
							 new_pop[i].geom,new_pop[i].rset,
							 mat_lib, new_pop[i].ident, adv_print=True)
							 
		# Create inputs to the job scheduler to define run parameters
        idents.append(new_pop[i].ident)
        run_particles.append(new_pop[i].rset.nps)
		
		# Delete?
        for m in range(9, len(new_pop[i].geom.matls)):
            if new_pop[i].geom.matls[m]==eta_params.fissile_mat:
                print("Line 82 Coeus.py")
                sys.exit()
    logger.info('Gen {} {} finished at {} sec\n'.format(history.tline[-1].g,
	            step, time.time() - start_time))
    return idents, run_particles

"""!
Runs the transport code for each operator provided a set of population members
to be evaluated.

@param args: \e object \n
	User input arguments for the job scheduler needed for run_transport. \n
@param algo: \e string \n
	String indicating the name of the algorithm being used. \n
@param updateGen: \e integer \n
	Flag used to update the generation number. \n
@param updateFeval: \e integer \n
	Flag used to update the number of function evaluation. \n
@param objFunc: \e string \n
	The objective function used to calculate the fitness. \n
@param radCode: \e string \n
	String indicating the name of the radiation transport code to be used. \n
"""
def run_MCNP_on_algo(args, algo, updateGen, updateFeval, objFunc,
                     radCode='MCNP'):
    global stats, logger, history, ids, particles, pop, new_pop, eta_params
    global mat_lib, mcnp_set
    
    slurmArgs=[args.qos, args.account, args.partition, args.timeout]
    
    if len(ids)>0:
        Run_Transport(ids, *slurmArgs, nps=particles, code='mcnp6.mpi')
        logger.info('Finished running MCNP at {} sec\n'.format(time.time() -
                                                               start_time))
    
        # Calculate Fitness
        Calc_Fitness(ids, new_pop, objFunc, eta_params.min_fiss,
                     eta_params.max_weight)
        (changes,feval)=Pop_Update(pop, new_pop, mcnp_set.nps, slurmArgs,
                                  eta_params, mat_lib, Run_Transport, rr=False) 
        pop=history.update(pop, updateGen, updateFeval)
        stats.update(algo,(changes, updateFeval + feval))
    
#-----------------------------------------------------------------------------#

def main():
    global stats, logger, history, start_time, ids, particles, pop, new_pop
    global eta_params, mat_lib, mcnp_set, objFunc
    
    # Start Clock
    start_time=time.time()    
    
    # Set Run directory path    
    rundir=os.path.abspath(os.path.join(os.path.abspath(os.getcwd()),
                                        os.pardir))+'/Results/Population/'
        
    # Set logging options
    if os.path.exists('../Results'):
        if os.path.isfile('../Results/logfile.log'):
            os.remove('../Results/logfile.log')
    else:
        os.mkdir('../Results')
    logger = logging.getLogger('Coeus')
    fh = logging.FileHandler('../Results/logfile.log')
    formatter = logging.Formatter(
                       '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.setLevel(logging.INFO)
    logger.info('Started Coeus:\n')

    # Create the output folder
    if os.path.exists('../Results/Population')==False:    
        os.mkdir('../Results/Population')
        
    # Set print options to print full numpy arrays
    np.set_printoptions(threshold=sys.maxsize)
        
    # Parse the input arguments
    logger.info('Reading inputs and initializing settings:')
    parser = argparse.ArgumentParser()
    parser.add_argument('--r', nargs='?', default='n',
                        help='Boolean indicator for if an initial population is supplied.  This initial population must be in the form of MCNP input decks in the Coeus standard directories.  Options are y or n.  [default = n]')
    parser.add_argument('--inp', nargs='?', 
                        default='../Inputs/user_inputs.txt',
                        help='The name and path for the user inputs file location. The format is a space delimited file with keyword arguments.  For more details, see the UserInputs class.  [default = ../Inputs/user_inputs.txt]')
    parser.add_argument('--eta', nargs='?', 
                        default='../Inputs/eta_constraints.csv',
                        help='The name and path for the file containing the ETA design constraints. The format is a comma delimited key word input file. All keywords are optional.  Non-specified keywords will default to preset program values. [default = ../Inputs/eta_constraints.csv]')
    parser.add_argument('--gs', nargs='?', 
                        default='../Inputs/gnowee_settings.csv',
                        help='The name and path for the file containing the Gnowee search settings. The format is a comma delimited key word input file. All keywords are optional.  Non-specified keywords will default to preset program values.   [default = ../Inputs/gnowee_settings.csv]')
    parser.add_argument('--adv', nargs='?', 
                        default='../Inputs/advantg_settings.csv',
                        help='The name and path for the file containing the advantg settings. The format is a comma delimited key word input file. All keywords are optional.  Non-specified keywords will default to preset program values. [default = ../Inputs/advantg_settings.csv]')
    parser.add_argument('--mcnp', nargs='?', 
                        default='../Inputs/mcnp_settings.csv',
                        help='The name and path for the file containing the mcnp settings. The format is a comma delimited key word input file. All keywords are optional.  Non-specified keywords will default to preset program values. [default = ../Inputs/mcnp_settings.csv]')
    parser.add_argument('--mat', nargs='?', 
                        default='../Inputs/eta_materials_compendium.csv',
                        help='The name and path for the file containing the materials to be included in the problem. The format is a comma delimited key word input file. [default = ../Inputs/eta_materials_compendium.csv]')
    parser.add_argument('--src', nargs='?', 
                        default='../Inputs/source.csv',
                        help='The name and path for the file containing the starting neutron source distribution. The format is a comma delimited key word input file. All keywords are optional. Non-specified keywords will default to preset program values. [default = ../Inputs/source.csv]')
    parser.add_argument('--log', nargs='?', default='INFO',
                        help='Valid levels are "DEBUG", "INFO", "WARNING", "ERROR", and "CRITICAL" in ascending order.')
    parser.add_argument('--qos', nargs='?', default='savio_lowprio')
    parser.add_argument('--account', nargs='?', default='co_nuclear')
    parser.add_argument('--partition', nargs='?', default='savio')
    parser.add_argument('--timeout', nargs='?', default='02:30:00')

  

    # Assign optional inputs to variables:
    args = parser.parse_args()  

    # Modify logging level based on user input
    try:
        logger.setLevel(args.log.upper())
        logger.info('\nLogger set to {} level.'.format(
                                                  logger.getEffectiveLevel()))
    except:
        logger.info('\nNo valid logger level specifed. Deault "INFO" used.')
    
    # Create ETA_Params object
    eta_params=ETA_Parameters()
       
    # Test path for user input file.  Create the object if file exists.
    if os.path.isfile(args.inp): 
        logger.info("\nLoading input file located at: {}".format(args.inp))
        inputs = UserInputs(coeusInputPath=args.inp)
        objFunc = inputs.read_coeus_settings()
        logger.debug("{}".format(str(objFunc)))
    else:
        logger.info("\nNo user supplier input file located.  Program default \
                    values to be used instead.")
        
    # Test path for constraint file. Call read_constraint if file exists
    if os.path.isfile(args.eta): 
        logger.info("\nLoading ETA constraints file located at: {}".format(args.eta))
        ETA_Parameters.read_constraints(eta_params,args.eta)
    else:
        logger.info("\nNo user supplier ETA constraints file located.  Program default values to be used instead.")
        
        
    # Create Gnowee Settings object
    g_set=Gnowee_Settings()
        
    # Test path for Gnowee settings file. Call read_settings if file exists
    if os.path.isfile(args.gs): 
        logger.info("\nLoading Gnowee settings file located at: {}".format(args.gs))
        Gnowee_Settings.read_settings(g_set,args.gs)
    else:
        logger.info("\nNo user supplier Gnowee Search settings file located.  Program default values to be used instead.")
        
    # Create ADVANTG Settings object
    advantg_set=ADVANTG_Settings()
       
    # Test path for ADVANTG settings file. Call read_settings if file exists
    if os.path.isfile(args.adv): 
        logger.info("\nLoading ADVANTG settings file located at: {}".format(args.adv))
        ADVANTG_Settings.read_settings(advantg_set, args.adv)
    else:
        logger.info("\nNo user supplier ADVANTG settings file located.  Program default values to be used instead.")

            
    # Create MCNP Settings object
    mcnp_set=MCNP_Settings(eta_params)
     
    # Test path for MCNP settings file. Call read_settings if file exists
    if os.path.isfile(args.mcnp): 
        logger.info("\nLoading MCNP settings file located at: {}".format(args.mcnp))
        MCNP_Settings.read_settings(mcnp_set, args.mcnp)
    else:
        logger.info("\nNo user supplier MCNP settings file located.  Program default values to be used instead.")
        
    # Test path for source file. Call read_source if file exists
    if os.path.isfile(args.src): 
        logger.info("\nLoading source file located at: {}\n".format(args.src))
        MCNP_Settings.read_source(mcnp_set, args.src)
    else:
        logger.info("\nNo user supplier source file located.  Program default values to be used instead.\n")       
        
    # Build Materials Library
#    mat_lib=Build_Matlib(args.mat)
        
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
        Print_MCNP_Input(eta_params, objFunc.objective, pop[i].geom,pop[i].rset,mat_lib,i,adv_print=False)
        ids.append(i)
        if args.r=='y':
            particles.append(pop[i].rset.nps)
        else:
            particles.append(mcnp_set.nps)
        
    # Print ADVANTG input Files
    for i in ids:
        Print_ADVANTG_Input(eta_params,pop[i].geom,advantg_set,i,cluster=True)
    logger.info('Finished printing initial input files at {} sec\n'.format(time.time() - start_time))

    # Run ADVANTG (might need to comment out for mcnp)
    Run_Transport(ids,code='advantg', qos=args.qos, account=args.account, partition=args.partition, timeout=args.timeout)
    logger.info('Finished running ADVANTG at {} sec\n'.format(time.time() - start_time))

    # Run MCNP
    for i in ids:
        Print_MCNP_Input(eta_params, objFunc.objective, pop[i].geom,pop[i].rset,mat_lib,i,adv_print=True)
    Run_Transport(ids, args.qos, args.account, args.partition, args.timeout, nps=particles, code='mcnp6.mpi')
    logger.info('Finished running MCNP at {} sec\n'.format(time.time() - start_time))

    # Calculate Fitness
    Calc_Fitness(ids, pop, objFunc, eta_params.min_fiss, eta_params.max_weight)

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
    (ids,particles)=print_transport_input('Partial Inversion')
    run_MCNP_on_algo(args,"part_inv", 0, int(g_set.p), objFunc)

            
    # Iterate until termination criterion met
    converge = False
    while history.tline[-1].g <= g_set.gm and history.tline[-1].e <= g_set.em and converge==False:

        logger.info('Generation {} with {} function evaluations completed started at {} sec\n'.format(history.tline[-1].g,history.tline[-1].e,time.time() - start_time))    
        
        ######## Levy flight permutation of materials ########
        new_pop=Mat_Levy_Flights(pop, mat_lib, mod_rat, g_set, [eta_params.fissile_mat,'Au'])
        (ids,particles)=print_transport_input("Levy flight permutation of materials")
        run_MCNP_on_algo(args,"mat_levy", 0, int(g_set.p*g_set.fl), objFunc)

            
        ######## Levy flight permutation of cells ########
        new_pop=Cell_Levy_Flights(pop,eta_params,g_set)      
        (ids,particles)=print_transport_input("Levy flight permutation of cells")
        run_MCNP_on_algo(args,"cell_levy", 0, int(g_set.p*g_set.fl), objFunc)

            
        ######## Elite_Crossover ########
        new_pop=Elite_Crossover(pop,mod_rat,eta_params,mat_lib,g_set,[eta_params.fissile_mat,'Au'])
        (ids,particles)=print_transport_input("Elite Crossover")
        run_MCNP_on_algo(args,"elite_cross", 0, 1, objFunc)
                
                
        ######## Mutate ########
        new_pop=Mutate(pop, eta_params, g_set)
        (ids,particles)=print_transport_input("Mutation Operator")
        run_MCNP_on_algo(args,"mutate", 0, int(g_set.p), objFunc)
          

        ######## Crossover ########
        new_pop=Crossover(pop,g_set)
        (ids,particles)=print_transport_input("Crossover")
        run_MCNP_on_algo(args,"crossover", 0, int(g_set.p*g_set.fe), objFunc)
            
            
        ######## 2-opt ########
        if eta_params.max_horiz >= 4:
            new_pop=Two_opt(pop,g_set)
            (ids,particles)=print_transport_input("2-opt")
            run_MCNP_on_algo(args,"two_opt", 0, int(g_set.p*g_set.fe), objFunc)
            
            
        ######## 3-opt ########
        if eta_params.max_horiz >= 6:
            new_pop=Three_opt(pop,g_set)
            (ids,particles)=print_transport_input("3-opt")
            run_MCNP_on_algo(args,"three_op", 0, int(g_set.p), objFunc)
            
            
        ######## Discard Cells ########
        new_pop=Discard_Cells(pop,mat_lib,g_set)
        (ids,particles)=print_transport_input("Discard Cells")
        run_MCNP_on_algo(args,"discard", 1, int(g_set.p*g_set.fd), objFunc)
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
            Print_MCNP_Input(eta_params,objFunc.objective,pop[i].geom,pop[i].rset,mat_lib,i,adv_print=False)
            ids.append(i)

        # Print ADVANTG input Files
        for i in ids:
            Print_ADVANTG_Input(eta_params,pop[i].geom,advantg_set,i,cluster=True)

        # Run ADVANTG
        Run_Transport(ids,code='advantg', qos=args.qos, account=args.account, partition=args.partition, timeout=args.timeout)

    #Determine execution time    
    logger.info('The optimization history is:{}\n'.format(history.tline))
    logger.info('Total run time was {} sec'.format(time.time() - start_time))

if __name__ == "__main__":
    main()
