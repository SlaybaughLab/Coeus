#!/usr/bin/env python

#######################################################################################################
#
# Program : Coeus_local.py.
#
# Contains : The main program for Coeus. Coeus is a fully functional and encapsulated design tool 
#            for Energy Tuning Assemblies.  Designs are optimized to set design criteria using contraints
#            and spectrum objective functions using the Gnowee metaheuristic algorithm.  ADVANTG and MCNP
#            are used to perform the radiation transport calculations.  
#
#            This version runs on single cores to enable desktop testing of the functions and methods. 
#            It is not useful to finalize an ETA design. 
#
# Author : James Bevins
#
# Last Modified: 24Oct16
#
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
from Utilities import Run_Transport_Threads, Event, Meta_Stats

import time
import shutil
import logging

import os
import sys
sys.path.insert(0,os.path.abspath(os.getcwd())+'/Sampling')
 

def main(**kwargs):
    
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
    def run_MCNP_on_algo(algo, update_gen, update_feval):
        global ids, particles, pop

        if len(ids)>0:
            Run_Transport_Threads(ids,code='mcnp6')
            logger.info('Finished running MCNP at {} sec\n'.format(time.time() - start_time))

            # Calculate Fitness
            Calc_Fitness(ids, new_pop, eta_params.spectrum[:,1], eta_params.min_fiss, eta_params.max_weight)
            (changes,feval)=Pop_Update(pop, new_pop, mcnp_set.nps, eta_params, mat_lib, Run_Transport, rr=False) 
            pop=history.update(pop, update_gen, update_feval)
            stats.update(algo,(changes, update_feval + feval)) 
#-------------------------------------------------------------------------------------------------------------# 

    """
    Entry point for the Coeus program.  
    
    All inputs are optional.  The program will load run inputs in the following order:
    
    1) User specified path
    
    2) Run directory defaults (Note: naming convention for files to be loaded from run directory must follow
    default naming convention shown in parameters)
    
    3) Preset program defaults.  
    
    This order will be followed for each of the settings files so that some may be ommitted if desired. 
    
    Parameters
    ==========
    restart : boolean
        Idicator whether program is to restart from a set population
    obj_path ('o' or 'obj' for short): str
        The name and path for the objective function (spectra) file location. The format is a comma delimited spectrum
        with the first column being the upper bin boundaries and the second column the bin flux/fluence.
        [default = ../Inputs/obj_spectrum.csv]
    eta_constraints_path ('e' or 'eta' for short): str
        The name and path for the file containing the ETA design constraints. The format is a comma delimited key word input 
        file. All keywords are optional.  Non-specified keywords will default to preset program values. 
        [default = ../Inputs/eta_constraints.csv]
    gnowee_settings_path ('g' or 'gs' for short): str
        The name and path for the file containing the Gnowee settings. The format is a comma delimited key word 
        input file. All keywords are optional.  Non-specified keywords will default to preset program values. 
        [default = ../Inputs/gnowee_settings.csv]
    advantg_settings_path ('a' or 'advantg' for short): str
        The name and path for the file containing the partisn settings. The format is a comma delimited key word input file.
        All keywords are optional.  Non-specified keywords will default to preset program values. 
        [default = ../Inputs/partisn_settings.csv]
    mcnp_settings_path ('mc' or 'mcnp' for short): str
        The name and path for the file containing the mcnp settings. The format is a comma delimited key word input file.
        All keywords are optional.  Non-specified keywords will default to preset program values. 
        [default = ../Inputs/mcnp_settings.csv]
    material_lib_path ('m' or 'material' for short): str
        The name and path for the file containing the starting neutron source distribution. The format is a comma delimited 
        key word input file. All keywords are optional.  Non-specified keywords will default to preset program values. 
        [default = ../Inputs/eta_materials_compendium.csv]
    source_path ('s' or 'source' for short): str
        The name and path for the file containing the starting neutron source distribution. The format is a comma delimited 
        key word input file. All keywords are optional.  Non-specified keywords will default to preset program values. 
        [default = ../Inputs/source.csv]
    restart: boolean
        Optional input to indicate the that search process will start with a preinitialized population.  All members of 
        the population must have full initialization inputs to work. 
    
    """
    
    start_time=time.time()     #Start Clock
    rundir=os.path.abspath(os.path.join(os.path.abspath(os.getcwd()),os.pardir))+"/Results/Population/"
    
    # Set logging options
    if os.path.isfile('{}/Results/logfile.log'.format(os.path.abspath(os.path.join(os.path.abspath(os.getcwd()), os.pardir)))):
        os.remove('{}/Results/logfile.log'.format(os.path.abspath(os.path.join(os.path.abspath(os.getcwd()), os.pardir))))
    logger = logging.getLogger('Coeus')
    hdlr = logging.FileHandler('{}/Results/logfile.log'.format(os.path.abspath(os.path.join(os.path.abspath(os.getcwd()), os.pardir))))
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr) 
    logger.setLevel(logging.INFO)
    logger.info('Started Coeus\n')

    # Set print options to print full numpy arrays
    np.set_printoptions(threshold='nan')
    
    # Read optional inputs:
    logger.info('Reading inputs and initializing settings:')
    if 'restart' in kwargs:
        restart = kwargs.get('restart')
    else :
        restart=False
    
    if 'obj_path' in kwargs:
        obj_path = kwargs.get('obj_path')
    elif "obj" in kwargs:
        obj_path = kwargs.get('obj')
    elif "o" in kwargs:
        obj_path = kwargs.get('o')
    else:
        obj_path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/Inputs/obj_spectrum.csv'   
        
    if 'eta_constraints_path' in kwargs:
        eta_path = kwargs.get('eta_constraints_path')
    elif "eta" in kwargs:
        eta_path = kwargs.get('eta')
    elif "e" in kwargs:
        eta_path = kwargs.get('e')
    else:
        eta_path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/Inputs/eta_constraints.csv' 
        
    if 'gnowee_settings_path' in kwargs:
        gs_path = kwargs.get('gs_settings_path')
    elif "gs" in kwargs:
        gs_path = kwargs.get('gs')
    elif "g" in kwargs:
        gs_path = kwargs.get('g')
    else:
        gs_path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/Inputs/gnowee_settings.csv' 
        
    if 'advantg_settings_path' in kwargs:
        advantg_path = kwargs.get('advantg_settings_path')
    elif "advantg" in kwargs:
        advantg_path = kwargs.get('advantg')
    elif "a" in kwargs:
        advantg_path = kwargs.get('a')
    else:
        advantg_path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/Inputs/advantg_settings.csv' 
        
    if 'mcnp_settings_path' in kwargs:
        mcnp_path = kwargs.get('mcnp_settings_path')
    elif "mcnp" in kwargs:
        mcnp_path = kwargs.get('mcnp')
    elif "mc" in kwargs:
        mcnp_path = kwargs.get('mc')
    else:
        mcnp_path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/Inputs/mcnp_settings.csv' 
        
    if 'materials_library_path' in kwargs:
        materials_library_path = kwargs.get('materials_library_path')
    elif "materials" in kwargs:
        materials_library_path = kwargs.get('materials')
    elif "m" in kwargs:
        materials_library_path = kwargs.get('m')
    else:
        materials_library_path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/Inputs/eta_materials_compendium.csv' 
        
    if 'source_path' in kwargs:
        source_path = kwargs.get('source_path')
    elif "source" in kwargs:
        source_path = kwargs.get('source')
    elif "s" in kwargs:
        source_path = kwargs.get('s')
    else:
        source_path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/Inputs/source.csv' 
    
    
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
    
    # Test path for gnowee settings file. Call read_settings if file exists
    if os.path.isfile(gs_path): 
        logger.info("\nLoading Gnowee settings file located at: {}".format(gs_path))
        Gnowee_Settings.read_settings(g_set,gs_path)
    else:
        logger.info("\nNo user supplier Gnowee settings file located.  Program default values to be used instead.")

        
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
        logger.info("\nLoading source file located at: {}".format(source_path))
        MCNP_Settings.read_source(mcnp_set, source_path)
    else:
        logger.info("\nNo user supplier source file located.  Program default values to be used instead.")  
    
    
    # Build Materials Library
    mat_lib=Build_Matlib(materials_library_path)
    
    
    # Initialize the population
    pop=[]
    # Create baseline ETA geometry based on ETA constraints
    base_eta=MCNP_Geometry()
    base_eta.init_geom(eta_params, mat_lib)
    if restart:
        for i in range(0,g_set.p):
            eta=MCNP_Geometry()
            nps=eta.read_geom(rundir+str(i)+"/ETA.inp", mat_lib)
            pop.append(Parent(i, eta_params, eta, g_set, mcnp_set, mat_lib, [eta_params.fissile_mat,'Au'], i,build_geom=False))
            pop[-1].rset.nps=nps
    else:
        for i in range(0,g_set.p):
            pop.append(Parent(i, eta_params, base_eta, g_set, mcnp_set, mat_lib, [eta_params.fissile_mat,'Au'], i))
            pop[-1].geom.fin_geom(eta_params, mat_lib)
    logger.info('Finished reading inputs and initializing settings at {} sec\n'.format(time.time() - start_time)) 
    
    # Print MCNP input Files
    ids=[]
    particles=[]
    for i in range(0,g_set.p):
        Print_MCNP_Input(eta_params,pop[i].geom,pop[i].rset,mat_lib,i,adv_print=False)
        ids.append(i)
        if restart:
            particles.append(pop[i].rset.nps)
        else:
            particles.append(mcnp_set.nps)
        
    # Print ADVANTG input Files
    for i in ids:
        Print_ADVANTG_Input(eta_params,pop[i].geom,advantg_set,i)
    logger.info('Finished printing initial files at {} sec\n'.format(time.time() - start_time)) 
    
    # Run ADVANTG
#    Run_Transport_Threads(ids,code='advantg')
    logger.info('Finished running ADVANTG at {} sec\n'.format(time.time() - start_time))
    
    # Run MCNP
    for i in ids:
        Print_MCNP_Input(eta_params,pop[i].geom,pop[i].rset,mat_lib,i,adv_print=True)
    Run_Transport_Threads(ids,code='mcnp6')
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
    run_MCNP_on_algo("part_inv", 0, int(g_set.p))
    
    # Iterate until termination criterion met
    converge = False
    while history.tline[-1].g <= g_set.gm and history.tline[-1].e <= g_set.em and converge==False:
        logger.info('Generation {} started at {} sec\n'.format(history.tline[-1].g,time.time() - start_time))
        
        ######## Levy flight permutation of materials ########
        new_pop=Mat_Levy_Flights(pop, mat_lib, mod_rat, g_set, [eta_params.fissile_mat,'Au'])
        (ids,particles)=print_MCNP_input_files( "Levy flight permutation of materials")
        run_MCNP_on_algo("mat_levy", 0, int(g_set.p*g_set.fl))

        
        ######## Levy flight permutation of cells ########
        new_pop=Cell_Levy_Flights(pop,eta_params,g_set)
        (ids,particles)=print_MCNP_input_files("Levy flight permutation of cells")
        run_MCNP_on_algo("cell_levy", 0, int(g_set.p*g_set.fl))


        ######## Elite_Crossover ########
        new_pop=Elite_Crossover(pop,mod_rat,eta_params,mat_lib,g_set,[eta_params.fissile_mat,'Au'])
        (ids,particles)=print_MCNP_input_files("Elite_Crossover")
        run_MCNP_on_algo("elite_cross", 0, 1) 
        
        
        ######## Mutate ########
        new_pop=Mutate(pop, eta_params, g_set)
        (ids,particles)=print_MCNP_input_files("Mutation Operator")
        run_MCNP_on_algo("mutate", 0, int(g_set.p))
        
        
        ######## Crossover ########
        new_pop=Crossover(pop,g_set)
        (ids,particles)=print_MCNP_input_files("Crossover")
        run_MCNP_on_algo("crossover", 0, int(g_set.p*g_set.fe))
        

        ######## 2-opt ########
        if eta_params.max_horiz >= 4:
            new_pop=Two_opt(pop,g_set)
            (ids,particles)=print_MCNP_input_files("2-opt")
            run_MCNP_on_algo("two_opt", 0, int(g_set.p*g_set.fe))
        
        
        ######## 3-opt ########
        if eta_params.max_horiz >= 6:
            new_pop=Three_opt(pop,g_set)
            (ids,particles)=print_MCNP_input_files("3-opt")
            run_MCNP_on_algo("three_op", 0, int(g_set.p))


        ######## Discard Cells ########
        new_pop=Discard_Cells(pop,mat_lib,g_set)
        (ids,particles)=print_MCNP_input_files("Discard Cells")
        run_MCNP_on_algo("discard", 1, int(g_set.p*g_set.fd))
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
    
        ###########################################33
        #if history.tline[-1].g>10:
        #converge=True
        ###########################################33
        
        ######## Update weight window maps ########
        # Print MCNP input Files
        ids=[]
        for i in range(0,g_set.p):
            Print_MCNP_Input(eta_params,pop[i].geom,pop[i].rset,mat_lib,i,adv_print=False)
            ids.append(i)

        # Print ADVANTG input Files
        for i in ids:
            Print_ADVANTG_Input(eta_params,pop[i].geom,advantg_set,i)

        # Run ADVANTG
#        Run_Transport_Threads(ids,code='advantg')

    #Determine execution time
    logger.info('The optimization history is:{}\n\n'.format(history.tline))
    logger.info('Total run time was {} sec'.format(time.time() - start_time))   
        
if __name__ == '__main__':
    main()
