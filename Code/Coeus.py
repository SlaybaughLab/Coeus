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

@author James Bevins

@date 14Jun19
"""

import time
import shutil
import logging
import os
import sys
import argparse

import numpy as np

from ETA_Utilities import ETA_Parameters

from Gnowee_Utilities import Gnowee_Settings, Parent, Calc_Fitness, Pop_Update
from Gnowee_Utilities import Timeline

from Metaheuristics import Mat_Levy_Flights, Cell_Levy_Flights
from Metaheuristics import Elite_Crossover, Partial_Inversion, Two_opt
from Metaheuristics import Crossover, Three_opt, Discard_Cells, Mutate

from ADVANTG_Utilities import ADVANTG_Settings, Print_ADVANTG_Input

# Delete in near future.  Maybe modify Build_Matlib for mat library
from NuclearData import Build_Matlib, Calc_Moderating_Ratio

from MCNP_Utilities import MCNP_Settings, MCNP_Geometry, Print_MCNP_Input

from Utilities import Run_Transport, Meta_Stats

from UserInputs import UserInputs

#-----------------------------------------------------------------------------#
# Global variables initial definition
stats = Meta_Stats()
logger = logging.getLogger('Coeus')
history = Timeline()
startTime = time.time()
ids, particles, pop, newPop = [], [], [], []
etaParams = ETA_Parameters()
mcnpSet = MCNP_Settings()
matLib = Build_Matlib()

#-----------------------------------------------------------------------------#
# Local Function definitions
def _print_transport_input(step, objFunc, radCode='MCNP'):
    """!
    Prints the radiation transport input files given a set of Gnowee generated
    new parameters.

    @param step: \e string \n
    	The current operator name. \n
    @param objFunc: \e object \n
    	An ObjectiveFunction object used to calculate the fitness. \n
    @param radCode: \e string \n
    	Indicates the name of the radiation transport code to be used. \n
    """
    global logger, history, startTime, newPop, etaParams, newPop, matib
    idents, runParticles = [], []

	# Loop over updated population and print MCNP input files
    for i in range(0, len(newPop)):
        if radCode == "MCNP":
            Print_MCNP_Input(etaParams, objFunc.objective,
                            newPop[i].geom, newPop[i].rset,
							matLib, newPop[i].ident, advPrint=True)

		# Create inputs to the job scheduler to define run parameters
        idents.append(newPop[i].ident)
        runParticles.append(newPop[i].rset.nps)

    logger.info('Gen {} {} finished at {} sec\n'.format(history.tline[-1].g,
	            step, time.time() - startTime))
    return idents, runParticles

def _run_transport_on_algo(args, algo, updateGen, updateFeval, objFunc,
                     radCode='MCNP'):
    """!
    Runs the transport code for each operator provided a set of population
    members to be evaluated.

    @param args: \e object \n
    	User input arguments for the job scheduler needed for run_transport. \n
    @param algo: \e string \n
    	String indicating the name of the algorithm being used. \n
    @param updateGen: \e integer \n
    	Flag used to update the generation number. \n
    @param updateFeval: \e integer \n
    	Flag used to update the number of function evaluation. \n
    @param objFunc: \e object \n
    	An ObjectiveFunction object used to calculate the fitness. \n
    @param radCode: \e string \n
    	String indicating the name of the radiation transport code to be used.
        \n
    """
    global stats, history, logger, ids, particles, pop, newPop, etaParams
    global mcnpSet, matib

    slurmArgs = [args.qos, args.account, args.partition, args.timeout]

    if not ids:
        if radCode == "MCNP":
            Run_Transport(ids, *slurmArgs, nps=particles, code='mcnp6.mpi')
        logger.info('Finished running MCNP at {} sec\n'.format(time.time() -
                                                               startTime))

        # Calculate Fitness
        Calc_Fitness(ids, newPop, objFunc, etaParams.min_fiss,
                     etaParams.max_weight)
        (changes, feval) = Pop_Update(pop, newPop, mcnpSet.nps, slurmArgs,
                                  etaParams, matLib, Run_Transport, rr=False)
        pop = history.update(pop, updateGen, updateFeval)
        stats.update(algo, (changes, updateFeval + feval))

#-----------------------------------------------------------------------------#
def main():
    """!
    Runs the Coeus optimization algorithm.

    Need to add robust description.
    """
    global stats, logger, history, startTime, ids, particles, pop, newPop
    global etaParams, mcnpSet, matib

    # Set Run directory path
    rundir = os.path.abspath(os.path.join(os.path.abspath(os.getcwd()),
                                        os.pardir))+'/Results/Population/'

    # Set logging options
    if os.path.exists('../Results'):
        if os.path.isfile('../Results/logfile.log'):
            os.remove('../Results/logfile.log')
    else:
        os.mkdir('../Results')
    fh = logging.FileHandler('../Results/logfile.log')
    formatter = logging.Formatter(
                       '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.setLevel(logging.INFO)
    logger.info('Started Coeus:\n')

    # Create the output folder
    if not os.path.exists('../Results/Population'):
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
    parser.add_argument('--qos', nargs='?', default='savio_lowprio',
                        help='The Quality of Service (QOS) for a specified \
                             account.')
    parser.add_argument('--account', nargs='?', default='co_nuclear',
                        help='Job submission account for all of the slave \
                        jobs.')
    parser.add_argument('--partition', nargs='?', default='savio',
                        help='Job submission partition for all of the slave \
                        jobs.')
    parser.add_argument('--timeout', nargs='?', default='02:30:00',
                        help='Job timout for all of the slave jobs.')

    # Assign optional inputs to variables:
    args = parser.parse_args()

    # Modify logging level based on user input
    try:
        logger.setLevel(args.log.upper())
        logger.info('\nLogger set to {} level.'.format(
                                                  logger.getEffectiveLevel()))
    except ValueError:
        logger.info('\nNo valid logger level specifed. Deault "INFO" used.')

    # Test path for user input file.  Create the object if file exists.
    if os.path.isfile(args.inp):
        logger.info("\nLoading input file located at: {}".format(args.inp))
        inputs = UserInputs(coeusInputPath=args.inp)
        objFunc = inputs.read_coeus_settings()
    else:
        logger.info("\nNo user supplier input file located.  Program default \
                    values to be used instead.")

    # Test path for constraint file. Call read_constraint if file exists
    if os.path.isfile(args.eta):
        logger.info("\nLoading ETA constraints file located at: {}".format(
                                                                     args.eta))
        ETA_Parameters.read_constraints(etaParams, args.eta)
    else:
        logger.info("\nNo user supplier ETA constraints file located.  \
                    Program default values to be used instead.")

    # Create Gnowee Settings object
    gSet = Gnowee_Settings()

    # Test path for Gnowee settings file. Call read_settings if file exists
    if os.path.isfile(args.gs):
        logger.info("\nLoading Gnowee settings file located at: {}".format(
                                                                      args.gs))
        Gnowee_Settings.read_settings(gSet, args.gs)
    else:
        logger.info("\nNo user supplier Gnowee Search settings file located. \
                    Program default values to be used instead.")

    # Create ADVANTG Settings object
    advantgSet = ADVANTG_Settings()

    # Test path for ADVANTG settings file. Call read_settings if file exists
    if os.path.isfile(args.adv):
        logger.info("\nLoading ADVANTG settings file located at: {}".format(
                                                                     args.adv))
        ADVANTG_Settings.read_settings(advantgSet, args.adv)
    else:
        logger.info("\nNo user supplier ADVANTG settings file located.  \
                                   Program default values to be used instead.")


    # Create MCNP Settings object
    mcnpSet = MCNP_Settings(etaParams)

    # Test path for MCNP settings file. Call read_settings if file exists
    if os.path.isfile(args.mcnp):
        logger.info("\nLoading MCNP settings file located at: {}".format(
                                                                    args.mcnp))
        MCNP_Settings.read_settings(mcnpSet, args.mcnp)
    else:
        logger.info("\nNo user supplier MCNP settings file located.  \
                                   Program default values to be used instead.")

    # Test path for source file. Call read_source if file exists
    if os.path.isfile(args.src):
        logger.info("\nLoading source file located at: {}\n".format(args.src))
        MCNP_Settings.read_source(mcnpSet, args.src)
    else:
        logger.info("\nNo user supplier source file located.  \
                                 Program default values to be used instead.\n")

    # Build Materials Library
#    matLib = Build_Matlib(args.mat)

    # Create baseline ETA geometry based on ETA constraints
    baseEta = MCNP_Geometry()
    baseEta.init_geom(etaParams, matLib)
    if args.r == 'y':
        for i in range(0, gSet.p):
            eta = MCNP_Geometry()
            nps = eta.read_geom(rundir+str(i)+"/ETA.inp", matLib)
            pop.append(Parent(i, etaParams, eta, gSet, mcnpSet, matLib,
                              [etaParams.fissile_mat, 'Au'], i,
                              build_geom=False))
            pop[-1].rset.nps = nps
    else:
        for i in range(0, gSet.p):
            pop.append(Parent(i, etaParams, baseEta, gSet, mcnpSet,
                              matLib, [etaParams.fissile_mat, 'Au'], i))
            pop[-1].geom.fin_geom(etaParams, matLib)
    logger.info('Finished reading inputs and initializing settings in {} sec\
                '.format(time.time() - startTime))

    # Print initial MCNP input Files
    for i in range(0, gSet.p):
        Print_MCNP_Input(etaParams, objFunc.objective, pop[i].geom,
                         pop[i].rset, matLib, i, advPrint=False)
        ids.append(i)
        if args.r == 'y':
            particles.append(pop[i].rset.nps)
        else:
            particles.append(mcnpSet.nps)

    # Print ADVANTG input Files
    for i in ids:
        Print_ADVANTG_Input(etaParams, pop[i].geom, advantgSet, i,
                            cluster=True)
    logger.info('Finished printing initial input files at {} sec\n'.format(
                                                     time.time() - startTime))

    # Run ADVANTG
    Run_Transport(ids, code='advantg', qos=args.qos, account=args.account,
                  partition=args.partition, timeout=args.timeout)
    logger.info('Finished running ADVANTG at {} sec\n'.format(time.time() -
                                                              startTime))

    # Run MCNP
    for i in ids:
        Print_MCNP_Input(etaParams, objFunc.objective, pop[i].geom,
                         pop[i].rset, matLib, i, advPrint=True)
    Run_Transport(ids, args.qos, args.account, args.partition, args.timeout,
                  nps=particles, code='mcnp6.mpi')
    logger.info('Finished running MCNP at {} sec\n'.format(time.time() -
                                                          startTime))

    # Calculate Fitness
    Calc_Fitness(ids, pop, objFunc, etaParams.min_fiss, etaParams.max_weight)

    # Save the output files
    for c in pop:
        shutil.copyfile(rundir+str(c.ident)+'/tmp/ETA.out',
                        rundir+str(c.ident)+'/ETA.out')

    # Create and store first event in timeline and Meta_Stats
    stats.write(header=True)
    pop = history.update(pop, 1, gSet.p)
    logger.info('Calculated fitness, saved files, and added to timeline at {} \
                sec\n'.format(time.time() - startTime))

    #! Modify to remove PyNE dependence
    # Calculate Moderating ratios
#    modRat = Calc_Moderating_Ratio(matLib)

    ######## Partial Inversion ########
#    newPop = Partial_Inversion(pop, mod_rat, matLib, gSet)
    (ids, particles) = _print_transport_input('Partial Inversion', objFunc)
    _run_transport_on_algo(args, "part_inv", 0, int(gSet.p), objFunc)

    # Iterate until termination criterion met
    converge = False
    while history.tline[-1].g <= gSet.gm and history.tline[-1].e <= gSet.em \
          and not converge:

        logger.info('Generation {} with {} function evaluations completed \
                    started at {} sec\n'.format(history.tline[-1].g,
                                                history.tline[-1].e,
                                                time.time() - startTime))

        ######## Levy flight permutation of materials ########
        newPop = Mat_Levy_Flights(pop, matLib, modRat, gSet,
                                 [etaParams.fissile_mat, 'Au'])
        (ids, particles) = _print_transport_input("Levy flight permutation of \
                                              materials", objFunc)
        _run_transport_on_algo(args, "mat_levy", 0, int(gSet.p*gSet.fl),
                               objFunc)

        ######## Levy flight permutation of cells ########
        newPop = Cell_Levy_Flights(pop, etaParams, gSet)
        (ids, particles) = _print_transport_input("Levy flight permutation of \
                                              cells", objFunc)
        _run_transport_on_algo(args, "cell_levy", 0, int(gSet.p*gSet.fl),
                               objFunc)

        ######## Elite_Crossover ########
        newPop = Elite_Crossover(pop, modRat, etaParams, matLib, gSet,
                                [etaParams.fissile_mat, 'Au'])
        (ids, particles) = _print_transport_input("Elite Crossover", objFunc)
        _run_transport_on_algo(args, "elite_cross", 0, 1, objFunc)

        ######## Mutate ########
        newPop = Mutate(pop, etaParams, gSet)
        (ids, particles) = _print_transport_input("Mutation Operator", objFunc)
        _run_transport_on_algo(args, "mutate", 0, int(gSet.p), objFunc)

        ######## Crossover ########
        newPop = Crossover(pop, gSet)
        (ids, particles) = _print_transport_input("Crossover", objFunc)
        _run_transport_on_algo(args, "crossover", 0, int(gSet.p*gSet.fe),
                               objFunc)

        ######## 2-opt ########
        if etaParams.max_horiz >= 4:
            newPop = Two_opt(pop, gSet)
            (ids, particles) = _print_transport_input("2-opt", objFunc)
            _run_transport_on_algo(args, "two_opt", 0, int(gSet.p*gSet.fe),
                             objFunc)

        ######## 3-opt ########
        if etaParams.max_horiz >= 6:
            newPop = Three_opt(pop, gSet)
            (ids, particles) = _print_transport_input("3-opt", objFunc)
            _run_transport_on_algo(args, "three_op", 0, int(gSet.p), objFunc)

        ######## Discard Cells ########
        newPop = Discard_Cells(pop, matLib, gSet)
        (ids, particles) = _print_transport_input("Discard Cells", objFunc)
        _run_transport_on_algo(args, "discard", 1, int(gSet.p*gSet.fd), objFunc)
        stats.write()

        ######## Test Convergence ########
        # Test generational convergence
        if history.tline[-1].g > gSet.sl:
            if history.tline[-1].g > history.tline[-2].g+gSet.sl:
                converge = True
                logger.info("Generational Stall {}".format(
                                                       str(history.tline[-1])))
            elif (history.tline[-2].f-history.tline[-1].f)/ \
                 history.tline[-2].f < gSet.ct:
                if history.tline[-1].g > history.tline[-2].g+gSet.sl:
                    converge = True
                    logger.info("Generational Convergence")

        # Test fitness convergence
        if abs((history.tline[-1].f-gSet.of)/gSet.of) <= gSet.ot:
            converge = True
            logger.info("Fitness Convergence {}".format(str(history.tline[-1])))

        ######## Update weight window maps ########
        # Print MCNP input Files
        ids = []
        for i in range(0, gSet.p):
            Print_MCNP_Input(etaParams, objFunc.objective,
                             pop[i].geom, pop[i].rset, matLib, i,
                             advPrint=False)
            ids.append(i)

        # Print ADVANTG input Files
        for i in ids:
            Print_ADVANTG_Input(etaParams, pop[i].geom, advantgSet, i,
                                cluster=True)

        # Run ADVANTG
        Run_Transport(ids, code='advantg', qos=args.qos, account=args.account,
                      partition=args.partition, timeout=args.timeout)

    #Determine execution time
    logger.info('The optimization history is:{}\n'.format(history.tline))
    logger.info('Total run time was {} sec'.format(time.time() - startTime))

if __name__ == "__main__":
    main()
