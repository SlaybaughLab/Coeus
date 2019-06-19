"""!
@file Utilities.py
@package Coeus

@defgroup Utilities Utilities

@brief Utility functions for the Coeus program.

@author James Bevins

@date 18June19
"""

import os
import time
import bisect
import logging

import multiprocessing as mp
import subprocess as sub
import numpy as np

module_logger = logging.getLogger('Coeus.Utilities')

#-----------------------------------------------------------------------------#
class Switch(object):
    """!
    @ingroup Utilities
    Creates a switch class object to switch between cases. \n
    """

    def __init__(self, value):
        """!
        Constructor to build the Switch class.

        @param self: \e pointer \n
            The Switch class pointer. \n
        @param value: \e string \n
            The selector value. \n
        """

        ## @var value \e string Case selector value.
        self.value = value
        ## @var fall \e boolean Truth case based on match.
        self.fall = False

    def __iter__(self):
        """!
        Return the match method once, then stop.

        @param self: \e pointer \n
            The Switch class pointer. \n
        """
        yield self.match
        raise StopIteration

    def match(self, *args):
        """!
        Indicate whether or not to enter a case suite.

        @param args: \e list \n
            List of arguments to match with. \n
        """
        if self.fall or not args:
            return True
        elif self.value in args:
            self.fall = True
            return True
        else:
            return False

#-----------------------------------------------------------------------------#
class Event(object):
    """!
    @ingroup Utilities
    An event object representing a snapshot in the optimization process.
    """

    # Creates an event object representing a snapshot in the optimization
    #process
    # @return None
    def __init__(self, generation, evaluations, fitness, nps, ident):
        """!
        Constructor to build the Event class.

        @param self: \e pointer \n
            The Event pointer. \n
        @param generation: \e integer \n
            The generation the design was arrived at. \n
        @param evaluations: \e integer \n
            The number of fitness evaluations done to obtain this design. \n
        @param fitness: \e float \n
           he assessed design fitness. \n
        @param nps: \e float \n
            The number of particles run for that event. \n
        @param ident: \e integer \n
            The identifty of the current top solution. \n
        """

        assert generation >= 0, 'The # of generations cannot be negative.'
        assert evaluations >= 0, 'The # of evaluations cannot be negative.'
        assert isinstance(generation, int), 'Generation must be of type int.'
        assert isinstance(evaluations, int), 'Evaluations must be of type int.'
        assert isinstance(fitness, float), 'Fitness must be of type float.'

        ## @var generation \e integer The generation the design was arrived at.
        self.g = generation
        ## @var evaluations \e integer The number of fitness evaluations done
        ## to obtain this design.
        self.e = evaluations
        ## @var fitness \e float The assessed design fitness.
        self.f = fitness
        ## @var nps \e float The number of particles run for that event.
        self.n = nps
        ## @var ident  \e integer The identifty of the current top solution.
        self.i = ident

    def __repr__(self):
        """!
        Event class param print function.

        @param self: \e pointer \n
            The Event pointer. \n
        """
        return "Event instance({0}, {1}, {2}, {3}, {4})".format(self.g, self.e,
                                                                self.f, self.n,
                                                                self.i)

    def __str__(self):
        """!
        Human readable Event print function.

        @param self: \e pointer \n
            The Event pointer. \n
        """
        header = ["\nEvent Instance:"]
        header += ["Generation # = {}".format(self.g)]
        header += ["Number of Evaluations = {}".format(self.e)]
        header += ["Fitness = {}".format(self.f)]
        header += ["NPS = {}".format(self.n)]
        header += ["Identity = {}".format(self.i)]
        header = "\n".join(header)+"\n"
        return header

#-----------------------------------------------------------------------------#
class WeightedRandomGenerator(object):
    """!
    @ingroup Utilities
    Defines a class of weights to be used to select number of instances in
    array randomly with linear weighting.
    """

    def __init__(self, weights):
        """!
        Constructor to build the Event class.

        @param self: \e pointer \n
            The WeightedRandomGenerator pointer. \n
        @param weights: \e array \n
            The array of weights (Higher = more likely to be selected). \n
        """

        self.totals = []
        running_total = 0

        for w in weights:
            running_total += w
            self.totals.append(running_total)

    # @return  The randomly selected index of the weights array
    def next(self):
        """!
        @param self: \e pointer \n
            The WeightedRandomGenerator pointer. \n

        @return \e array The next randomly selected index of the weights array.
             \n
        """
        rnd = np.random.rand() * self.totals[-1]
        return bisect.bisect_right(self.totals, rnd)

    def __call__(self):
        """!
        @param self: \e pointer \n
            The WeightedRandomGenerator pointer. \n

        @return \e array The next randomly selected index of the weights array.
             \n
        """
        return self.next()

#-----------------------------------------------------------------------------#
class MetaStats(object):
    """!
    @ingroup Utilities
    Stores and prints effectiveness stats for each metaheuristic search method.
    """

    def __init__(self, mat_levy=(0, 0), cell_levy=(0, 0), elite_cross=(0, 0),
                 part_inv=(0, 0), mutate=(0, 0), two_opt=(0, 0),
                 crossover=(0, 0), three_op=(0, 0), discard=(0, 0),
                 fname="../Results/meta_stats.txt"):
        """!
        Constructor to build the MetaStats class.

        @param self: \e pointer \n
            The MetaStats pointer. \n
        @param mat_levy: \e tuple \n
            Contains the changes and total number of function evaluations for
            the Mat_Levy_Flights function. \n
        @param cell_levy: \e tuple \n
            Contains the changes and total number of function evaluations for
            the Cell_Levy_Flights function. \n
        @param elite_cross: \e tuple \n
            Contains the changes and total number of function evaluations for
            the Mutate_Mats function. \n
        @param part_inv: \e tuple \n
            Contains the changes and total number of function evaluations for
            the Partial_Inversion function. \n
        @param mutate: \e tuple \n
            Contains the changes and total number of function evaluations for
            the Mutate function. \n
        @param two_opt: \e tuple \n
            Contains the changes and total number of function evaluations for
            the 2-opt function. \n
        @param crossover: \e tuple \n
            Contains the changes and total number of function evaluations for
            the Crossover function. \n
        @param three_op: \e tuple \n
            Contains the changes and total number of function evaluations for
            the Three_opt function. \n
        @param discard: \e tuple \n
            Contains the changes and total number of function evaluations for
            the Discard function. \n
        @param fname: \e string \n
            The file name for the output results. \n
        """

        ## @var \e dictionary String name of each algorithms
        self.algorithms = {"mat_levy": mat_levy, "cell_levy": cell_levy,
                    "elite_cross": elite_cross, "part_inv":part_inv,
                    "mutate": mutate, "two_opt":two_opt,
                    "crossover":crossover, "three_op": three_op,
                    "discard": discard}
        ## @var \e string Name and path of the file to store the timeline for
        ## post processing
        self.fname = fname

        if os.path.isfile(self.fname):
            os.remove(self.fname)

    def __repr__(self):
        """!
        MetaStats class param print function.

        @param self: \e pointer \n
            The MetaStats pointer. \n
        """
        return '{:8d},{};  {:14d},{};  {:11d},{};  {:12d},{};  {:12d},{}; ' \
                '{:12d},{};  {:11d},{};  {:9d},{};  {:6d},{};\n'.format(
                self.algorithms["mat_levy"][0], self.algorithms["mat_levy"][1],
                self.algorithms["cell_levy"][0],
                self.algorithms["cell_levy"][1],
                self.algorithms["elite_cross"][0],
                self.algorithms["elite_cross"][1],
                self.algorithms["part_inv"][0], self.algorithms["part_inv"][1],
                self.algorithms["mutate"][0], self.algorithms["mutate"][1],
                self.algorithms["two_opt"][0], self.algorithms["two_opt"][1],
                self.algorithms["crossover"][0],
                self.algorithms["crossover"][1],
                self.algorithms["three_op"][0], self.algorithms["three_op"][1],
                self.algorithms["discard"][0], self.algorithms["discard"][1])

    def __str__(self):
        """!
        Human readable MetaStats print function.

        @param self: \e pointer \n
            The MetaStats pointer. \n
        """
        header = ["Mat_Levy_Flights  Cell_Levy_Flights  Elite_Crossover  \
                  Partial_Inversion  Mutate  Two_opt  Crossover  Three_opt  \
                  Discard"]
        header += repr(self)
        header = "\n".join(header)+"\n"
        s = header
        return s

    def update(self, alg, val):
        """!
        Adds val tuples to the algorithm arg's tuples.

        @param self: \e pointer \n
            The MetaStats object pointer. \n
        @param alg: \e string \n
            The name of the algorithm selected. \n
        @param val: \e tuple \n
            The tuple value to be added. \n
        """

        self.algorithms[alg] = (self.algorithms[alg][0]+val[0],
                                self.algorithms[alg][1]+val[1])

    def write(self, header=False):
        """!
        Adds val tuples to the algorithm arg's tuples.

        @param header: \e boolean \n
            Print header if true. \n
        """

        try:
            with open(self.fname, "a") as f:
                if header:
                    f.write(str(self))
                else:
                    f.write(repr(self))

            f.close()

        except IOError as e:
            module_logger.error("I/O error({0}): {1}".format(e.errno,
                                                             e.strerror))

        # Test that the file closed
        assert f.closed, "File did not close properly."

#-----------------------------------------------------------------------------#
def run_transport(lst, batchArgs, nps=[], code='mcnp6'):
    """!
    Build a Slurm Batch script using the Jobs Array feature to run transport
    calculations.

    @param lst: \e list \n
        List of parent identifier numbers to be ran. \n
    @param batchArgs: \e list \n
        List of batch job submission args. \n
    @param nps: \e list \n
        List of number of particles to run per code thread instance. If left
        blank, calculation will be performed to assign all availiable cpus
        evenly. \n
    @param code: \e string \n
        An indicator for which code to run; options = 'mcnp', 'mcnp6',
        'mcnp6.mpi', 'advantg'. [Default = 'mcnp6'] \n
    """
    module_logger.debug("In Run Transport, the lst input = {}, nps = {}, and \
                        code is = {}".format(lst, nps, code))

    # Start Clock
    start_time = time.time()
    runFiles = []

    # Ensure log directories are ready and clean old files
    path = os.path.abspath(os.path.join(os.path.abspath(os.getcwd()),
                                        os.pardir))
    if not os.path.isdir("{}/logs/".format(path)):
        os.mkdir("{}/logs/".format(path))
    elif os.listdir("{}/logs/".format(path)):
        sub.Popen("rm {}/logs/*".format(path), cwd=path, stdout=sub.PIPE,
                  shell=True)

    # Ensure output directories are ready and clean old files
    for i in lst:
        if not os.path.isdir("{}/Results/Population/{}/tmp/".format(path, i)):
            os.mkdir("{}/Results/Population/{}/tmp/".format(path, i))
        else:
            sub.Popen("rm -rf {}/Results/Population/{}/tmp/*".format(path, i),
                      cwd=path, stdout=sub.PIPE, shell=True)

    # Define number of tasks to assign to each run
    cores = mp.cpu_count()
    module_logger.info("The number of cores avaliable is = {}".format(cores))

    # Run MCNP
    tasks = []
    if code in ["mcnp", "mcnp6", "mcnp6.mpi"]:
        if nps == []:
            for item in lst:
                tasks.append(cores)
        else:
            for n in nps:
                if n <= 1E6:
                    tasks.append(cores)
                elif n <= 1E7:
                    tasks.append(cores*4)
                elif n <= 1E8:
                    tasks.append(cores*12)
                elif n > 1E8:
                    tasks.append(cores*14)
                else:
                    module_logger.error('\nThe nps condition wasn't covered. '
                                        'NPS = {}'.format(n))

        module_logger.debug('Number of Cores = {}'.format(cores))
        module_logger.debug('Number of Tasks = {}\n'.format(tasks))

        # Determine unique numbers of tasks to set number of batch files
        task_set = sorted(set(tasks), reverse=True)
        module_logger.info('Unique Task Identifiers = {}\n'.format(task_set))
        module_logger.debug('lst = {}\n'.format(lst))
        module_logger.debug('nps = {}\n'.format(nps))
        module_logger.debug('tasks = {}\n'.format(tasks))

        for t in task_set:

            # Build sub list
            subLst = []
            for i in range(0, len(tasks)):
                if tasks[i] == t:
                    subLst.append(lst[i])
            module_logger.info('For {} tasks, the sub list is = {}'.format(t,
                               subLst))

            # Build batch
            if (t < cores and len(subLst)%2 == 0) or t >= cores:
                runFiles.append(build_batch(subLst, t, code, *batchArgs))
            else:
                runFiles.append(build_batch(subLst[0:-1], t, code, *batchArgs))
                runFiles.append(build_batch([subLst[-1]], t, code, *batchArgs,
                                             suf="a"))

    # Run ADVANTG
    elif code == "advantg":
        # Build batch
        fname = build_batch(lst, cores, code, *batchArgs)

        # Copy files into correct run directory
        for i in lst:
            if os.path.isfile("{}/{}".format(os.path.abspath(os.getcwd()),
                              fname)):
                sub.Popen("cp {} {}".format(os.path.abspath(os.getcwd()) +
                          "/"+fname,
                          path+"/Results/Population/"+str(i)+"/tmp/"+fname),
                          cwd=path, stdout=sub.PIPE, shell=True)
            else:
                module_logger.info('{}/{} doesnt exist.'.format(
                                   os.path.abspath(os.getcwd()), fname))

            runFiles.append(fname)

            if os.path.isfile("{}/Results/Population/{}/runCADIS.adv".format(
                              path, str(i))):
                sub.Popen("cp {} {}".format(path+"/Results/Population/"
                          +str(i)+"/runCADIS.adv", path+"/Results/Population/"
                          +str(i)+"/tmp/runCADIS.adv"), cwd=path,
                          stdout=sub.PIPE, shell=True)
            else:
                module_logger.info("{}/Results/Population/{}/runCADIS.adv "
                                   "doesn't exist. ".format(path, str(i)))
    else:
        module_logger.warning('Unknown code ({}) specified. Please try again.'
                              '\n'.format(code))

    job_id_list = []

    # Execute batch
    # runFiles should contains the second ID part of mcnp jobs
    module_logger.info("The runFiles are: {}".format(runFiles))
    for i in range(0, len(runFiles)):
        cmd = "sbatch {}".format(runFiles[i])
        if code == 'advantg':
            rundir = path+"/Results/Population/"+str(lst[i])+"/tmp/"
            jobOut = sub.Popen(cmd, cwd=rundir, stdin=sub.PIPE,
                               stdout=sub.PIPE, stderr=sub.PIPE,
                               shell=True).communicate()[0].strip().split()
            module_logger.info("ADVANTG job submission communication: {}"
                               "".format(jobOut))
        elif code in ["mcnp", "mcnp6", "mcnp6.mpi"]:
            jobOut = sub.Popen(cmd, cwd=os.path.abspath(os.getcwd()),
                               stdin=sub.PIPE, stdout=sub.PIPE,
                               stderr=sub.PIPE,
                               shell=True).communicate()[0].strip().split()
            module_logger.info("MCNP job submission communication: {}"
                               "".format(jobOut))
        if jobOut:
            job_id_list.append(jobOut[3])

    # Monitor for completion
    time.sleep(10)
    module_logger.info("job ids: {}".format(job_id_list))
    def monitor():
        output = []
        for jobid in job_id_list:
            job = sub.Popen("squeue | grep " + jobid, cwd=path,
                            stdout=sub.PIPE,
                            shell=True).communicate()[0].strip().split()
            if job:
                output.append(job[0])
        return output

    output = monitor()
    module_logger.info("monitor output={}\n".format(output))
    while output:
        output = monitor()
        module_logger.debug("\n\n\nLen(full_out)={}, Line 1 of Squeue output "
                            "= {}".format(len(output), output))
        time.sleep(1)

    # Copy ADVANTG generated inputs to correct directory
    if code == 'advantg':
        time.sleep(10)
        for i in lst:
            path = os.path.abspath(os.path.join(os.path.abspath(os.getcwd()),
                                   os.pardir)) + \
                                   "/Results/Population/"+str(i)+"/"
            sub.Popen("cp {} {}".format(path+"tmp/output/wwinp", path+"wwinp"),
                      cwd=path, stdout=sub.PIPE, shell=True)
            sub.Popen("cp {} {}".format(path+"tmp/output/inp_edits.txt",
                      path+"inp_edits.txt"), cwd=path, stdout=sub.PIPE,
                      shell=True)
            sub.Popen("rm -rf {}tmp/*".format(path), cwd=path,
                      stderr=sub.STDOUT, stdout=sub.PIPE, shell=True)

    module_logger.info('Total transport time was {} sec'.format(time.time() -
                                                                start_time))

#-----------------------------------------------------------------------------#
def build_batch(lst, tasks, code, qos, account, partition, timeout, scheduler,
                suf=""):
    """!
    Build a Slurm Batch script using the Jobs Array feature to run transport
    calculations.

    @param lst: \e list \n
        List of parent identifier numbers to be ran. \n
    @param tasks: \e list \n
       Number of tasks to run per code thread instance. \n
    @param code: \e string \n
        An indicator for which code to run; options = 'mcnp6', 'mcnp6.mpi',
        'advantg'. [Default = 'mcnp6'] \n
    @param qos: \e string \n
        Quality of service to be used for job submission. \n
    @param account: \e string \n
        Account to be used for job submission. \n
    @param partition: \e string \n
        Partition to be used for job submission. \n
    @param timeout: \e string \n
        Timeout time to be used for job submission. \n
    @param scheduler: \e string \n
        The job scheduler used; options are 'slurm'. \n
    @param suf: \e string \n
        Optional string identifier suffix to be added at end of file. \n

    @return \e string Filename for the batchfile created.
    """
    #Sort the identifiers
    lst = sorted(lst)

    # Determine whether to use weight windows
    ww = "wwinp=../Results/Population/$SLURM_ARRAY_TASK_ID/wwinp"
    for i in lst:
        if not os.path.isfile("../Results/Population/{}/wwinp".format(str(i))):
            ww = ""
            module_logger.debug('Running without weight windows.')
            break

    # Determine whether to specify tasks
    t_str = ''
    if code in ["mcnp", "mcnp6", "mcnp6.mpi"]:
        t_str = 'tasks {}'.format(tasks)

    # Set filename
    path = "{}/".format(os.path.abspath(os.getcwd()))
    fname = "run{}_{}{}.sh".format(code, tasks, suf)

    if scheduler.lower() == 'slurm':
        try:
            with open(path+fname, "w") as f:
                f.write("#!/bin/sh\n\n")
                f.write("#SBATCH --time=" + timeout +"\n")
                f.write("# Job name:\n")

                if code in ["mcnp", "mcnp6", "mcnp6.mpi"]:
                    f.write("#SBATCH --job-name=mc{}\n".format(tasks))
                elif code == "advantg":
                    f.write("#SBATCH --job-name=adv{}\n".format(tasks))

                f.write("# Partition:\n")
                f.write("#SBATCH --partition=" + partition + "\n")
                f.write("# QoS:\n")

                f.write("#SBATCH --qos=" + str(qos) + "\n")
                f.write("# Account:\n")
                f.write("#SBATCH --account=" + account + "\n")

                f.write("# Processors:\n")
                f.write("#SBATCH --ntasks={}\n".format(tasks))

                if code in ["mcnp", "mcnp6", "mcnp6.mpi"]:
                    if code == "mcnp6.mpi":
                        code = "mpirun mcnp6.mpi"
                    elif code == "mcnp":
                        code = "mcnp6"

                    f.write("#SBATCH --output=../logs/arrayJob_%A_%a.out\n")
                    f.write("#SBATCH --error=../logs/arrayJob_%A_%a.err\n")
                    f.write("# Array:\n")
                    f.write("#SBATCH --array={}\n\n".format(
                                                ",".join(str(l) for l in lst)))
                    f.write("module load openmpi\n")
                    wd = "../Results/Population/$SLURM_ARRAY_TASK_ID/"
                    f.write("{0} i={1}ETA.inp o={1}tmp/ETA.out \
                            run={1}tmp/runtpe {2} {3}\n".format(code, wd, ww,
                                                                t_str))

                elif code == "advantg":
                    f.write("#SBATCH --output=slurm_%j.out\n")
                    f.write("#SBATCH --error=slurm_%j.err\n")
                    f.write("{} runCADIS.adv\n".format(code))

            # Close the file
            f.close()

        except IOError as e:
            module_logger.error("I/O error({0}): {1}".format(e.errno,
                                                             e.strerror))

    # Test that the file closed
    assert f.closed, "File did not close properly."

    module_logger.debug('Built {}'.format(path+fname))

    return fname

#-----------------------------------------------------------------------------#
def to_Norm(spectrum):
    """!
    Normalizes a MCNP tallied flux.

    @param spectrum: \e array \n
        An MCNP tally input spectrum. \n

    @return \e array The output normalized flux spectrum. \n
    """

    flux = np.zeros(len(spectrum[:, 0]))
    return flux/np.sum(flux)

#-----------------------------------------------------------------------------#
def to_NormDiff(spectrum):
    """!
    Converts a MCNP tallied flux to a Normalized Differential flux.

    @param spectrum: \e array \n
        An MCNP tally input spectrum. \n

    @return \e array The output normalized differential flux spectrum. \n
    """

    # Initialize variables
    diff = np.zeros(len(spectrum[:, 0]))

    # Calculate the differential flux
    diff[0] = (spectrum[0, 1])/(spectrum[0, 0])
    for i in range(1, len(spectrum[:, 0])):
        diff[i] = (spectrum[i, 1])/(spectrum[i, 0]-spectrum[i-1, 0])

    # Calculate the normalized differential flux
    return diff/np.sum(diff)
