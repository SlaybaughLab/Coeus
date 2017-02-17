
## Module : Utilities.py
#
## Contains : Utility functions for the Coeus program.
#
## Author : James Bevins
#
## Last Modified: 17Oct16

import logging
module_logger = logging.getLogger('Coeus.Utilities')

from threading import Thread

import multiprocessing as mp
import subprocess as sub
import numpy as np

import os
import glob
import time
import sys
import threading
import shutil
import bisect



## Creates a switch class object to switch between cases
class Switch(object):
    
    ## The constructor.
    # @param value selector value
    def __init__(self, value):
        ## string case selector value
        self.value = value
        ## boolean based on match
        self.fall = False

    ## Return the match method once, then stop
    def __iter__(self):
        yield self.match
        raise StopIteration


    ## PrintIndicate whether or not to enter a case suite
    # @param args list of arguments to match with
    def match(self, *args):
        if self.fall or not args:
            return True
        elif self.value in args: 
            self.fall = True
            return True
        else:
            return False   
        
## Creates a Thread class object to run command line programs in parallel
class Cmd_Thread(Thread):

    ## The constructor
    # @param cwdir Current working directory path
    # @param cmd The command line input to be executed
    def __init__(self,cwdir,cmd):
        Thread.__init__(self)
        ## Current working directory path
        self.cwdir = cwdir
        ## The command line input to be executed
        self.cmd=cmd
        
    def __repr__(self):
        return "Thread instance({0}, {1})".format(self.cwdir, self.cmd)
    
    def __str__(self):
        header = ["\nThread Instance:"]
        header += ["The current working directory is = {}".format(self.cwdir)]
        header += ["The cmd line input is = {}".format(self.cmd)]
        header ="\n".join(header)+"\n"
        return header

    ## Run Thread in local working directory
    def run(self):             
        t=sub.Popen(self.cmd,cwd=self.cwdir,shell=True)
        t.communicate()
        if t.returncode !=0:
            module_logger.error("The thread {} did not execute properly.".format(self.getName()))
            sys.exit()

## Creates a Thread class object to run functions without returns in parallel.
class FuncThread(Thread):

    ## The constructor
    # @param target The function to be executed
    # @param args The functions arguments
    def __init__(self, target, *args):
        ## The function to be executed
        self._target = target
        ## The functions arguments
        self._args = args
        Thread.__init__(self)
 
    def run(self):
        self._target(*self._args)      


## Creates a Thread class object to run functions containing returns in parallel.
class FuncThreadWithReturn(Thread):

    # @param *args The functions arguments
    def __init__(self, *args, **kwargs):
        super(FuncThreadWithReturn, self).__init__(*args, **kwargs)

        self._return = None

    def run(self):
        if self._Thread__target is not None:
            self._return = self._Thread__target(*self._Thread__args, **self._Thread__kwargs)

    def join(self, *args, **kwargs):
        super(FuncThreadWithReturn, self).join(*args, **kwargs)

        return self._return    
    
## Runs a multi-threaded transport calculation. Doesn't work for clusters.
# @param lst list of parent identifier number to be ran
# @param tasks Number of tasks to run per code thread instance. If left blank, calculation will be performed to assign all
# availiable cpus evenly
# @param code  [Default = 'mcnp6'] An indicator for which code to run  (options = 'mcnp6', 'mcnp6.mpi', 'advantg')
def Run_Transport_Threads(lst,tasks=0,code='mcnp6'):
    #Start Clock
    start_time=time.time()
    # Initialize thread list
#   thread_lst=[]
    processes=[]
    
    # Define number of threads to run at once
    cores=mp.cpu_count()
    tasks=cores/len(lst)
    if tasks==0:
        tasks=1
    module_logger.debug("\n\nNumber of Cores = {}".format(cores))
    module_logger.debug("Number of Tasks = {}\n\n".format(tasks))
    
    # Create thread for each current solution
    for i in lst:

        # Define path to current parent run directory
        path=os.path.abspath(os.path.join(os.path.abspath(os.getcwd()),os.pardir))+"/Results/Population/"+str(i)+"/tmp"
          
        # Create subdirectory for each parent if it doesn't exist
        if os.path.isdir(path)==False:
            os.mkdir(str(path))
            
        # Set the cmd line input
        if code.strip().lower()=='mcnp6':
            # Clean up run directory
            if os.path.isfile("{}/ETA.out".format(path)):
                os.remove("{}/ETA.out".format(path))
            if os.path.isfile("{}/runtpe".format(path)):
                os.remove("{}/runtpe".format(path))
            if os.path.isfile("{}/../wwinp".format(path)):
                cmd="mcnp6 i=../ETA.inp o=ETA.out wwinp=../wwinp tasks {}".format(tasks)
            else:
                cmd="mcnp6 i=../ETA.inp o=ETA.out tasks {}".format(tasks)
        elif code.strip().lower()=='mcnp6.mpi':
            # Clean up run directory
            if os.path.isfile("{}/ETA.out".format(path)):
                os.remove("{}/ETA.out".format(path))
            if os.path.isfile("{}/runtpe".format(path)):
                os.remove("{}/runtpe".format(path))
            if os.path.isfile("{}/../wwinp".format(path)):
                cmd="mpirun mcnp6.mpi i=../ETA.inp o=ETA.out wwinp=../wwinp"
            else:
                cmd="mpirun mcnp6.mpi i=../ETA.inp o=ETA.out"
        elif code.strip().lower()=='advantg':
            # Clean up run directory
            shutil.rmtree(path)
            os.mkdir(str(path))
            cmd="advantg ../runCADIS.adv"
        else:
            module_logger.warning("Unknown radiation transport code specified: {}.  Only ADVANTG, MCNP6, and MCNP6.mpi are valid options.".format(code))
        
        processes.append(sub.Popen(cmd,cwd=path,shell=True))
        
    # Create processes    
    for p in processes:
        p.wait()
    
    # Copy ADVANTG generated inputs to correct directory
    if code=='advantg':
        for i in lst:
            path=''
            path=os.path.abspath(os.path.join(os.path.abspath(os.getcwd()),os.pardir))+"/Resuts/Population/"+str(i)+"/"
            shutil.copyfile(path+'tmp/output/wwinp',path+'wwinp')
            shutil.copyfile(path+'tmp/output/inp_edits.txt',path+'inp_edits.txt')
            shutil.rmtree(path+'tmp/output')
            shutil.rmtree(path+'tmp/model')
            shutil.rmtree(path+'tmp/adj_solution')
    
    module_logger.info('Total transport time was {} sec'.format(time.time() - start_time))
    
##  A callable function to execute a command line program.
#   @param cwdir Current working directory path
#   @param cmd The command line input to be executed
def Run_CmdLine(cmd,cwdir):
    try:    
        #proc = subprocess.Popen(['/home/pyne-user/MCNP/MCNP_CODE/bin/mcnp6', 'i=../ETA.inp o=ETA.out'],cwd=cwdir,stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)#,shell=True)
        proc = subprocess.Popen(cmd,cwd=cwdir,stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=True)
        output = proc.communicate()
    except:
        module_logger.error('Simulation not even started')
    
    print output     
        
  
##Runs a multi-threaded transport calculation. Doesn't work for clusters.
# @param lst list of parent identifier numbers to be ran
# @param tasks Number of tasks to run per code thread instance. If left blank, calculation will be performed to assign all
#    availiable cpus evenly
# @param code [Default = 'mcnp6'] An indicator for which code to run  (options = 'mcnp6', 'mcnp6.mpi', 'advantg')
def Run_Transport_PP(lst,tasks=0,code='mcnp6'):
    #Start Clock
    start_time=time.time()

    # Initialize job queue
    jobs=[]
    ppservers=("*",)
    job_server = pp.Server(ppservers=ppservers) 
    
    # Determine computational resources available
    module_logger.info('The number of servers is: {}'.format(ppservers))
    module_logger.info('The number of cpus is: {}'.format(job_server.get_ncpus()))
    module_logger.info('The number of nodes is: {}'.format(job_server.get_active_nodes()))
    
    # Define number of tasks to assign to each run
    cores=mp.cpu_count()
    tasks=cores/len(lst)
    if tasks==0:
        tasks=1
    module_logger.debug("\n\nNumber of Cores = {}".format(cores))
    module_logger.debug("Number of Tasks = {}\n\n".format(tasks))

    # Create thread for each current solution
    for i in lst:

        # Define path to current parent run directory
        path=''
        path=os.path.abspath(os.path.join(os.path.abspath(os.getcwd()), os.pardir))+"/Results/Population/"+str(i)+"/tmp"
          
        # Create subdirectory for each parent if it doesn't exist
        if os.path.isdir(path)==False:
            os.mkdir(str(path))
            
        # Set the cmd line input
        if code.strip().lower()=='mcnp6':
            # Clean up run directory
            if os.path.isfile("{}/ETA.out".format(path)):
                os.remove("{}/ETA.out".format(path))
            if os.path.isfile("{}/runtpe".format(path)):
                os.remove("{}/runtpe".format(path))
            if os.path.isfile("{}/../wwinp".format(path)):
                cmd="mcnp6 i=../ETA.inp o=ETA.out wwinp=../wwinp tasks {}".format(tasks)
            else:
                cmd="mcnp6 i=../ETA.inp o=ETA.out tasks {}".format(tasks)
        elif code.strip().lower()=='mcnp6.mpi':
            # Clean up run directory
            if os.path.isfile("{}/ETA.out".format(path)):
                os.remove("{}/ETA.out".format(path))
            if os.path.isfile("{}/runtpe".format(path)):
                os.remove("{}/runtpe".format(path))
            if os.path.isfile("{}/../wwinp".format(path)):
                cmd="mpirun mcnp6.mpi i=../ETA.inp o=ETA.out wwinp=../wwinp"
            else:
                cmd="mpirun mcnp6.mpi i=../ETA.inp o=ETA.out"
        elif code.strip().lower()=='advantg':
            # Clean up run directory
            shutil.rmtree(path)
            os.mkdir(str(path))
            cmd="advantg ../runCADIS.adv"
        else:
            module_logger.warning("Unknown radiation transport code specified: {}.  Only ADVANTG, MCNP6, and MCNP6.mpi are valid options.".format(code))
        
        jobs.append(job_server.submit(Run_CmdLine,(cmd,path,),modules=("subprocess","time",)))
        
    job_server.wait()
    job_server.destroy()
    
    # Copy ADVANTG generated inputs to correct directory
    if code=='advantg':
        for i in lst:
            path=''
            path=os.path.abspath(os.path.join(os.path.abspath(os.getcwd()), os.pardir))+"/Results/Population/"+str(i)+"/"
            shutil.copyfile(path+'tmp/output/wwinp',path+'wwinp')
            shutil.copyfile(path+'tmp/output/inp_edits.txt',path+'inp_edits.txt')
            shutil.rmtree(path+'tmp/output')
            shutil.rmtree(path+'tmp/model')
            shutil.rmtree(path+'tmp/adj_solution')
    
    module_logger.info("Job server stats: {}".format(job_server.print_stats()))
    module_logger.info('Total transport time was {} sec'.format(time.time() - start_time))
    
## Build a Slurm Batch script using the Jobs Array feature to run transport calculations. 
# @param  lst list of parent identifier numbers to be ran
# @param nps list of number of particles to run per code thread instance. If left blank, calculation will be performed to assign all
#    availiable cpus evenly
# @param code [Default = 'mcnp6'] An indicator for which code to run  (options = 'mcnp6', 'mcnp6.mpi', 'advantg')
def Run_Transport(lst, qos, account, partition, timeout, nps=[],code='mcnp6'):
    module_logger.debug("In Run Transport, the lst input = {}, nps = {}, and code is = {}".format(lst,nps,code))
    
    # Start Clock
    start_time=time.time()
    run_files=[]    
    
    # Ensure output directories are ready and clean old files
    path=os.path.abspath(os.path.join(os.path.abspath(os.getcwd()), os.pardir))
    if os.path.isdir("{}/logs/".format(path))==False:
        os.mkdir("{}/logs/".format(path))
    elif os.listdir("{}/logs/".format(path)):
        sub.Popen("rm {}/logs/*".format(path),cwd=path,stdout=sub.PIPE,shell=True)
    for i in lst:
        if os.path.isdir("{}/Results/Population/{}/tmp/".format(path,i))==False:
            os.mkdir("{}/Results/Population/{}/tmp/".format(path,i))
        else:
            sub.Popen("rm -r {}/Results/Population/{}/tmp/*".format(path,i),cwd=path,stdout=sub.PIPE,shell=True)
    
    if code=="mcnp6" or code=="mcnp6.mpi":
        # Define number of tasks to assign to each run
        cores=mp.cpu_count()
        if nps==[]:
            tasks=cores/len(lst)
            if tasks==0:
                tasks=1
        else:
            tasks=[]
            for n in nps:
                if n<=1E6:
                    tasks.append(cores)
                elif n<=1E7:
                    tasks.append(cores*4)
                elif n<=1E8:
                    tasks.append(cores*12)
                elif n>1E8:
                    tasks.append(cores*14)
                else:
                    module_logger.error("\nThe nps condition wasn't covered. NPS = {}".format(n))

        module_logger.debug("Number of Cores = {}".format(cores))
        module_logger.debug("Number of Tasks = {}\n".format(tasks))

        # Determine unique numbers of tasks to set number of batch files
        task_set=sorted(set(tasks),reverse=True)
        module_logger.info("Unique Task Identifiers = {}\n".format(task_set))
        module_logger.debug("lst = {}\n".format(lst))
        module_logger.debug("nps = {}\n".format(nps))
        module_logger.debug("tasks = {}\n".format(tasks))

        for t in task_set:

            # Build sub list
            sub_lst=[]
            for i in range(0,len(tasks)):
                if tasks[i] == t:
                    sub_lst.append(lst[i])
            module_logger.info("For {} tasks, the sub list is = {}".format(t,sub_lst))
            
            # Build batch
            if (t < 20 and len(sub_lst)%2==0) or t >= 20:
                run_files.append(Build_Batch(sub_lst,t,code, qos, account, partition, timeout))
            else:
                run_files.append(Build_Batch(sub_lst[0:-1],t,code, qos, account, partition, timeout))
                run_files.append(Build_Batch([sub_lst[-1]],t,code, qos, account, partition, timeout, suf="a"))
                
    elif code=="advantg":
        # Build batch
        fname=Build_Batch(lst,20,code, qos, account, partition, timeout)
        
        # Copy files into correct run directory
        for i in lst:
            if os.path.isfile("{}/{}".format(os.path.abspath(os.getcwd()),fname)):
                sub.Popen("cp {} {}".format(os.path.abspath(os.getcwd())+"/"+fname,path+"/Results/Population/"+str(i)+"/tmp/"+fname), cwd=path, stdout=sub.PIPE, shell=True)
            else: 
                module_logger.info("{}/{} doesnt exist. ".format(os.path.abspath(os.getcwd()),fname))
                
            run_files.append(fname)
            
            if os.path.isfile("{}/Results/Population/{}/runCADIS.adv".format(path,str(i))):
                sub.Popen("cp {} {}".format(path+"/Results/Population/"+str(i)+"/runCADIS.adv",path+"/Results/Population/"+str(i)+"/tmp/runCADIS.adv"), cwd=path, stdout=sub.PIPE, shell=True)
            else:
                module_logger.info("{}/Results/Population/{}/runCADIS.adv doesn't exist. ".format(path,str(i)))
    else:
        module_logger.warning("Unknown code ({}) specified. Please try again. \n".format(code))
        
    # Execute batch
    main_jobid=sub.Popen("squeue | grep jbevins",cwd=path,stdout=sub.PIPE,shell=True).communicate()[0].strip().split()[0]
    module_logger.debug("main_jobid={}\n".format(main_jobid))
    for i in range(0,len(run_files)): # run_files should contains the second ID part of mcnp jobs
        cmd="sbatch {}".format(run_files[i])
        if code == 'advantg':
            rundir=path+"/Results/Population/"+str(lst[i])+"/tmp/"
            jobOut=sub.Popen(cmd,cwd=rundir,stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE,shell=True).communicate()
            module_logger.info("ADVANTG job submission communication: {}".format(jobOut))
        else:

            sub.Popen(cmd,cwd=os.path.abspath(os.getcwd()),stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE,shell=True)
 
    # Monitor for completion
    time.sleep(15)
    output=sub.Popen("squeue | grep jbevins",cwd=path,stdout=sub.PIPE,shell=True).communicate()[0]
    while output.strip().split()[0] != main_jobid or len(output.split()) > 8:
        output=sub.Popen("squeue | grep jbevins",cwd=path,stdout=sub.PIPE,shell=True).communicate()[0]
        if output.strip().split()[0] == main_jobid and len(output.split()) <= 8:
            time.sleep(1)
            output=sub.Popen("squeue | grep jbevins",cwd=path,stdout=sub.PIPE,shell=True).communicate()[0]
        module_logger.debug("\n\n\nLen(full_out)={}, Line 1 of Squeue output = {}".format(len(output),output))
        time.sleep(1)
    
    # Copy ADVANTG generated inputs to correct directory
    if code=='advantg':
        for i in lst:
            time.sleep(5)
            path=os.path.abspath(os.path.join(os.path.abspath(os.getcwd()), os.pardir))+"/Results/Population/"+str(i)+"/"
            sub.Popen("cp {} {}".format(path+"tmp/output/wwinp",path+"wwinp"),cwd=path,stdout=sub.PIPE,shell=True)
            sub.Popen("cp {} {}".format(path+"tmp/output/inp_edits.txt",path+"inp_edits.txt"), cwd=path,stdout=sub.PIPE,shell=True)
            sub.Popen("rm -rf {}tmp/*".format(path,i),cwd=path,stderr=sub.STDOUT,stdout=sub.PIPE,shell=True)
            
    module_logger.info('Total transport time was {} sec'.format(time.time() - start_time))
        

## Build a Slurm Batch script using the Jobs Array feature to run transport calculations. 
# @param  lst list of parent identifier numbers to be ran
# @param tasks Number of tasks to run per code thread instance
# @param code [Default = 'mcnp6'] An indicator for which code to run  (options = 'mcnp6', 'mcnp6.mpi', 'advantg')
# @param suf Optional string identifier suffix to be added at end of file
# @return Filename for the batchfile created
def Build_Batch(lst,tasks,code, qos, account, partition, timeout, suf=""):
    #Sort the identifiers
    lst=sorted(lst)       
    
    # Determine whether to use weight windows
    ww="wwinp=../Results/Population/$SLURM_ARRAY_TASK_ID/wwinp"
    for i in lst:
        if os.path.isfile("{}/wwinp".format(os.path.abspath(os.path.join(os.path.abspath(os.getcwd()), os.pardir))+"/Results/Population/"+str(i)))==False:
            ww=""
            module_logger.debug('Running without weight windows.')
            break
    
    # Determine whether to specify tasks
    t_str=''
    if code=='mcnp6':
        t_str='tasks {}'.format(tasks)
        
    # Set filename
    path="{}/".format(os.path.abspath(os.getcwd()))
    fname="run{}_{}{}.sh".format(code,tasks,suf)
            
    try:
        with open(path+fname, "w") as f:
            f.write("#!/bin/sh\n\n")
            f.write("#SBATCH --time=" + timeout +"\n")
            f.write("# Job name:\n")
            
            if code == "mcnp6.mpi" or code=="mcnp6":
                f.write("#SBATCH --job-name=mc{}\n".format(tasks))
            elif code == "advantg":
                f.write("#SBATCH --job-name=adv{}\n".format(tasks))
                
            f.write("# Partition:\n")
            f.write("#SBATCH --partition=" + partition + "\n") 
            f.write("# QoS:\n")
            
            f.write("#SBATCH --qos=" + qos + "\n")
            f.write("# Account:\n")
            f.write("#SBATCH --account=" + account + "\n")

            f.write("# Processors:\n")
            f.write("#SBATCH --ntasks={}\n".format(tasks))

            if code == "mcnp6.mpi" or code=="mcnp6":
                if code == "mcnp6.mpi":
                    code="mpirun mcnp6.mpi"
                    
                f.write("#SBATCH --output=../logs/arrayJob_%A_%a.out\n")
                f.write("#SBATCH --error=../logs/arrayJob_%A_%a.err\n")
                f.write("# Array:\n")
                f.write("#SBATCH --array={}\n\n".format(",".join(str(l) for l in lst))) 
                f.write("module load openmpi\n")
                f.write("{} i=../Results/Population/$SLURM_ARRAY_TASK_ID/ETA.inp o=../Results/Population/$SLURM_ARRAY_TASK_ID/tmp/ETA.out run=../Results/Population/$SLURM_ARRAY_TASK_ID/tmp/runtpe {} {}\n".format(code,ww,t_str))
                    
            elif code == "advantg":
                f.write("#SBATCH --output=slurm_%j.out\n")
                f.write("#SBATCH --error=slurm_%j.err\n")
                f.write("{} runCADIS.adv\n".format(code))
                
        # Close the file
        f.close()

    except IOError as e:
        module_logger.error("I/O error({0}): {1}".format(e.errno, e.strerror)) 

    # Test that the file closed
    assert f.closed==True, "File did not close properly."
    
    module_logger.debug('Built {}'.format(path+fname))
        
    return fname

## Normalizes a MCNP tallied flux 
# @param spectrum array Teh input flux spectrum
# @return result array the output normalized differential flux spectrum   
def to_Norm(spectrum):
    
    flux=np.zeros(len(spectrum[:,0]))
    result=flux/np.sum(flux)
        
    return result

## Converts a MCNP tallied flux to a Normalized Differential flux
# @param spectrum The input flux spectrum
# @return The output normalized differential flux spectrum
def to_NormDiff(spectrum):
    # Initialize variables
    diff=np.zeros(len(spectrum[:,0]))
    intdiff=np.zeros(len(spectrum[:,0]))
    normdiff=np.zeros(len(spectrum[:,0]))
    result=np.zeros(len(spectrum[:,0]))
    
    # Calculate the differential flux
    diff[0]=(spectrum[0,1])/(spectrum[0,0])
    for i in range(1,len(spectrum[:,0])):
        diff[i]=(spectrum[i,1])/(spectrum[i,0]-spectrum[i-1,0])
    
    # Calculate the normalized differential flux 
    result=diff/np.sum(diff)
        
    return result
    
## Calculate the U-optimality 
# @param c the candidate design
# @param d the objective design
# @return The u-optimality design based fitness
def Uopt(c,d):
    assert len(c)==len(d), "The length of the candidate and objective design must be equal in Uopt."  
   
    return np.sum(abs(d-c))

    
## Calculate the U-optimality 
# @param c the candidate design
# @param d the objective design
# @return The least-squares design based fitness
def LeastSquares(c,d):
    
    assert len(c)==len(d), "The length of the candidate and objective design must be equal in LeastSquares."  
   
    return np.sum((d-c)**2)

## Calculates the relative least squares.  Assumes a normalized input candidate and objective spectrum to simplify  
#    calculation (i.e. sum of bins should equal 1). Not valid for unnormalized spectra.  
# @param c the candidate design
# @param o the objective design
# @return The least-squares design based fitness
def RelativeLeastSquares(c,o):  
    assert len(c)==len(o), "The length of the candidate and objective design must be equal in RelativeLeastSquares."  
    rls=(o-c)**2/o
    
    # For bins with no tally results, project the fitness 
    loc=len(c)-6
    while c[loc]!=0.0:
        loc-=1
        if loc == -1:
            break
    rls[0:loc+1]=np.array([np.average(rls[loc+1:loc+4])]*(loc+1))
    return np.sum(rls) 

    
## an event object representing a snapshot in the optimization process
class Event:

    ## Creates an event object representing a snapshot in the optimization process
    # @return None
    def __init__(self,generation,evaluations,fitness,nps,ident):
        assert generation >= 0, "The number of generations cannot be negative."
        assert evaluations >= 0, "The number of evaluations cannot be negative."
#        assert fitness >= 0, "The fitness cannot be negative."
        assert isinstance(generation, int)==True, 'Generation must be of type int.'
        assert isinstance(evaluations, int)==True, 'Evaluations must be of type int.'
        assert isinstance(fitness, float)==True, 'Fitness must be of type float.'
        ## The generation the design was arrived at
        self.g=generation
        ## The number of fitness evaluations done to obtain this design
        self.e=evaluations
        ## The assessed design fitness
        self.f=fitness 
        ## The number of particles run for that event
        self.n=nps
        ## The identifty of the current top solution
        self.i=ident
        
    def __repr__(self):
        return "Event instance({0}, {1}, {2}, {3}, {4})".format(self.g, self.e, self.f, self.n, self.i)
    
    def __str__(self):
        header = ["\nEvent Instance:"]
        header += ["Generation # = {}".format(self.g)]
        header += ["Number of Evaluations = {}".format(self.e)]
        header += ["Fitness = {}".format(self.f)]
        header += ["NPS = {}".format(self.n)]
        header += ["Identity = {}".format(self.i)]
        header ="\n".join(header)+"\n"
        return header

## Defines a class of weights to be used to select number of instances in array randomly with
#    linear weighting. 
class WeightedRandomGenerator(object):

    # @param self Current instance of the class
    # @param weights The array of weights (Higher = more likely to be selected)
    # @return  The randomly selected index of the weights array
    def __init__(self, weights):
        self.totals = []
        running_total = 0

        for w in weights:
            running_total += w
            self.totals.append(running_total)

    # @return  The randomly selected index of the weights array
    def next(self):
        rnd = np.random.rand() * self.totals[-1]
        return bisect.bisect_right(self.totals, rnd)

    def __call__(self):
        return self.next()
        
## Stores and prints effectiveness stats for each metaheuristic search method.
class Meta_Stats():
    
    ## Initializer
    # @param mat_levy tuple contains the changes and total number of function evaluations for the Mat_Levy_Flights function
    # @param cell_levy tuple contains the changes and total number of function evaluations for the Cell_Levy_Flights function
    # @param elite_cross tuple contains the changes and total number of function evaluations for the Mutate_Mats function
    # @param part_inv tuple contains the changes and total number of function evaluations for the Partial_Inversion function
    # @param mutate tuple contains the changes and total number of function evaluations for the Mutate function
    # @param two_opt tuple contains the changes and total number of function evaluations for the 2-opt function
    # @param crossover tuple contains the changes and total number of function evaluations for the Crossover function
    # @param three_op tuple contains the changes and total number of function evaluations for the Three_opt function
    # @param discard tuple contains the changes and total number of function evaluations for the Discard function
    def __init__(self, mat_levy=(0,0), cell_levy=(0,0), elite_cross=(0,0), part_inv=(0,0), mutate=(0,0), two_opt=(0,0), crossover=(0,0), three_op=(0,0), discard=(0,0), fname=os.path.abspath(os.path.join(os.getcwd(), os.pardir))+"/Results/meta_stats.txt"):
        ## dictionary string name of each algorithms
        self.algorithms = {"mat_levy": mat_levy, "cell_levy": cell_levy, 
                    "elite_cross": elite_cross, "part_inv":part_inv, "mutate": mutate, 
                    "two_opt":two_opt, "crossover":crossover,"three_op": three_op, 
                    "discard": discard}
        ## str Name and path of the file to store the timeline for post processing
        self.fname=fname

        if os.path.isfile(self.fname)==True:
            os.remove(self.fname)
        
    def __repr__(self):
        return "{:8d},{};  {:14d},{};  {:11d},{};  {:12d},{};  {:12d},{};  {:12d},{};  {:11d},{};  {:9d},{};  {:6d},{};\n".format(self.algorithms["mat_levy"][0], self.algorithms["mat_levy"][1], self.algorithms["cell_levy"][0], self.algorithms["cell_levy"][1], self.algorithms["elite_cross"][0], self.algorithms["elite_cross"][1], self.algorithms["part_inv"][0], self.algorithms["part_inv"][1], self.algorithms["mutate"][0], self.algorithms["mutate"][1], self.algorithms["two_opt"][0], self.algorithms["two_opt"][1], self.algorithms["crossover"][0], self.algorithms["crossover"][1], self.algorithms["three_op"][0], self.algorithms["three_op"][1], self.algorithms["discard"][0], self.algorithms["discard"][1])
    
    
    def __str__(self):
        header = ["Mat_Levy_Flights  Cell_Levy_Flights  Elite_Crossover  Partial_Inversion  Mutate  Two_opt  Crossover  Three_opt  Discard"]
        header += repr(self)
        header ="\n".join(header)+"\n"
        s = header
        return s
    
    ## Adds val tuples to the algorithm arg's tuples
    # @param alg str name of the algorithm selected
    # @param val tuple value to be added
    def update(self, alg, val):
        self.algorithms[alg] = (self.algorithms[alg][0]+val[0], self.algorithms[alg][1]+val[1])
    
    ## Create and open input file 
    def write(self,header=False):
        try:
            with open(self.fname, "a") as f:  
                if header:
                    f.write(str(self))
                else:
                    f.write(repr(self))

            f.close()

        except IOError as e:
            module_logger.error("I/O error({0}): {1}".format(e.errno, e.strerror))   

        # Test that the file closed
        assert f.closed==True, "File did not close properly."
        
