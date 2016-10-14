#######################################################################################################
#
# Module : Gnowee_Utilities.py
#
# Contains : Gnowee methods and the functions and methods to support Gnowee hybrid metaheuristic 
#            optimization algorithms.
#
# Author : James Bevins
#
# Last Modified: 19Aug16
#
#######################################################################################################

import logging
module_logger = logging.getLogger('Coeus.Gnowee_Utilities')

import shutil
import time
import os
import sys
sys.path.insert(0,os.path.abspath(os.getcwd())+'/Sampling')

import copy as cp
import numpy as np

from SamplingMethods import Initial_Samples
from MCNP_Utilities import MCNP_Surface, MCNP_Cell, Read_Tally_Output, Read_MCNP_Output, Print_MCNP_Input
from Utilities import Switch, to_NormDiff, RelativeLeastSquares, FuncThreadWithReturn, Event
from math import tan, radians, log
from random import random

#---------------------------------------------------------------------------------------#    
class Gnowee_Settings:
    """
    Creates an object representing the settings for the optimization algorithm
   
    Attributes
    ==========
    population : integer
        The number of parents in each generation 
        [Default: 25]
    initial_sampling : string
        The method used to sample the phase space and create the initial population 
        Valid('random','nolh','nolh-rp','nolh-cdr',and 'lhc') 
        [Default: 'LHC']
    frac_discovered : float
        Discovery probability 
        [Default: 0.25]
    frac_elite : float
        Elite fraction of population
        [Default: 0.20]
    frac_levy : float
        Fraction of Levy population engaging in Levy flights in a given generation
        [Default: 0.40]
    max_gens : integer
        The maximum number of generations to search 
        [Default: 10000]
    feval_max : integer
        The maximum number of objective function evaluations 
        [Default: 100000]
    conv_tol : float
        The minimum change of the best objective value before the search
        terminates 
        [Default: 1e-6]
    stall_iter_limit : integer
        The maximum number of genrations to search without a decrease
        exceeding conv_tol 
        [Default: 200]       
    optimal_fitness : float
        The best know fitness value for the problem considered 
        [Default: 0.0]
    opt_conv_tol : float
        The maximum deviation from the best know fitness value before the search
        terminates 
        [Default: 1e-2]
    alpha : float
        Levy exponent - defines the index of the distribution and controls scale properties of the stochastic process
        [Default: 1.5]
    gamma : float
        Gamma - Scale unit of process for Levy flights 
        [Default: 1.0]
    n : integer
        Number of independent variables - can be used to reduce Levy flight variance 
        [Default: 1]
    scaling_factor : scalar
        Step size scaling factor used to adjust Levy flights to length scale of system 
        [Default: 10]
    Returns
    =======
    None
    """
        
        
    def __init__(self,population=25,initial_sampling='lhc',frac_discovered=0.25,frac_elite=0.20, frac_levy=0.4,
                 max_gens=10000, feval_max=100000, conv_tol=1e-6, stall_iter_limit=200, optimal_fitness=0.01,
                 opt_conv_tol=1e-2,alpha=1.5, gamma=1.0,n=1,scaling_factor=10.0):          
        self.p=population                    
        self.s=initial_sampling          
        self.fd=frac_discovered            
        self.fe=frac_elite              
        self.fl=frac_levy  
        self.gm=max_gens                 
        self.em=feval_max               
        self.ct=conv_tol                 
        self.sl=stall_iter_limit        
        self.of=optimal_fitness         
        self.ot=opt_conv_tol         
        self.a=alpha                   
        self.g=gamma                   
        self.n=n                       
        self.sf=scaling_factor
        
    def __repr__(self):
        return "Gnowee Settings({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {})".format(self.p, self.s,
                self.fd, self.fe, self.fl, self.gm, self.em, self.ct, self.sl, self.of, self.ot, self.a, self.g, self.n, \
                self.sf)
    
    
    def __str__(self):
        header = ["\nGnowee Optimization Settings:"]
        header += ["Population Size = {}".format(self.p)]
        header += ["Initial sampling method = {}".format(self.s)]
        header += ["Discovery Fraction = {}".format(self.fe)]
        header += ["Elite fraction = {}".format(self.fe)]
        header += ["Levy fraction = {}".format(self.fl)]
        header += ["Maximum number of genrations = {}".format(self.gm)]
        header += ["Maximum number of function evaluations = {}".format(self.em)]
        header += ["Stall convergence tolerance = {}".format(self.ct)]
        header += ["Stall iteration limit = {}".format(self.sl)]
        header += ["Optimal fitness = {}".format(self.of)]
        header += ["Optimal convergence tolerance = {}".format(self.ot)]
        header += ["Levy exponent = {}".format(self.a)]
        header += ["Levy scale unit = {}".format(self.g)]
        header += ["Levy independent variables = {}".format(self.n)]
        header += ["Step size scaling factor = {}".format(self.sf)]
        header ="\n".join(header)+"\n"
        s = header
        return s
    
    
    def read_settings(self, filename):
        """Parses a Gnowee settings csv input file. 
        The key word options are:
            Population Size
            Initial Sampling Method
            Discovery Fraction
            Elite Fraction
            Levy Fraction
            Max Generations
            Max Function Evaluations
            Stall Convergence Tolerance
            Stall Iteration Limit
            Optimal Fitness
            Optimal Convergence Tolerance
            Levy Alpha
            Levy Gamma
            Levy Indepentent Variables
            Step Size Scaling Factor
        """

        # Open file
        try: 
            self.f = open(filename, 'r') 
            
            # Read the file line by line and store the values in the ETA_Params object
            for line in self.f:
                split_list=line.split(',')
                for case in Switch(split_list[0].strip().lower()):
                    if case('Population Size'.lower()): 
                        self.p=int(split_list[1].strip())
                        break
                    if case('Initial Sampling Method'.lower()): 
                        self.s=split_list[1].strip()
                        break
                    if case('Discovery Fraction'.lower()): 
                        self.fd=float(split_list[1].strip())
                        break
                    if case('Elite Fraction'.lower()): 
                        self.fe=float(split_list[1].strip())
                        break
                    if case('Levy Fraction'.lower()): 
                        self.fl=float(split_list[1].strip())
                        break
                    if case('Max Generations'.lower()): 
                        self.gm=int(split_list[1].strip())
                        break
                    if case('Max Function Evaluations'.lower()): 
                        self.em=int(split_list[1].strip())
                        break
                    if case('Stall Convergence Tolerance'.lower()): 
                        self.ct=float(split_list[1].strip())
                        break
                    if case('Stall Iteration Limit'.lower()): 
                        self.sl=int(split_list[1].strip())
                        break
                    if case('Optimal Fitness'.lower()): 
                        self.of=float(split_list[1].strip())
                        break
                    if case('Optimal Convergence Tolerance'.lower()): 
                        self.ot=float(split_list[1].strip())
                        break
                    if case('Levy Alpha'.lower()): 
                        self.a=float(split_list[1].strip())
                        break
                    if case('Levy Gamma'.lower()): 
                        self.g=float(split_list[1].strip())
                        break
                    if case('Levy Indepentent Variables'.lower()): 
                        self.n=int(split_list[1].strip())
                        break
                    if case('Step Size Scaling Factor'.lower()): 
                        self.sf=float(split_list[1].strip())
                        break
                    if case(): # default, could also just omit condition or 'if True'
                        module_logger.warning("A user input was found in the Gnowee settings file that does not match the allowed input types ({}): \
                            Population Size, Initial Sampling Method, Discovery Fraction, Elite Fraction, \
                            Levy Fraction, Max Generations, \
                            Max Function Evaluations, Stall Convergence Tolerance, Stall Iteration Limit, Optimal Fitness, \
                            Optimal Convergence Tolerance, Levy Alpha, Levy Gamma, Levy Indepentent Variables, \
                            Step Size Scaling Factor".format(split_list[0].strip()))
        
            # Close the file
            self.f.close()
        except IOError as e:
            module_logger.error("I/O error({0}): {1}".format(e.errno, e.strerror))    
       
        # Test that the file closed
        assert self.f.closed==True, "File did not close properly."

#-------------------------------------------------------------------------------------------------------------#  
class Parent:
    """
    Creates a parent object representing a current design
   
    Attributes
    ==========
    identifier : integer
        A set identifier tying a parent to a folder set
    fit : scalar
        The assessed design fitness
    geom : object
        An object containing the design geometry variables 
    rset : MCNP_Settings object
        An object representing the settings for running the MCNP radiation trasport code. Contains the source, physics, 
        and tally information.
    fixed_mats : integer
        Number of fixed materials in the geometry.  This accounts for structural materials and foils that shouldn't change.
    Returns
    =======
    None
    """
    
    def __init__(self,identifier,eta,geometry,GS,mcnp,mats,ex,i=0,fitness=1E15,build_geom=True):
        """

        Parameters
        ==========
        identifier : integer
            A set identifier tying a parent to a folder set
        eta : ETA parameters object
            An object that contains all of the constraints required to initialize the geometry
        geometry : MCNP_Geometry object
            The geometry for running the MCNP radiation trasport code. Contains the surfaces, cells, and material information
        GS : Gnowee Settings object
            An object representing the settings for the optimization algorithm
        mcnp : MCNP_Settings object
            An object representing the settings for running the MCNP radiation trasport code. Contains the source, physics, 
            and tally information.
        mats : dict of material objects
            A dictionary of the material objects from which ETA materials can be selected
            [Default: {}]
        ex : list    
            A list of materials to be excluded

        Optional
        ========
        i : integer
            Parent indexed location for initial LHC sampling purposes
        fitness : float
            The assessed design fitness
        build_geom : boolean
            Used to indicate if the geometry needs to be build for a new parent.  
            Turned off if the complete geometry is passed in.

        Returns
        =======
        """

        # Initialize the parent
        self.ident=identifier
        self.geom=cp.deepcopy(geometry)
        self.rset=cp.deepcopy(mcnp)
        self.fit=fitness 
        
        # After initial initialization set sampling method to random
        GS.s='random'
        
        # Determine number of fixed mats
        self.fixed_mats=1    # ETA structural material
        tmp=[]
        tmp.append(eta.struct_mat)
        if eta.ds_mat not in tmp:
            self.fixed_mats+=1
            tmp.append(eta.ds_mat)
        if eta.nas_mat not in tmp:
            self.fixed_mats+=1
            tmp.append(eta.nas_mat)
        if eta.toad_mat not in tmp:
            self.fixed_mats+=1
            tmp.append(eta.toad_mat)
        for m in eta.nas_mat_f:
            if m not in tmp:
                self.fixed_mats+=1
                tmp.append(m)
        for m in eta.toad_mat_f:
            if m not in tmp:
                self.fixed_mats+=1
                tmp.append(m)   
        if eta.holder_mat not in tmp:
            self.fixed_mats+=1
            tmp.append(eta.holder_mat) 
        if eta.h_fill_mat not in tmp:
            self.fixed_mats+=1
            tmp.append(eta.h_fill_mat) 
        
        if build_geom==True:
            foil_cells=[]
        
            # Build Master Key list for materials library
            master_keys=mats.keys()

            # Randomly sample the NAS location.  The sampled dimensions are [axial, vertical macros, horizontal macros]
            lb=np.array([max(eta.tcc_dist+eta.t_w+2*eta.t_nas+sum(eta.t_nas_f)+0.203,eta.tcc_dist+eta.t_w+(eta.r_nas-eta.r_f+eta.t_w)/tan(radians(90-eta.theta))),0])
            ub=np.array([eta.snout_dist-eta.t_c-2*eta.t_nas-sum(eta.t_nas_f)-0.203,1])
            samples=Initial_Samples(lb,ub,GS.s,GS.p)
        
            # Determine location of NAS and first foil
            z_nas=samples[i,0]
            del_r_nas=eta.r_nas
            del_z_nas=2*eta.t_nas+sum(eta.t_nas_f)+0.203   #0.203cm is Toad thickness
            z=samples[i,0]+eta.t_nas
                
            # Build Surfaces and Cells for each foil
            for f in range(0,len(eta.nas_mat_f)):
                # Determine del_z and del r
                del_z=eta.t_nas_f[f]
                del_r=eta.r_nas_f
                
                # Find index in mat lib for foil material
                ind=next((j for j, item in enumerate(self.geom.matls) if item == eta.nas_mat_f[f]), -1)

                # Add foil surfaces and cells
                foil_cells.append(self.geom.cells[-1].name+1)
                self.geom.add_surf(MCNP_Surface(self.geom.surfaces[-1].name+1,"RCC",\
                                                vx=0.0,vy=0.0,vz=z,\
                                                hx=0.0,hy=0.0,hz=del_z,\
                                                r=del_r,comment="{} foil".format(eta.nas_mat_f[f])))
                self.geom.add_cell(MCNP_Cell(foil_cells[-1],ind+1,"mass",mats[self.geom.matls[ind]].density,\
                                             "-{}".format(self.geom.surfaces[-1].name), \
                                             (1,0), "{} foil".format(eta.nas_mat_f[f])))

                # Update currentz location
                z+=eta.t_nas_f[f]

                # Check to see if Toad follows
                if eta.toad_loc.strip()==eta.nas_mat_f[f].strip():
                    # Determine toad and first foil locations
                    z_toad=z
                    del_r=eta.r_toad
                    del_r_toad=2.5
                    del_z_toad=0.203
                    z+=0.05

                    # Add TOAD foils
                    for t in range(0,len(eta.toad_mat_f)):
                        # Determine del_z
                        del_z=eta.t_toad[t]

                        # Find index in mat lib for foil material
                        ind=next((j for j, item in enumerate(self.geom.matls) if item == eta.toad_mat_f[t]), -1)

                        # Add foil surfaces and cells
                        foil_cells.append(self.geom.cells[-1].name+1)
                        self.geom.add_surf(MCNP_Surface(self.geom.surfaces[-1].name+1,"RCC",\
                                                       vx=0.0,vy=0.0,vz=z,\
                                                       hx=0.0,hy=0.0,hz=del_z,\
                                                       r=del_r,comment="{} foil".format(eta.toad_mat_f[t])))
                        self.geom.add_cell(MCNP_Cell(foil_cells[-1],ind+1,"mass",mats[self.geom.matls[ind]].density,\
                                                    "-{}".format(self.geom.surfaces[-1].name), \
                                                     (1,0), "{} foil".format(eta.toad_mat_f[t])))

                        # Update current z location
                        z+=eta.t_toad[t]
                
                    # Find index in mat lib for TOAD casing material
                    ind=next((j for j, item in enumerate(self.geom.matls) if item == eta.toad_mat), -1)

                    # Exclude TOAD Foils
                    exclude="" 
                    for s in foil_cells[-len(eta.toad_mat_f):]:
                        exclude+=" #{}".format(s) 

                    # Add TOAD surface and cell   
                    foil_cells.append(self.geom.cells[-1].name+1)
                    toad_surf=self.geom.surfaces[-1].name+1
                    self.geom.add_surf(MCNP_Surface(toad_surf,"RCC",\
                                                    vx=0.0,vy=0.0,vz=z_toad,\
                                                    hx=0.0,hy=0.0,hz=del_z_toad,\
                                                    r=del_r_toad,comment="TOAD"))
                    self.geom.add_cell(MCNP_Cell(foil_cells[-1],ind+1,"mass",mats[self.geom.matls[ind]].density,\
                                                 "-{} {}".format(self.geom.surfaces[-1].name, exclude), \
                                                 (1,0), "TOAD"))

                    # Update current z location
                    z=z_toad+del_z_toad

            # Find index in mat lib for NAS casing material
            ind=next((j for j, item in enumerate(self.geom.matls) if item == eta.nas_mat), -1)

            # ADD NAS Casing
            # Exclude NAS Foils
            exclude="" 
            for s in foil_cells:
                exclude+=" #{}".format(s) 

            # Add NAS surface and cell   
            foil_cells.append(self.geom.cells[-1].name+1)
            self.geom.add_surf(MCNP_Surface(self.geom.surfaces[-1].name+1,"RCC",\
                                            vx=0.0,vy=0.0,vz=z_nas,\
                                            hx=0.0,hy=0.0,hz=del_z_nas,\
                                            r=del_r_nas,comment="NAS"))
            self.geom.add_cell(MCNP_Cell(foil_cells[-1],ind+1,"mass",mats[self.geom.matls[ind]].density,\
                                         "-{} {}".format(self.geom.surfaces[-1].name, exclude), \
                                         (1,0), "NAS"))   

            # Add holder fill
            ind=next((j for j, item in enumerate(self.geom.matls) if item == eta.h_fill_mat), -1)
            self.geom.add_surf(MCNP_Surface(self.geom.surfaces[-1].name+1,"RCC",\
                                            vx=0.0,vy=0.0,vz=z_nas,\
                                            hx=0.0,hy=0.0,hz=del_z_nas,\
                                            r=eta.r_o-eta.t_w,comment="Holder Fill"))
            self.geom.add_cell(MCNP_Cell(self.geom.cells[-1].name+1,ind+1,"mass",mats[self.geom.matls[ind]].density,\
                                         "(-{0} {1} -{2}):(-{0} {1} -{3})".format(self.geom.surfaces[-1].name, \
                                         self.geom.surfaces[-2].name,  self.geom.surfaces[2].name,  self.geom.surfaces[4].name),\
                                         (1,0), "Holder Fill"))   
        
            # Add holder 
            ind=next((j for j, item in enumerate(self.geom.matls) if item == eta.holder_mat), -1)
            holder_surf=self.geom.surfaces[-1].name+1
            self.geom.add_surf(MCNP_Surface(holder_surf,"RCC",\
                                            vx=0.0,vy=0.0,vz=z_nas-eta.t_h,\
                                            hx=0.0,hy=0.0,hz=del_z_nas+2*eta.t_h,\
                                            r=eta.r_o-eta.t_w,comment="Holder"))
            self.geom.add_cell(MCNP_Cell(self.geom.cells[-1].name+1,ind+1,"mass",mats[self.geom.matls[ind]].density,\
                                         "(-{0} {1} -{2}):(-{0} {1} -{3})".format(self.geom.surfaces[-1].name, \
                                         self.geom.surfaces[-2].name,  self.geom.surfaces[2].name,  self.geom.surfaces[4].name),\
                                         (1,0), "Holder"))   
        

            # Place the vertical and horizontal layers
            mat_lb=np.ones(eta.max_vert+eta.max_horiz)
            mat_ub=np.ones(eta.max_vert+eta.max_horiz)*(len(mats)-1)
            mat_samp=Initial_Samples(mat_lb,mat_ub,'random',1)
            mat_samp=np.round(mat_samp,0).astype(int)
            mat_samp=mat_samp.ravel()
        
            # Place Vertical layers
            exclude=""
            vert_cells=[]
            for l in range(eta.max_vert):
                # Sample Starting Location
                eta_lb=np.array([0, eta.tcc_dist+eta.t_w])                        # [r,z]
                eta_ub=np.array([eta.r_o-eta.t_w, eta.snout_dist-eta.t_c])  
                loc_samp=Initial_Samples(eta_lb,eta_ub,"random",1) 
                loc_samp=loc_samp.ravel()

                # Sample delta r and z
                lb=np.array([0,0])                        # [del_r,del_z]
                ub=np.array([eta.r_o-eta.t_w-loc_samp[0], eta.snout_dist-eta.t_c-loc_samp[1]])  
                del_samp=Initial_Samples(lb,ub,"random",1) 
                del_samp=del_samp.ravel()

                # Check Boundries 
                if any(loc_samp+del_samp>eta_ub):
                    module_logger.warning("Vertical layer placement is outside of ETA boundaries.")

                # Exclude materials in exlusion list
                while master_keys[mat_samp[l]] in ex:
                    tmp=Initial_Samples(mat_lb,mat_ub,'random',1).ravel()
                    mat_samp[l]=np.round(tmp[l],0).astype(int)

                # Add material
                self.geom.add_matls(mats,master_keys[mat_samp[l]])
                ind=next((i for i, item in enumerate(self.geom.matls) if item == master_keys[mat_samp[l]]), -1)

                # Create Surfaces
                self.geom.add_surf(MCNP_Surface(self.geom.surfaces[-1].name+1,"RCC",\
                                                vx=0.0,vy=0.0,vz=loc_samp[1],\
                                                hx=0.0,hy=0.0,hz=del_samp[1],\
                                                r=loc_samp[0],comment="vert #{}".format(l+1)))
                self.geom.add_surf(MCNP_Surface(self.geom.surfaces[-1].name+1,"RCC",\
                                                vx=0.0,vy=0.0,vz=loc_samp[1],\
                                                hx=0.0,hy=0.0,hz=del_samp[1],\
                                                r=del_samp[0]+loc_samp[0],comment="vert #{}".format(l+1)))

                # Create Cell
                vert_cells.append(self.geom.cells[-1].name+1) 
                self.geom.add_cell(MCNP_Cell(vert_cells[-1],ind+1,"mass",mats[self.geom.matls[ind]].density,\
                                            "({0} -{1} -{2} {4} {5}):({0} -{1} -{3} {4}\
                                            {5})".format(self.geom.surfaces[-2].name, \
                                            self.geom.surfaces[-1].name, self.geom.surfaces[2].name, \
                                            self.geom.surfaces[4].name, holder_surf, exclude), (1,0), comment='vert'))   
                exclude+=" #{}".format(vert_cells[-1])


            # Place Horizontal layers
            z=eta.tcc_dist+eta.t_ds
            for l in range(eta.max_vert,eta.max_vert+eta.max_horiz):         
                # Sample delta z 
                lb=np.array([0])                        # [del_z]
                ub=np.array([eta.snout_dist-eta.t_c-z])  
                del_z=Initial_Samples(lb,ub,"random",1) 
                del_z=del_z.ravel()*(1.0/(eta.max_vert+eta.max_horiz-l))

                # Check Boundries 
                if z+del_z[0]>eta.snout_dist-eta.t_c:
                    module_logger.warning("Horizontal layer placement is outside of ETA boundaries.")

                # Exclude materials in exlusion list
                while master_keys[mat_samp[l]] in ex:
                    tmp=Initial_Samples(mat_lb,mat_ub,'random',1).ravel()
                    mat_samp[l]=np.round(tmp[l],0).astype(int)

                # Add material
                self.geom.add_matls(mats,master_keys[mat_samp[l]])
                ind=next((i for i, item in enumerate(self.geom.matls) if item == master_keys[mat_samp[l]]), -1)

                # Create Surfaces
                if l==eta.max_vert:
                    self.geom.add_surf(MCNP_Surface(self.geom.surfaces[-1].name+1,"PZ",d=z,comment="horiz start"))
                self.geom.add_surf(MCNP_Surface(self.geom.surfaces[-1].name+1,"PZ",d=z+del_z[0],comment="horiz #{}".format(l-eta.max_vert+1)))

                # Update current z location
                z+=del_z[0]

                # Create Cell            
                self.geom.add_cell(MCNP_Cell(self.geom.cells[-1].name+1,ind+1,"mass",\
                                            mats[self.geom.matls[ind]].density,\
                                            "({0} -{1} -{2} {4} {5}):({0} -{1} -{3} {4} \
                                            {5})".format(self.geom.surfaces[-2].name,\
                                            self.geom.surfaces[-1].name, self.geom.surfaces[2].name,\
                                            self.geom.surfaces[4].name, holder_surf, exclude), (1,0), comment='horiz'))  
           
        
        # Set tallies
        m=next((j for j, item in enumerate(self.geom.matls) if item == eta.fissile_mat), -1) +1
        c=next((j for j, item in enumerate(self.geom.cells) if item.m == m), -1)    
        c=self.geom.cells[c].name
        self.rset.set_tallies(c, m)
        
    def __repr__(self):
        return "Parent Design Object({}, {}, {}, {})".format(self.ident, self.geom, self.fit, self.rset)
    
    
    def __str__(self):
        header = ["\nParent Design Object:"]
        header += ["Identifier = {}".format(self.ident)]
        header += ["Fitness = {}".format(self.fit)]
        header += ["Geometry = {}".format(self.geom)]
        header += ["NPS = {}".format(self.rset)]
        header ="\n".join(header)+"\n"
        s = header
        return s
       
#-------------------------------------------------------------------------------------------------------------#  
def Calc_Fitness(ids, pop, obj, min_fiss=0, max_w=1000):
    """
    Print the generated MCNP input deck to file 
   
    Parameters
    ==========
    ids : list of integers
        The parents that need to have fitness solutions calculated
    pop : list of parent objects
        The population and their design features
    obj : array
        The objective spectrum

    Optional
    ========
    min_fiss : float
        A constraint specifying the minimum number fo fissions. Implemented as a soft constraint.
        [Default = 0]
    max_w : float
        A constraint specifying the maximum weight of the assembly.  Implemented as a hard constraint.
        
    Returns
    =======
    None
    """  
    rundir=os.path.abspath(os.path.join(os.path.abspath(os.getcwd()),os.pardir))+"/Results/Population/"
    
    for i in ids:
        index = next((c for c, parent in enumerate(pop) if parent.ident == i), -1)
        (tally,fissions,weight)=Read_MCNP_Output(rundir+str(i)+'/tmp/ETA.out', '24', '14')
        tally=to_NormDiff(tally)
        tmp_fit=RelativeLeastSquares(tally,obj)
        module_logger.debug("Parent ID # {} has fitness = {} from RLS".format(i,tmp_fit))
        
        # Check constraints
        weight=weight/1000         # conversion to kg
        if fissions[0] == 0.0:
            module_logger.warning("WARNING: No fissions occured for the ETA design in parent #{}".format(i))
            tmp_fit+=1E15
        elif fissions[0] > 0 and fissions[0] < min_fiss:
            tmp_fit+= 0.1*(min_fiss/fissions[0]-1)
        elif fissions[0] > min_fiss:
            tmp_fit-= 0.01*(fissions[0]/min_fiss-1)
        module_logger.debug("fissions[0] = {} and min_fiss = {} ".format(fissions[0],min_fiss))
        module_logger.debug("Parent ID # {} has fitness = {} from RLS+fissions".format(i,tmp_fit))
        
        if weight > max_w:
            tmp_fit+=1E15
        module_logger.info("Parent ID # {} has fitness = {} from RLS+fissions+weight".format(i,tmp_fit))
        
        # Save fitness
        pop[index].fit=tmp_fit  
    
#-------------------------------------------------------------------------------------------------------------#  
def Pop_Update(old, new, nps, eta=None, mats=None, run=None, rr=False):
    """
    Updates the population based on the assessed fitness values.  
   
    Parameters
    ==========
    old : list of parent objects
        The current population and their design features
    new : list of parent objects
        The proposed population and their design features
    nps : float
        The baseline number of NPS particles specified 
    Optional
    ========
    eta : ETA parameters object
        An object that contains all of the constraints required to initialize the geometry
    mats : dict of material objects
        A dictionary of the material objects from which ETA materials can be selected
    run : function
        A function that runs the radiation transport calculations
    rr : boolean
        Indicator if random replacement is used to update the list
        
    Returns
    =======
    changes : Integer
        The number of accepted changes
    feval : Integer
        The number of function evaluations performed for increasing the particles
    """
    
    #Test input values for consistency
    if run != None:
        assert hasattr(run, '__call__'), 'Invalid function handle'
    
    # Initialize lists for updated nps
    ids_1E8=[]
    ids_1E7=[]
    ind_1E8=[]
    ind_1E7=[]
    
    changes=0
    path=os.path.abspath(os.path.join(os.path.abspath(os.getcwd()),os.pardir))+"/Results/Population/"
    
    for i in range(0,len(new)):
        # Determine appropriate comparison index
        if rr==False:
            ind = next((c for c, parent in enumerate(old) if parent.ident == new[i].ident), -1)
            old_ident=old[ind].ident
        elif rr==True:
            ind=int(random()*len(old))
            old_ident=old[ind].ident
            
        # Compare the old and new parent.  Replace if the new one is better
        if new[i].fit<old[ind].fit:
            changes+=1
            old[ind]=cp.deepcopy(new[i])
            old[ind].ident=old_ident
            if old[ind].fit <= 0.135 and old[ind].rset.nps<nps*100:
                old[ind].rset.nps=nps*100
                ids_1E8.append(old[ind].ident)
                ind_1E8.append(ind)
            elif old[ind].fit <= 0.5 and old[ind].rset.nps<nps*10:
                old[ind].rset.nps=nps*10
                ids_1E7.append(old[ind].ident)
                ind_1E7.append(ind)
    
            # Save the output file
            else:
                shutil.copyfile(path+str(old[ind].ident)+'/tmp/ETA.out', path+str(old[ind].ident)+'/ETA.out')

    if eta != None and mats != None and run != None:
        if len(ids_1E7)!=0:
            module_logger.info("For 1E7, the file ids are = {}".format(ids_1E7))
            for i in ind_1E7:
                Print_MCNP_Input(eta,old[i].geom,old[i].rset,mats,old[i].ident,adv_print=True)
            run(ids_1E7,[1E7]*len(ids_1E7),code='mcnp6.mpi')
            Calc_Fitness(ids_1E7, old, eta.spectrum[:,1], eta.min_fiss, eta.max_weight)
        if len(ids_1E8)!=0:
            module_logger.info("For 1E8, the file ids are = {}".format(ids_1E8))
            for i in ind_1E8:
                Print_MCNP_Input(eta,old[i].geom,old[i].rset,mats,old[i].ident,adv_print=True)
            run(ids_1E8,[1E8]*len(ids_1E8),code='mcnp6.mpi')
            Calc_Fitness(ids_1E8, old, eta.spectrum[:,1], eta.min_fiss, eta.max_weight)
            
        for i in ids_1E7+ids_1E8:
            shutil.copyfile(path+str(i)+'/tmp/ETA.out',path+str(i)+'/ETA.out')
    
    return changes,len(ids_1E7)+len(ids_1E8)

#-------------------------------------------------------------------------------------------------------------#  
class Timeline():
    """
    An object that stores event objects to track optimization progress.
   
    Attributes
    ==========
    tline : list of event objects
        Contains a list of event objects detailing the optimization history
    fname : str
        Name and path of the file to store the timeline for post processing

    Optional
    ========
        
    Returns
    =======
    None    
    """

    def __init__(self,tline=[],fname=os.path.abspath(os.path.join(os.getcwd(), os.pardir))+"/Results/timeline.txt"):
        self.tline=tline
        self.fname=fname
        if os.path.isfile(self.fname)==True:
            os.remove(self.fname)
        
    def __repr__(self):
        return "{}  {}  {:6e}  {:.2e}  {}\n".format(self.tline[-1].g, self.tline[-1].e, self.tline[-1].f, self.tline[-1].n, self.tline[-1].i)
    
    
    def __str__(self):
        header = ["\nGenerations  Evaluations  Fitness  NPS  ID"]
        header ="\n".join(header)+"\n"
        s = header + "\n".join(["{:5d}  {:5d}  {:8.4f}  {:.2e}  {}".format(t.g, t.e, t.f, t.n, t.i) for t in self.tline])
        return s
    
    def update(self, pop, gen, feval):
        #Sort the population according to fitness
        pop.sort(key=lambda x: x.fit)
        
        module_logger.info("\nAfter sorting:")
        for n in pop:
            module_logger.info("Parent ID # {} has fitness = {}".format(n.ident,n.fit))
        module_logger.info("\n".format(n.ident,n.fit))    
        
        # Check for results folders existence
        if os.path.exists('{}/Results/History'.format(path))==False:    
            os.mkdir(('{}/Results/History'.format(path))
    
        #Store history on timeline if new optimal design found
        path=os.path.abspath(os.path.join(os.getcwd(), os.pardir))
        if len(self.tline)<1:
            self.tline.append(Event(0,feval,pop[0].fit,pop[0].rset.nps,pop[0].ident)) 
            self.write()
            # Save the input file in results folder
            shutil.copyfile(path+'/Results/Population/'+str(pop[0].ident)+'/ETA.inp',\
                            path+'/Results/History/ETA_{}.inp'.format(self.tline[-1].e))
            # Save the output file in results folder
            shutil.copyfile(path+'/Results/Population/'+str(pop[0].ident)+'/ETA.out',\
                            path+'/Results/History/ETA_{}.out'.format(self.tline[-1].e))
            # Save the wwinp file in results folder
            if os.path.isfile(path+'/Results/Population/'+str(pop[0].ident)+'/wwinp'):
                shutil.copyfile(path+'/Results/Population/'+str(pop[0].ident)+'/wwinp',\
                            path+'/Results/History/wwinp')
        elif pop[0].fit< self.tline[-1].f:
            self.tline.append(Event(self.tline[-1].g+gen,self.tline[-1].e+feval,pop[0].fit,pop[0].rset.nps,\
                                    pop[0].ident))
            self.write()
            # Save the input file in results folder
            shutil.copyfile(path+'/Results/Population/'+str(pop[0].ident)+'/ETA.inp',\
                            path+'/Results/History/ETA_{}.inp'.format(self.tline[-1].e))
            # Save the output file in results folder
            shutil.copyfile(path+'/Results/Population/'+str(pop[0].ident)+'/ETA.out',\
                            path+'/Results/History/ETA_{}.out'.format(self.tline[-1].e))
            # Save the wwinp file in results folder
            if os.path.isfile(path+'/Results/Population/'+str(pop[0].ident)+'/wwinp'):
                shutil.copyfile(path+'/Results/Population/'+str(pop[0].ident)+'/wwinp',\
                            path+'/Results/History/wwinp')
        else:
            self.tline[-1].e+=feval
            self.tline[-1].g+=gen
        
        return pop
    
    def write(self):
        # Create and open input file 
        try:
            with open(self.fname, "a") as f:  
                f.write(repr(self))

            # Close the file
            f.close()

        except IOError as e:
            module_logger.error("I/O error({0}): {1}".format(e.errno, e.strerror))   

        # Test that the file closed
        assert f.closed==True, "File did not close properly."
        
#---------------------------------------------------------------------------------------# 
def Rejection_Bounds(parent,child,stepsize,lb,ub,S,change_count=0):
    """
    Application of problem boundaries to generated solutions. If a solution is outside of the 
    bounds, the step is rejected and that particular value reverts to the previous solution.
   
    Parameters
    ==========
    parent : array
        The current system designs
    child : array
        The proposed new system designs
    stepsize : float
        The Levy flight stepsize
    lb : array
        The lower bounds of the design variable(s)
    ub : array
        The upper bounds of the design variable(s)
    S : Gnowee Settings Object    
        An object representing the settings for the Gnowee optimization algorithm
   
    Optional
    ======== 
    change_count : integer
        Counter to track the number of solutions that occur outside of problem boundaries.  
        Can be used to diagnose too large or small of alpha
        (Default: 0)
   
    Returns
    =======
    new : array
        The new system designs that are within problem boundaries
    """
    
    assert len(child)==len(lb), 'Child and lb best have different # of design variables in Rejection_Bounds function.'
    assert len(ub)==len(lb), 'Boundaries best have different # of design variables in Rejection_Bounds function.'
            
    for i in range(0,len(child),1): 
        change_count=0
        while child[i]<lb[i] or child[i]>ub[i]:
            if change_count >8:
                module_logger.info("Stubborn Child:{},{},{},{},{}".format(child[i], lb[i], ub[i], child[i]<lb[i], child[i]>ub[i]))
                sys.exit()
            elif change_count >=6: 
                if ub[i] < 0:
                    ub[i]=abs(ub[i])
                child=Simple_Bounds(child,lb,ub)
                change_count+=1
            elif change_count >= 5:
                child[i]=cp.copy(parent[i])
                change_count+=1
            else:
                stepsize[i]=stepsize[i]/2.0  
                child[i]=child[i]-stepsize[i]
                change_count+=1
        if child[i]<0.0:
            module_logger.info("Negative Child:{},{},{},{},{}".format(child[i], lb[i], ub[i], child[i]<lb[i], child[i]>ub[i]))
            sys.exit()
    module_logger.debug("Change Count = {}".format(change_count)) 
    return child

#---------------------------------------------------------------------------------------#
def Simple_Bounds(tmp,lb,ub,change_count=0):
    """
    Application of problem boundaries to generated solutions
   
    Parameters
    ==========
    tmp : array
        The proposed new system designs
    lb : array
        The lower bounds of the design variable(s)
    ub : array
        The upper bounds of the design variable(s)
   
    Optional
    ========
    change_count : integer
        Counter to track the number of solutions that occur outside of problem boundaries.  
        Can be used to diagnose too large or small of alpha
        (Default: 0)
   
    Returns
    =======
    tmp : array
        The new system designs that are within problem boundaries
    """
    
    assert len(tmp)==len(lb), 'Tmp and lb best have different # of design variables in Simple_Bounds function.'
    assert len(ub)==len(lb), 'Boundaries best have different # of design variables in Simple_Bounds function.'
            
    #Check consistency of bounds
    for i in range(len(lb)):
        if lb[i]>ub[i] and lb[i]>0.0 and ub[i]>0.0:
            tb=lb[i]
            lb[i]=cp.copy(ub[i])
            ub[i]=tb
        if ub[i]<0.0:
            ub[i]=cp.copy(lb[i])
        if lb[i]<0.0:
            lb[i]=0.0000001
            
    #Apply bounds; update to boundary if out of bounds
    for i in range(0,len(tmp),1): 
        if tmp[i]<lb[i]:
            module_logger.debug("tmp[{0}] = {1}, lb[{0}] = {2} ".format(i,tmp[i],lb[i]))
            tmp[i]=cp.copy(lb[i])
            module_logger.debug("Changing LB: {}, {}".format(tmp[i], lb[i]))
            change_count+=1
            
        elif tmp[i]>ub[i]:
            module_logger.debug("tmp[{0}] = {1}, ub[{0}] = {2} ".format(i,tmp[i],ub[i]))
            tmp[i]=cp.copy(ub[i])
            module_logger.debug("Changing UB: {}, {}".format(tmp[i], ub[i]))
            change_count+=1
        if tmp[i]<0.0:
            module_logger.info("Negative tmp[{}]: {},{},{},{},{}".format(i, tmp[i], lb[i], ub[i], tmp[i]<lb[i], tmp[i]>ub[i]))
    
   
        module_logger.debug("Change Count = {}".format(change_count)) 
    return tmp
