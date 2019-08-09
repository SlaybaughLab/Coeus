#######################################################################################################
#
# Module : ADVANTG_Utilities.py
#
# Contains : Routines to build input decks for the ADVANTG radiation transport code.  
#
# Author : James Bevins
#
# Last Modified: 17Aug16
#
#######################################################################################################

import logging
module_logger = logging.getLogger('Coesh.ADVANTG_Utilities')

import os

import Utilities as util
import numpy as np
import copy as cp
import multiprocessing as mp

from math import ceil, floor, sqrt

#---------------------------------------------------------------------------------------#    
class ADVANTG_Settings:
        
    ## Creates a object representing the settings for running the ADVANTG deterministic radiation trasport code calculations.
    def __init__(self, lib="dplus", method="cadis", outputs="mcnp silo", tnum=24, pt_src="True", mix_tol=0.01, pn_order=1, eta_x=0.5, eta_y=0.5, eta_z =0.5, foil_x=0.25, foil_y=0.25, foil_z=0.05, ext_spacing=1.0):  
        ## string The multi-group library used
        # [Default: "dplus"]
        self.lib=lib
        ## string The solution method for ADVANTG (CADIS or DX)
        #  [Default: "cadis"]
        self.method=method
        ## string The output files to be produced
        # [Default: "mcnp"]
        self.outputs=outputs
        ## int The tally number for calculating the adjoint flux
        # [Default: 24]
        self.tnum=tnum
        ## string Whether or not the source should be treated as a point source in deterministic transport
        # [Default: False]
        self.pt_src=pt_src
        ## string The material mix tolerance fraction.  Controls the precision of mixed cells.
        # [Default: 0.01]
        self.mix_tol=mix_tol
        ## integer The scattering order
        # [Default: 1]
        self.pn_order=pn_order
        ## float The spacing of the mesh intervals in the x (radial) axis in cm in the ETA.  
        # [Default: 0.5 cm]
        self.eta_x=eta_x
        ## float The spacing of the mesh intervals in the y (radial) axis in cm in the ETA.  
        # [Default: 0.5 cm]
        self.eta_y=eta_y
        ## float The spacing of the mesh intervals in the z (axial) axis in cm in the ETA.  
        # [Default: 0.5 cm]
        self.eta_z=eta_z
        ## float The spacing of the mesh intervals in the x (radial) axis in cm near the foil.  
        # [Default: 0.25 cm]
        self.foil_x=foil_x
        ## float The spacing of the mesh intervals in the y (radial) axis in cm near the foil.  
        # [Default: 0.25 cm]
        self.foil_y=foil_y
        ## float The spacing of the mesh intervals in the z (axial) axis in cm near the foil.  
        # [Default: 0.05 cm]
        self.foil_z=foil_z
        ## float The spacing of the mesh intervals in x,y, and z axis in cm outside the ETA.  
        # [Default: 1 cm]
        self.ext=ext_spacing          
        
        
    def __repr__(self):
        return "ADVANTG Settings({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {})".format(self.lib,self.method,self.outputs,self.tnum,self.pt_src,self.mix_tol,self.pn_order,self.eta_x,self.eta_y,self.eta_z,self.foil_x,self.foil_y,self.foil_z,self.ext)
    
    
    def __str__(self):
        header = ["\nADVANTG Program Settings:"]
        header += ["Multi-Group Library = {}".format(self.lib)]
        header += ["Solution Method = {}".format(self.method)]
        header += ["Outputs = {}".format(self.outputs)]
        header += ["Adjoint Tally Number = {}".format(self.tnum)]
        header += ["Force Point Source = {}".format(self.pt_src)]
        header += ["Material Mix Tolerance = {}".format(self.mix_tol)]
        header += ["Scattering Order = {}".format(self.pn_order)]
        header += ["ETA X Spacing Interval = {}".format(self.eta_x)]
        header += ["ETA Y Spacing Interval = {}".format(self.eta_y)]
        header += ["ETA Z Spacing Interval = {}".format(self.eta_z)]
        header += ["Foil X Spacing Interval = {}".format(self.foil_x)]
        header += ["Foil Y Spacing Interval = {}".format(self.foil_y)]
        header += ["Foil Z Spacing Interval = {}".format(self.foil_z)]
        header += ["External Spacing Interval = {}".format(self.ext)]
        header ="\n".join(header)+"\n"
        s = header
        return s
    
    ## Parses a ADVANTG settings csv input file. 
    # The key word options are:
    #     Library
    #     Method
    #     Outputs
    #     Tally Number
    #     Point Source
    #     Material Mix Tolerance
    #     Scattering Order
    #     ETA X Spacing Interval
    #     ETA Y Spacing Interval
    #     ETA Z Spacing Interval
    #     Foil X Spacing Interval
    #     Foil Y Spacing Interval
    #     Foil Z Spacing Interval
    #     External Spacing Interval
    def read_settings(self, filename):
    
        # Open file
        try: 
            self.f = open(filename, 'r') 
            
            # Read the file line by line and store the values in the Advantg_Settings object
            for line in self.f:
                split_list=line.split(',')
                for case in util.Switch(split_list[0].strip().lower()):
                    if case('Library'.lower()): 
                        self.lib=split_list[1].strip()
                        break
                    if case('Method'.lower()): 
                        self.method=split_list[1].strip()
                        break
                    if case('Outputs'.lower()): 
                        self.outputs=split_list[1].strip()
                        break
                    if case('Tally Number'.lower()): 
                        self.tnum=int(split_list[1].strip())
                        break
                    if case('Point Source'.lower()): 
                        self.pt_src=split_list[1].strip()
                        break
                    if case('Material Mix Tolerance'.lower()): 
                        self.mix_tol=float(split_list[1].strip())
                        break
                    if case('Scattering Order'.lower()): 
                        self.pn_order=int(split_list[1].strip())
                        break
                    if case('ETA X Spacing Interval'.lower()): 
                        self.eta_x=float(split_list[1].strip())
                        break
                    if case('ETA Y Spacing Interval'.lower()): 
                        self.eta_y=float(split_list[1].strip())
                        break
                    if case('ETA Z Spacing Interval'.lower()): 
                        self.eta_z=float(split_list[1].strip())
                        break
                    if case('Foil X Spacing Interval'.lower()): 
                        self.foil_x=float(split_list[1].strip())
                        break
                    if case('Foil Y Spacing Interval'.lower()): 
                        self.foil_y=float(split_list[1].strip())
                        break
                    if case('Foil Z Spacing Interval'.lower()): 
                        self.foil_z=float(split_list[1].strip())
                        break
                    if case('External Spacing Interval'.lower()): 
                        self.ext=float(split_list[1].strip())
                        break
                    if case('/'): 
                        break
                    if case(): # default, could also just omit condition or 'if True'
                        module_logger.warning("A user input was found in the PartiSn settings file that does not match the allowed input types ({}) : Library,Method,Outputs,Tally Number,Point Source,Material Mix Tolerance,Scattering Order,ETA X Spacing Interval,ETA Y Spacing Interval,ETA Z Spacing Interval,Foil X Spacing Interval,Foil Y Spacing Interval,Foil Z Spacing Interval,External Spacing Interval".format(split_list[0].strip()))
        
            # Close the file
            self.f.close()
        except IOError as e:
            module_logger.error("I/O error({0}): {1}".format(e.errno, e.strerror))   
            module_logger.error("File not found was: {0}".format(filename)) 
       
        # Test that the file closed
        assert self.f.closed==True, "File did not close properly."
        
## Print the generated MCNP input deck to file 
# @param eta [ETA parameters object] An object that contains all of the constraints required to initialize the geometry
# @param geom [MCNP_Geometry object] The geometry for running the MCNP radiation trasport code. Contains the surfaces, cells, and material information
# @param S [ADVANTG_Settings object] An object representing the settings for running the ADVANTG radiation trasport code.  
# @param num integer The current cuckoo number being generated
# @param cluster boolean (optional) An indicator to change the file to run on a cluster using Run_Transport function and slurm job submission   
def Print_ADVANTG_Input(eta,geom,S,num,cluster=False):
 
    path=os.path.abspath(os.path.join(os.path.abspath(os.getcwd()), os.pardir))
        
    # Delete previous input file if present
    if os.path.exists("{}/Results/Population/{}/".format(path,num)):
        if os.path.isfile("{}/Results/Population/{}/runCADIS.adv".format(path,num)):
            os.remove("{}/Results/Population/{}/runCADIS.adv".format(path,num))
    else:
        os.mkdir("{}/Results/Population/{}/".format(path,num))

    # Create and open input file 
    try:
        with open("{}/Results/Population/{}/runCADIS.adv".format(path,num), "w") as inp_file:  

            # Write the run control information
            inp_file.write("{:25s} {}\n".format("model","mcnp"))
            inp_file.write("{:25s} {}\n".format("method",S.method))
            inp_file.write("{:25s} {}\n".format("outputs",S.outputs))
            inp_file.write("\n")
            inp_file.write("{:25s} {}\n".format("mcnp_input","../ETA.inp"))
            inp_file.write("{:25s} {}\n".format("mcnp_input_template","../ETA.inp"))
            inp_file.write("{:25s} {}\n".format("mcnp_tallies",S.tnum))
            inp_file.write("{:25s} {}\n".format("mcnp_force_point_source",S.pt_src))
            inp_file.write("{:25s} {}\n".format("mcnp_mix_tolerance",S.mix_tol))
            inp_file.write("\n")
            inp_file.write("{:25s} {}\n".format("anisn_library",S.lib.lower()))
            inp_file.write("\n")
            inp_file.write("{:25s} {}\n".format("denovo_pn_order",S.pn_order))
            inp_file.write("\n")

            # Run parallel if on cluster
            if cluster == True:
                cores=mp.cpu_count()
                inp_file.write('denovo_x_blocks    {}'.format(int(floor(sqrt(cores))))+'\n')
                inp_file.write('denovo_y_blocks    {}'.format(int(ceil(sqrt(cores))))+'\n')
            
            # Write the x mesh information
            xmesh=[]
            xints=[]
            xmesh.append(-(eta.r_o+2*S.ext))
            xmesh.append(-eta.r_o)
            xints.append(int(ceil((xmesh[-1]-xmesh[-2])/S.ext)))
            xmesh.append(-eta.r_toad)
            xints.append(int(ceil((xmesh[-1]-xmesh[-2])/S.eta_x)))
            xmesh.append(eta.r_toad)
            xints.append(int(ceil((xmesh[-1]-xmesh[-2])/S.foil_x)))
            xmesh.append(eta.r_o)
            xints.append(int(ceil((xmesh[-1]-xmesh[-2])/S.eta_x)))
            xmesh.append(eta.r_o+2*S.ext)
            xints.append(int(ceil((xmesh[-1]-xmesh[-2])/S.ext)))
                
            inp_file.write('{:25s} '.format("mesh_x")+ ' '.join('{}  '.format(k) for k in xmesh)+'\n')
            inp_file.write('{:25s} '.format("mesh_x_ints")+ ' '.join('{}  '.format(k) for k in xints)+'\n\n')
            
            # Write the y mesh information
            ymesh=[]
            yints=[]
            ymesh.append(-(eta.r_o+2*S.ext))
            ymesh.append(-eta.r_o)
            yints.append(int(ceil((ymesh[-1]-ymesh[-2])/S.ext)))
            ymesh.append(-eta.r_toad)
            yints.append(int(ceil((ymesh[-1]-ymesh[-2])/S.eta_y)))
            ymesh.append(eta.r_toad)
            yints.append(int(ceil((ymesh[-1]-ymesh[-2])/S.foil_y)))
            ymesh.append(eta.r_o)
            yints.append(int(ceil((ymesh[-1]-ymesh[-2])/S.eta_y)))
            ymesh.append(eta.r_o+2*S.ext)
            yints.append(int(ceil((ymesh[-1]-ymesh[-2])/S.ext)))
                
            inp_file.write('{:25s} '.format("mesh_y")+ ' '.join('{}  '.format(k) for k in ymesh)+'\n')
            inp_file.write('{:25s} '.format("mesh_y_ints")+ ' '.join('{}  '.format(k) for k in yints)+'\n\n')
            
            # Write the z mesh information
            zmesh=[]
            zints=[]
            zmesh.append(-S.ext)
            zmesh.append(eta.tcc_dist)
            zints.append(int(ceil((zmesh[-1]-zmesh[-2])/S.ext)))
            ind=next((i for i,item in enumerate(geom.surfaces) if item.c == "TOAD"), -1)
            zmesh.append(geom.surfaces[ind].vz-0.25)
            zints.append(int(ceil((zmesh[-1]-zmesh[-2])/S.eta_z)))
            zmesh.append(geom.surfaces[ind].vz+0.4503)
            zints.append(int(ceil((zmesh[-1]-zmesh[-2])/S.foil_z)))
            zmesh.append(eta.snout_dist+eta.t_m)
            zints.append(int(ceil((zmesh[-1]-zmesh[-2])/S.eta_z)))
            zmesh.append(eta.snout_dist+eta.t_m+S.ext*2)
            zints.append(int(ceil((zmesh[-1]-zmesh[-2])/S.ext)))
                
            inp_file.write('{:25s} '.format("mesh_z")+ ' '.join('{}  '.format(k) for k in zmesh)+'\n')
            inp_file.write('{:25s} '.format("mesh_z_ints")+ ' '.join('{}  '.format(k) for k in zints)+'\n')
            
        # Close the file
        inp_file.close()
    
    except IOError as e:
        module_logger.error("I/O error({0}): {1}".format(e.errno, e.strerror))   

    # Test that the file closed
    assert inp_file.closed==True, "File did not close properly."