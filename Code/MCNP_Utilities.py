#######################################################################################################
#
# Module : MCNP_Utilities.py
#
# Contains : Routines to build input decks for the MCNP radiation transport code.  
#
# Author : James Bevins
#
# Last Modified: 31Oct16
#
#######################################################################################################

import logging
module_logger = logging.getLogger('Coeus.MCNP_Utilities')

import bisect
import sys
import os

import Utilities as util
import numpy as np
import copy as cp

from math import ceil
from math import sin, cos, tan, atan, radians

class MCNP_Settings:

    ## Creates a object representing the settings for running the MCNP radiation trasport code.
    def __init__(self,physics="MODE n\n",nps=1E6, tally="", source=[]): 
        ## str The physics cards for run parameters
        self.phys=physics
        ## int The starting number of particles to run.  The number ran by the code will depend on generation and fitness.
        # [Default: 1E6]
        self.nps=nps
        ## array Stores the upper energy bin bounds and source strength for each bin
        # [default=[]]
        self.source=source  
        ## str The tallies for the problem.
        # [Default: ""]
        self.tally=tally
        
        
    def __repr__(self):
        return "MCNP Settings({}, {}, {}, {}, {})".format(self.phys, self.nps, self.tally, self.source)
    
    
    def __str__(self):
        header = ["\nMCNP Program Settings:"]
        header += ["Physics Cards = {}".format(self.phys)]
        header += ["Number of Source Particles = {}".format(self.nps)]
        header += ["Tallies = {}".format(self.tally)]
        header += ["\nSource spectra:"]
        header += ["Energy    Flux"]
        header ="\n".join(header)+"\n"
        s = header + "\n".join(["{0:<7}{1}".format(ebin, flux) for ebin, flux in self.source])
        return s
    
    ## Parses a MCNP settings csv input file. 
    #    The key word options are:
    #        Physics
    #        NPS
    def read_settings(self, filename):
        # Open file
        try: 
            self.f = open(filename, 'r') 
            
            # Initialize logic variables
            phys=False
            nps=False
            
            # Read the file line by line and store the values in the ETA_Params object
            for line in self.f:
                if phys==True:
                    if line.strip().lower()=='/':
                        phys=False
                    else:
                        self.phys=self.phys+line.strip()+'\n'
                if nps==True:
                    if line.strip().lower()=='/':
                        nps=False
                    else:
                        self.nps=float(line.strip())
                else:
                    for case in util.Switch(line.strip().lower()):
                        if case('Physics:'.lower()): 
                            phys=True
                            self.phys=""
                            break
                        if case('NPS:'.lower()): 
                            nps=True
                            break
                        if case(''): 
                            phys=False
                            nps=False
                            break
        
            # Close the file
            self.f.close()
        except IOError as e:
            module_logger.error( "I/O error({0}): {1}".format(e.errno, e.strerror)) 
            module_logger.error("File not found was: {0}".format(filename))  
       
        # Test that the file closed
        assert self.f.closed==True, "File ({}) did not close properly.".format(fname)
        
    ## Parses an source input csv file. 
    #    The first column contains the upper energy bin boundaries. 
    #    The second column contains the flux/fluence of the bin.
    def read_source(self, filename):

        # Open file
        try: 
            self.f = open(filename, 'r') 
            
            # Store the spectrum
            self.source=[]
            for line in self.f:
                split_list=line.split(',')
                self.source.append([float(split_list[0].strip()),float(split_list[1].strip())])
                
            # Close the file
            self.f.close()
        except IOError as e:
            module_logger.error("I/O error({0}): {1}".format(e.errno, e.strerror))
            module_logger.error("File not found was: {0}".format(filename))  
       
        # Test that the file closed
        assert self.f.closed==True, "File ({}) did not close properly.".format(fname)
    
    ## Sets the standard tallies to be used. 
    # @param cell int the cell for volume tallies
    # @param mat int the amterial number for reaction tallies
    def set_tallies(self, cell, mat):
        
        # Add Standard Tallies to User Defined Tallies
        self.tally+="FC14 Fission Reaction Rate (Fissions per cm^3 per src particle)\n"
        self.tally+="F14:n {}\n".format(cell)
        self.tally+="FM14  (-1 {} -6)     $Flux * atom density of material {} * sigma f\n".format(mat,mat)
        self.tally+="FC24 Uranium Flux Spectra (Number per cm^2 per src neutron)\n"
        self.tally+="F24:n {}\n".format(cell) 


class MCNP_Geometry:
    
    ## Creates the geometry for running the MCNP radiation trasport code.
    def __init__(self):  
        # [list of surface objects] Contains a list of all surface objects used in the design
        self.surfaces=[]
        # [list of cell objects] A list of all of the cell objects used in the design
        self.cells=[]
        # [list of material object keys] A list of the keys to material objects used in geometry.
        self.matls=[]
        
    def __repr__(self):
        return "MCNP geometry instance(There are {} cells, {} surfaces, and {} materials used.)".format(len(self.surfaces), len(self.cells), len(self.matls))
    
    
    def __str__(self):
        header = ["MCNP geometry instance properties:"]
        header += ["MCNP Surfaces:"]  
        for surf in self.surfaces:
            header += [str(surf)]
        header += ["MCNP Cells:"]  
        for cell in self.cells:
            header += [str(cell)]
        header += ["MCNP Materials:"]  
        for mat in self.matls:
            header += [str(mat)]
        header ="\n".join(header)+"\n"
        s = header 
        return s  
          
    ## Builds the inital surface list, cells dictionary, and materials list for the ETA geometry envelope
    # @param eta [ETA parameters object] An object that contains all of the constraints required to initialize the geometry
    # @param mats [dictionary of material objects] A materials library containing all relevant nulcear data required to run radiation transport codes.  
    #        Isotopic densities are in atoms/b-cm
    def init_geom(self, eta, mats):

        assert eta.r_f>0.0, 'The ETA face radius must be greater than zero'
        assert eta.theta>0.0, 'The ETA cone angle must be great than zero.'
        assert eta.r_o>=eta.r_f, 'The ETA outer radius must be greater than or equal to the face radius.'
        assert eta.tcc_dist>0.0, 'Distance to target chamber center must be positive and greater than zero.'
        assert eta.t_ds>0.0, 'The debris cover structural thickness must be greater than zero'
        assert eta.t_w>0.0, 'The ETA face structural thickness must be greater than zero'
        assert eta.t_c>0.0, 'The back cover thickness must be greater than zero.'  
        assert eta.t_m>0.0, 'The snout mounting thickness must be greater than zero.'        
        assert eta.snout_dist>eta.tcc_dist, 'Distance to the snout must be great thant the distance from the face to tcc.'

        # Add materials
        self.add_matls(mats, eta.struct_mat)
        if eta.ds_mat not in self.matls:
            self.add_matls(mats,eta.ds_mat)
        if eta.nas_mat not in self.matls:
            self.add_matls(mats,eta.nas_mat)
        for f in eta.nas_mat_f:
            if f not in self.matls:
                self.add_matls(mats,f)
        if eta.toad_mat not in self.matls:
            self.add_matls(mats,eta.toad_mat)
        for f in eta.toad_mat_f:
            if f not in self.matls:
                self.add_matls(mats,f)
        if eta.h_fill_mat not in self.matls:
            self.add_matls(mats,eta.h_fill_mat)
        if eta.holder_mat not in self.matls:
            self.add_matls(mats,eta.holder_mat)

        # Create the ETA shell surfaces
        h=eta.r_f*tan(radians(eta.theta))
        delh=(eta.r_f-eta.t_ds/sin(radians(eta.theta)))*tan(radians(eta.theta))
        self.add_surf(MCNP_Surface(500,"TRC",vx=0.0,vy=0.0,vz=eta.tcc_dist+h-delh,hx=0.0,hy=0.0,\
                      hz=delh,\
                      r1=0.00001,r2=eta.r_f-eta.t_ds/sin(radians(eta.theta)),comment="inner debris cover"))
        self.add_surf(MCNP_Surface(501,"TRC",vx=0.0,vy=0.0,vz=eta.tcc_dist,hx=0.0,hy=0.0,\
                      hz=h,\
                      r1=0.00001,r2=eta.r_f,comment="outer debris cover"))
        
        eta.tcc_dist+=h
        self.add_surf(MCNP_Surface(502,"TRC",vx=0.0,vy=0.0,vz=eta.tcc_dist+eta.t_ds,hx=0.0,hy=0.0,\
                      hz=(eta.r_o-eta.r_f)*tan(radians(eta.theta))-eta.t_ds,\
                      r1=eta.r_f+tan(radians(90.0-eta.theta))*eta.t_ds-eta.t_ds/sin(radians(eta.theta)),\
                                   r2=eta.r_o-eta.t_ds/sin(radians(eta.theta)),comment="inner cone"))
        self.add_surf(MCNP_Surface(503,"TRC",vx=0.0,vy=0.0,vz=eta.tcc_dist,hx=0.0,hy=0.0,\
                      hz=(eta.r_o-eta.r_f)*tan(radians(eta.theta)),\
                      r1=eta.r_f,r2=eta.r_o,comment="outer cone"))
        
        self.add_surf(MCNP_Surface(504,"RCC",vx=0.0,vy=0.0,vz=eta.tcc_dist+(eta.r_o-eta.r_f)*tan(radians(eta.theta)),\
                                   hx=0.0,hy=0.0,\
                                   hz=eta.snout_dist-eta.t_c-eta.tcc_dist-(eta.r_o-eta.r_f)*tan(radians(eta.theta)),\
                                   r=eta.r_o-eta.t_w,comment="inner cylinder"))
        self.add_surf(MCNP_Surface(505,"RCC",vx=0.0,vy=0.0,vz=eta.tcc_dist+(eta.r_o-eta.r_f)*tan(radians(eta.theta)),\
                                   hx=0.0,hy=0.0,\
                                   hz=eta.snout_dist-eta.t_c-eta.tcc_dist-(eta.r_o-eta.r_f)*tan(radians(eta.theta)),\
                                   r=eta.r_o,comment="outer cylinder"))
        
        self.add_surf(MCNP_Surface(506,"RCC",vx=0.0,vy=0.0,vz=eta.snout_dist-eta.t_c,\
                                   hx=0.0,hy=0.0,hz=eta.t_c,\
                                   r=eta.r_o,comment="cover"))
        
        self.add_surf(MCNP_Surface(507,"RCC",vx=0.0,vy=0.0,vz=eta.snout_dist,\
                                   hx=0.0,hy=0.0,hz=eta.t_m,\
                                   r=0.6*eta.r_o,comment="adapter"))
        
        # Create the ETA shell cells
        ind=next((i for i, item in enumerate(self.matls) if item == eta.ds_mat), -1)
        self.add_cell(MCNP_Cell(1,ind+1,"mass",mats[self.matls[ind]].density, "500 -501", (1,0)))
        self.add_cell(MCNP_Cell(2,ind+1,"mass",mats[self.matls[ind]].density, "502 -503", (1,0)))
        ind=next((i for i, item in enumerate(self.matls) if item == eta.struct_mat), -1)
        self.add_cell(MCNP_Cell(3,ind+1,"mass",mats[self.matls[ind]].density, "504 -505", (1,0)))
        self.add_cell(MCNP_Cell(4,ind+1,"mass",mats[self.matls[ind]].density, "-506", (1,0)))
        self.add_cell(MCNP_Cell(5,ind+1,"mass",mats[self.matls[ind]].density, "-507", (1,0)))
    
    ## Finishes the geometry by adding the filler and void cells, surfaces, and materials to the geometry
    # @param eta [ETA parameters object] An object that contains all of the constraints required to initialize the geometry
    # @param mats [dictionary of material objects] A materials library containing all relevant nulcear data required to run radiation transport codes.  
    #        Isotopic densities are in atoms/b-cm
    def fin_geom(self, eta, mats):
        assert eta.r_f>0.0, 'The ETA face radius must be greater than zero'
        assert eta.theta>0.0, 'The ETA cone angle must be great than zero.'
        assert eta.r_o>=eta.r_f, 'The ETA outer radius must be greater than or equal to the face radius.'
        assert eta.tcc_dist>0.0, 'Distance to target chamber center must be positive and greater than zero.'
        assert eta.t_w>0.0, 'The ETA face structural thickness must be greater than zero'
        assert eta.t_c>0.0, 'The back cover thickness must be greater than zero.' 
        assert eta.t_m>0.0, 'The snout mounting thickness must be greater than zero.'        
        assert eta.snout_dist>eta.tcc_dist, 'Distance to the snout must be great thant the distance from the face to tcc.'

        # Add materials
        if eta.fill_mat not in self.matls:
            self.add_matls(mats, eta.fill_mat)
            
        # Exclude vertical and foil cells in eta
        exclude=''
        for s in range(0,len(self.cells)-1):
            if self.cells[s].comment=="vert" or self.cells[s].comment=="foil":
                exclude+=" #{}".format(self.cells[s].name)  
        
        # Create the Inner Void cell in Debris Shield
        ind=next((i for i, item in enumerate(self.matls) if item == eta.fill_mat), -1)
        self.add_cell(MCNP_Cell(self.cells[-1].name+1,ind+1,"mass",mats[self.matls[ind]].density, "-{}".format(self.surfaces[0].name), (1,0)))
        
        # Create the Inner Void cell in ETA
        holder=next((i for i, item in enumerate(self.surfaces) if item.c.lower() == "Holder".lower()), -1)
        holder=self.surfaces[holder].name
        self.add_cell(MCNP_Cell(self.cells[-1].name+1,ind+1,"mass",mats[self.matls[ind]].density, \
                      "({0} -{1} {5} {4}):({0} -{2} {5} {4})".format(\
                       self.surfaces[-1].name, self.surfaces[2].name,\
                       self.surfaces[4].name,self.surfaces[5].name, exclude,holder), (1,0),comment="eta fill"))     
        
        # Create the outer colume surface
        self.add_surf(MCNP_Surface(self.surfaces[-1].name+1,"SO",r=eta.snout_dist*1.5,comment="kill radius"))
        
        # Create the Outer Volume cell
        self.add_cell(MCNP_Cell(self.cells[-1].name+1,0,"void",0.0,"{} {} {} {} {} -{}".format(self.surfaces[1].name,\
                       self.surfaces[3].name, self.surfaces[5].name, self.surfaces[6].name, self.surfaces[7].name, self.surfaces[-1].name), (1,0),comment="chamber fill"))
        
        # Create the Outer Kill cell
        self.add_cell(MCNP_Cell(self.cells[-1].name+1,0,"void",0.0,"{}".format(self.surfaces[-1].name), (0,0),comment="kill cell"))
        
    ## Adds new surface object to geometry surface list.
    # @param add A list of the surface objects to add
    def add_surf(self,adds):
        if isinstance(adds,list)==False:
            assert isinstance(adds, MCNP_Surface)==True, 'Surfaces in the MCNP geometry must be a MCNP_Surface instance.'
            if any(s.name==adds.name for s in self.surfaces): 
                module_logger.warning("WARNING: Surface {} already exists in this geometry.".format(adds.name))
            else:    
                self.surfaces.append(cp.deepcopy(adds))
            
        else:
            assert all(isinstance(x, MCNP_Surface) for x in adds)==True, 'Surfaces in the MCNP geometry must be a MCNP_Surface instance.'
            for y in adds:
                if any(s.name==y.name for s in self.surfaces): 
                    module_logger.warning("WARNING: Surface {} already exists in this geometry.".format(adds.name))
                else:    
                    self.surfaces.append(cp.deepcopy(y))
        
    ## Adds new cell object to geometry cells list.
    # @param adds A list of the cell objects to add
    def add_cell(self,adds):
        if isinstance(adds,list)==False:
            assert isinstance(adds, MCNP_Cell)==True, 'Cells in the MCNP geometry must be a MCNP_Surface instance.'
            if any(s.name==adds.name for s in self.cells): 
                module_logger.warning("WARNING: Surface {} already exists in this geometry.".format(adds.name))
            else:    
                self.cells.append(cp.deepcopy(adds))
            
        else:
            assert all(isinstance(x, MCNP_Cell) for x in adds)==True, 'Cells in the MCNP geometry must be a MCNP_Surface instance.'
            for y in adds:
                if any(s.name==y.name for s in self.cells): 
                    module_logger.warning("WARNING: Cell {} already exists in this geometry.".format(adds.name))
                else:    
                    self.cells.append(cp.deepcopy(y))
        
    ## Add materials to the matls list.  Checks for materials existing in the materials library.
    # @param mat_lib [dictionary of material objects] A materials library containing all relevant nulcear data required to run radiation transport codes
    # @param adds A list of the names of the materials to add to the matls list
    def add_matls(self,mat_lib,adds):
        # Add the material
        if isinstance(adds,list)==False:
            if adds in mat_lib.keys(): # and adds not in self.matls:
                self.matls.append(adds)
            else:
                module_logger.error("Material {} not found in the material library.".format(adds))

        else:
            for mat in adds:
                if mat in mat_lib.keys(): #and mat not in self.matls:
                    self.matls.append(mat)
                elif mat not in self.matls:
                    module_logger.error("Material {} not found in the material library.".format(mat))
    
    ## Builds the geometry object from an MCNP input file. Fairly specific to the current ETA design.
    # @param path String The path, including filename, to the MCNP output file to be read
    # @param mats [dictionary of material objects] A materials library containing all relevant nulcear data required to run radiation transport codes.  
    #        Isotopic densities are in atoms/b-cm
    # @return nps integer Number of particles for the MCNP run
    def read_geom(self, path, mats):
        assert isinstance(path, str)==True, 'Path must be of type str.'
      
        # Open input file 
        try:
            with open(path, "r") as f:
                module_logger.info("Importing ETA design at: {}".format(path))
                # Read the output file line by line
                for line in f:
                    # Find key word for start of Cells
                    if line =="c  Cell Cards  \n":
                        line=f.next()
                        line=f.next()
                        
                        while line!="\n":
                            while line.find('$')==-1:
                                line=line+f.next()
                            tmp=line.rstrip().split("$")
                            splt_lst=tmp[0].split()
                            splt_lst.append(tmp[1])
                            module_logger.debug("Cell list: {}".format(splt_lst))
                            if int(splt_lst[1])!=0:
                                if float(splt_lst[2])<0.0:
                                    self.add_cell(MCNP_Cell(int(splt_lst[0]),int(splt_lst[1]),"mass",\
                                       abs(float(splt_lst[2])), " ".join(splt_lst[3:-3]), (int(splt_lst[-3].split('=')[1]),\
                                       int(splt_lst[-2].split('=')[1])), comment=splt_lst[-1]))  
                                else:
                                    self.add_cell(MCNP_Cell(int(splt_lst[0]),int(splt_lst[1]),"atom",\
                                       abs(float(splt_lst[2])), " ".join(splt_lst[3:-3]), (int(splt_lst[-3].split('=')[1]),\
                                       int(splt_lst[-2].split('=')[1])), comment=splt_lst[-1]))  
                            else:
                                self.add_cell(MCNP_Cell(int(splt_lst[0]),int(splt_lst[1]),"void",0.0, \
                                       " ".join(splt_lst[2:-3]), (int(splt_lst[-3].split('=')[1]),\
                                       int(splt_lst[-2].split('=')[1])), comment=splt_lst[-1])) 
                            module_logger.debug("Cell after import: {}".format(str(self.cells[-1]))) 
                            
                            line=f.next() 
                            
                    elif line =="c  Surface Cards  \n":
                        line=f.next()
                        line=f.next()
                        
                        while line!="\n":
                            while line.find('$')==-1:
                                line=line+f.next()
                            tmp=line.rstrip().split("$")
                            splt_lst=tmp[0].split()
                            splt_lst.append(tmp[1])
                            module_logger.debug("Surface list: {}".format(splt_lst))
                            
                            # Determine surface type and add to surface list
                            if splt_lst[1].lower() == "so" or splt_lst[1].lower() == "cx" or splt_lst[1].lower() == "cy" or \
                                splt_lst[1].lower() == "cz":
                                self.add_surf(MCNP_Surface(int(splt_lst[0]), splt_lst[1], r=float(splt_lst[2]),\
                                                           comment=splt_lst[3]))    

                            # If specified surface typse is "PX", "PY", or "PZ", save appropriate optional arguments
                            elif splt_lst[1].lower() == "px" or splt_lst[1].lower() == "py" or splt_lst[1].lower() == "pz":
                                self.add_surf(MCNP_Surface(int(splt_lst[0]), splt_lst[1], d=float(splt_lst[2]),\
                                                           comment=splt_lst[3]))     

                            # If specified surface typse is RCC, save appropriate optional arguments
                            elif splt_lst[1].lower() == "rcc":
                                self.add_surf(MCNP_Surface(int(splt_lst[0]), splt_lst[1], vx=float(splt_lst[2]),\
                                                           vy=float(splt_lst[3]), vz=float(splt_lst[4]),\
                                                           hx=float(splt_lst[5]),hy=float(splt_lst[6]),\
                                                           hz=float(splt_lst[7]), r=float(splt_lst[8]),\
                                                           comment=splt_lst[9]))

                            # If specified surface typse is RPP, save appropriate optional arguments
                            elif splt_lst[1].lower() == "rpp":
                                self.add_surf(MCNP_Surface(int(splt_lst[0]), splt_lst[1], x_min=float(splt_lst[2]),\
                                                           x_max=float(splt_lst[3]), y_min=float(splt_lst[4]),\
                                                           y_max=float(splt_lst[5]), z_min=float(splt_lst[6]),\
                                                           z_max=float(splt_lst[7]),comment=splt_lst[8]))
                            # If specified surface typse is TRC, save appropriate optional arguments
                            elif splt_lst[1].lower() == "trc":
                                self.add_surf(MCNP_Surface(int(splt_lst[0]), splt_lst[1], vx=float(splt_lst[2]),\
                                                           vy=float(splt_lst[3]), vz=float(splt_lst[4]),\
                                                           hx=float(splt_lst[5]), hy=float(splt_lst[6]),\
                                                           hz=float(splt_lst[7]), r1=float(splt_lst[8]),\
                                                           r2=float(splt_lst[9]),comment=splt_lst[10]))
                            
                            module_logger.debug("Surface after import: {}".format(str(self.surfaces[-1]))) 
                            line=f.next()  
                                           
                    elif line=="c  Materials  \n":
                        while line!="c  Source":
                            line=f.next().rstrip()
                            if line.find('name:')!=-1:
                                self.add_matls(mats, line.split(':')[1].strip()) 
                                module_logger.debug("Imported material: {}".format(line.split(':')[1].strip()))
                                
                    elif line[0:3]=="NPS":
                        splt_lst=line.split()
                        nps=int(float(splt_lst[1]))
            # Close the file
            f.close()
    
        except IOError as e:
            module_logger.error("I/O error({0}): {1}".format(e.errno, e.strerror))
            module_logger.error("File not found was: {0}".format(path))

        # Test that the file closed
        assert f.closed==True, "File ({}) did not close properly.".format(path) 
    
        return nps

class MCNP_Surface:

    ## Creates a MCNP surface object.  Currently can handle SO, (PX,PY,PZ), (CX,CY,CZ), RCC, RPP, and TRC
    #    surfaces and macrobodies.  All others will throw an exception.  Attributes not used are specified as -1. 
    #    All atribute names follow those shown in Table 3.1 and Section 3.III.D in Volume II of the MCNP manual
    def __init__(self, name, s_type, r=-1, d=-1, x_min=-1, x_max=-1,  y_min=-1, y_max=-1, z_min=-1, z_max=-1,
                 vx=-1, vy=-1, vz=-1, hx=-1, hy=-1, hz=-1, r1=-1, r2=-1, comment=""):  
        ## Surface number
        self.name=name
        ## The type of MCNP surface.  Currently can specify "SO", ("PX","PY","PZ"), ("CX","CY","CZ"), "RCC", 
        # "RPP", and ("KX","KY","KY")
        self.s_type=s_type
        ## A radius in cm.  Used for the SO, CX, CY, CZ, and RCC surfaces
        # [Default: -1]
        self.r=r
        ## A location in cm.  Used for the PX, PY, and PZ surfaces
        # [Default: -1]
        self.d=d
        ## The minimum x location in cm.  Used for the RPP macrobody
        # [Default: -1]
        self.x_min=x_min
        ## The maximum x location in cm.  Used for the RPP macrobody
        # [Default: -1]
        self.x_max=x_max
        ## The minimum y location in cm.  Used for the RPP macrobody
        # [Default: -1]
        self.y_min=y_min
        ## The maximum y location in cm.  Used for the RPP macrobody
        # [Default: -1]
        self.y_max=y_max
        ## The minimum z location in cm.  Used for the RPP macrobody
        # [Default: -1]
        self.z_min=z_min
        ## The maximum z location in cm.  Used for the RPP macrobody
        # [Default: -1]
        self.z_max=z_max
        ## The x location of the center of the base in cm.  Used for the RCC and TRC macrobody
        # [Default: -1]
        self.vx=vx
        ## The y location of the center of the base in cm.  Used for the RCC and TRC macrobody
        # [Default: -1]
        self.vy=vy
        ## The z location of the center of the base in cm.  Used for the RCC and TRC macrobody
        # [Default: -1]
        self.vz=vz
        ## The change in x for the height vector in cm.  Used for the RCC and TRC macrobody
        # [Default: -1]
        self.hx=hx
        ## The change in y for the height vector in cm.  Used for the RCC and TRC macrobody
        # [Default: -1]
        self.hy=hy
        ## The change in z for the height vector in cm.  Used for the RCC and TRC macrobody
        # [Default: -1]
        self.hz=hz
        ## The radius of the lower cone base in cm.  Used for the TRC macrobody
        # [Default: -1]
        self.r1=r1
        ## The radius of the upper cone base in cm.  Used for the TRC macrobody
        # [Default: -1]
        self.r2=r2
        ## describing the surface feature.  Can be used to find the surface corresponding
        # to a particular geometric feature
        self.c=comment
        
        assert isinstance(self.name, int)==True, 'name must be of type int.'
        assert isinstance(self.s_type, str)==True, 's_types must be of type str.'
        assert isinstance(self.c, str)==True, 'comment must be of type str.'
        
        # If specified surface typse is "SO", check for appropriate optional arguments
        if self.s_type.lower() == "so" or self.s_type.lower() == "cx" or self.s_type.lower() == "cy" or \
           self.s_type.lower() == "cz":
            assert isinstance(self.r, float)==True, 'R must be of type float.'
            assert r>0, 'r must be greater than zero.'
            assert self.d==-1 and self.x_min==-1 and self.y_min==-1 and self.z_min==-1 and self.x_max==-1 \
                   and self.y_max==-1 and self.z_max==-1 and self.vx==-1 and self.vy==-1 and self.vz==-1 \
                   and self.hx==-1 and self.hy==-1 and self.hz==-1 and self.r1==-1 and self.r2==-1, \
                   'All non used attributes for surface type {} must not be specified.'.format(self.s_type)
        
        # If specified surface typse is "PX", "PY", or "PZ", check for appropriate optional arguments
        elif self.s_type.lower() == "px" or self.s_type.lower() == "py" or self.s_type.lower() == "pz":
            assert isinstance(self.d, float)==True, 'D must be of type float.'
            assert self.r==-1 and self.x_min==-1 and self.y_min==-1 and self.z_min==-1 and self.x_max==-1 \
                   and self.y_max==-1 and self.z_max==-1 and self.vx==-1 and self.vy==-1 and self.vz==-1 \
                   and self.hx==-1 and self.hy==-1 and self.hz==-1 and self.r1==-1 and self.r2==-1, \
                   'All non used attributes for surface type {} must not be specified.'.format(self.s_type)
       
        
        # If specified surface typse is RCC, check for appropriate optional arguments
        elif self.s_type.lower() == "rcc":
            assert isinstance(self.vx, float)==True, 'vx must be of type float.'
            assert isinstance(self.vy, float)==True, 'vy must be of type float.'
            assert isinstance(self.vz, float)==True, 'vz must be of type float.'
            assert isinstance(self.hx, float)==True, 'hx must be of type float.'
            assert isinstance(self.hy, float)==True, 'hy must be of type float.'
            assert isinstance(self.hz, float)==True, 'hz must be of type float.'
            assert isinstance(self.r, float)==True, 'R must be of type float.'
            assert self.r>0, 'r must be greater than zero.'
            assert self.d==-1 and self.x_min==-1 and self.y_min==-1 and self.z_min==-1 and self.x_max==-1 \
                   and self.y_max==-1 and self.z_max==-1 and self.r1==-1 and self.r2==-1, \
                   'All non used attributes for surface type {} must not be specified.'.format(self.s_type)
       
        
        # If specified surface typse is RPP, check for appropriate optional arguments
        elif self.s_type.lower() == "rpp":
            assert isinstance(self.x_min, float)==True, 'x_min must be of type float.'
            assert isinstance(self.x_max, float)==True, 'x_max must be of type float.'
            assert isinstance(self.y_min, float)==True, 'y_min must be of type float.'
            assert isinstance(self.y_max, float)==True, 'y_max must be of type float.'
            assert isinstance(self.z_min, float)==True, 'z_min must be of type float.'
            assert isinstance(self.z_max, float)==True, 'z_max must be of type float.'
            assert self.r==-1 and self.d==-1 and self.vx==-1 and self.vy==-1 and self.vz==-1 \
                   and self.hx==-1 and self.hy==-1 and self.hz==-1 and self.r1==-1 and self.r2==-1, \
                   'All non used attributes for surface type {} must not be specified.'.format(self.s_type)
       
        
        # If specified surface typse is TRC, check for appropriate optional arguments
        elif self.s_type.lower() == "trc":
            assert isinstance(self.vx, float)==True, 'vx must be of type float.'
            assert isinstance(self.vy, float)==True, 'vy must be of type float.'
            assert isinstance(self.vz, float)==True, 'vz must be of type float.'
            assert isinstance(self.hx, float)==True, 'hx must be of type float.'
            assert isinstance(self.hy, float)==True, 'hy must be of type float.'
            assert isinstance(self.hz, float)==True, 'hz must be of type float.'
            assert isinstance(self.r1, float)==True, 'r1 must be of type float.'
            assert isinstance(self.r2, float)==True, 'r2 must be of type float.'
            assert self.r1>0 and self.r2>0, 'r1 and r2 must be greater than zero.'
            assert self.r==-1 and self.d==-1 and self.x_min==-1 and self.y_min==-1 and self.z_min==-1 and self.x_max==-1 \
                   and self.y_max==-1 and self.z_max==-1, \
                   'All non used attributes for surface type {} must not be specified.'.format(self.s_type)
         
        # Catch all for non-covered surface types
        else:
            module_logger.error("An uknown surface type ({}) was specified.".format(s_type))
            sys.exit
            
            
    def __repr__(self):
        if self.s_type.lower() == "so" or self.s_type.lower() == "cx" or self.s_type.lower() == "cy" or \
           self.s_type.lower() == "cz":
            return "MCNP Surface({0}, {1}, r={2}, c={3})".format(self.name, self.s_type, self.r, self.c)
        
        # If specified surface typse is "PX", "PY", or "PZ", check for appropriate optional arguments
        elif self.s_type.lower() == "px" or self.s_type.lower() == "py" or self.s_type.lower() == "pz":
            return "MCNP Surface({0}, {1}, d={2}, c={3})".format(self.name, self.s_type, self.d, self.c)       
        
        # If specified surface typse is RCC, check for appropriate optional arguments
        elif self.s_type.lower() == "rcc":
            return "MCNP Surface({0}, {1}, vx={2}, vy={3}, vz={4}, hx={5}, hy={6}, hz={7}, r={8}, c={9})".format(self.name, \
                    self.s_type, self.vx, self.vy, self.vz, self.hx, self.hy, self.hz, self.r, self.c)
           
        # If specified surface typse is RPP, check for appropriate optional arguments
        elif self.s_type.lower() == "rpp":
            return "MCNP Surface({0}, {1}, x_min={2}, x_max={3}, y_min={4}, y_max={5}, z_min={6}, z_max={7}, c={8}"\
                    .format(self.name, self.s_type, self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, \
                            self.z_max, self.c)
        
        # If specified surface typse is TRC, check for appropriate optional arguments
        elif self.s_type.lower() == "trc":
            return "MCNP Surface({0}, {1}, vx={2}, vy={3}, vz={4}, hx={5}, hy={6}, hz={7}, r1={8}, r2={9}, c={10})"\
                   .format(self.name, self.s_type, self.vx, self.vy, self.vz, self.hx, self.hy, self.hz, \
                           self.r1, self.r2, self.c)
   
    
    def __str__(self):
        if self.s_type.lower() == "so" or self.s_type.lower() == "cx" or self.s_type.lower() == "cy" or \
           self.s_type.lower() == "cz":
            surf="{:3d}  {}  {:9.5f}  ${}\n".format(self.name, self.s_type, self.r, self.c)
       
        # If specified surface typse is "PX", "PY", or "PZ", check for appropriate optional arguments
        elif self.s_type.lower() == "px" or self.s_type.lower() == "py" or self.s_type.lower() == "pz":
            surf="{:3d}  {}  {:9.5f}  ${}\n".format(self.name, self.s_type, self.d, self.c)       
        
        # If specified surface typse is RCC, check for appropriate optional arguments
        elif self.s_type.lower() == "rcc":
            surf="{:3d}  {}  {:9.5f} {:9.5f} {:9.5f}  {:9.5f} {:9.5f} {:9.5f}  {:9.5f}  ${}\n"\
                    .format(self.name, self.s_type, self.vx, \
                    self.vy, self.vz, self.hx, self.hy, self.hz, self.r, self.c)
            
        # If specified surface typse is RPP, check for appropriate optional arguments
        elif self.s_type.lower() == "rpp":
            surf="{:3d}  {}  {:9.5f} {:9.5f} {:9.5f}  {:9.5f} {:9.5f} {:9.5f}  ${}\n".format(self.name, self.s_type, \
                    self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max, self.c)
        
        # If specified surface typse is TRC, check for appropriate optional arguments
        elif self.s_type.lower() == "trc":
            surf="{:3d}  {}  {:9.5f} {:9.5f} {:9.5f}  {:9.5f} {:9.5f} {:9.5f}  {:9.5f}  {:9.5f}  ${}\n"\
                    .format(self.name, self.s_type, self.vx, self.vy, self.vz, self.hx, self.hy, self.hz, \
                            self.r1, self.r2, self.c)
            
        # If the length approaches 80 columns, split over multiple lines
        if len(surf)>75:
            surf_list=surf.split(' ')
            surf=''
            tmp=''
            for i in surf_list:
                if len(tmp+i) < 75:
                    if tmp=='':
                        tmp=tmp+i
                    else:
                        tmp=tmp+' '+i
                else:
                    surf=surf+tmp+'\n    '  
                    tmp=' '+i
            surf+=tmp
            
        return surf
 

class MCNP_Cell:

    ##  Creates a MCNP cell object. 
    def __init__(self, name, mat, units, dens, geom, imp, comment=''):  
        
        assert isinstance(name, int)==True, 'name must be of type int.'
        assert isinstance(mat, int)==True, 'mat must be of type int.'
        assert isinstance(units, str)==True, 'units must be of type str.'
        assert isinstance(dens, float)==True, 'dens must be of type float.'
        assert isinstance(geom, str)==True, 'geom must be of type str.'
        assert isinstance(imp[0], int)==True and isinstance(imp[1], int)==True, 'all imp must be of type int.'
        assert isinstance(comment, str)==True, 'comment must be of type str.'
        
        ## int Cell number
        self.name=name
        ## int Material number.  0 for a void cell
        self.m=mat
        ## string Acceptable values are "atom", "mass", and "void".  Capitalization does not matter
        self.units=units
        ## float The density of the cell
        self.d=dens
        ## string Specification of the Boolean geometry of the cell. 
        self.geom=geom
        ## int tuple Specification of the importance of the regions for (neutron, photons)
        self.imp=imp
        ## string Comment describing the surface feature.  Can be used to find the surface corresponding
        # to a particular geometric feature
        # [Default='']
        self.comment=comment
               
    def __repr__(self):
        return "MCNP Cell:({0}, mat={1}, units={2}, density={3}, booleam geom={4}, n imp={5}, p_imp={6}, comment={7})"\
                .format(self.name, self.m, self.units, self.d, self.geom, self.imp[0], self.imp[1], self.comment)
    
    def __str__(self):
        if self.units.strip().lower()=="atom":
            cell="{}  {:2d}  {:.5e}  {}  imp:n={:1d} imp:p={:1d}  ${}\n".format(self.name, self.m, self.d, self.geom, self.imp[0], self.imp[1], self.comment)
        elif self.units.strip().lower()=="mass":
            cell="{}  {:2d}  {:.5e}  {}  imp:n={:1d} imp:p={:1d}  ${}\n".format(self.name, self.m, self.d*-1, self.geom, self.imp[0], self.imp[1], self.comment)
        elif self.units.strip().lower()=="void":
            cell="{}  {:2d}            {}  imp:n={:1d} imp:p={:1d} ${}\n".format(self.name, self.m, self.geom, self.imp[0], self.imp[1], self.comment)
        else:
            module_logger.error("Unknown value specified for density units type.  {} was specified.  Accepted values are atom, mass, and void.".format(self.units))
            
        # If the length approaches 80 columns, split over multiple lines
        if len(cell)>75:
            cell_list=cell.split(' ')
            cell=''
            tmp=''
            for i in cell_list:
                if len(tmp+i) < 75:
                    if tmp=='':
                        tmp=tmp+i
                    else:
                        tmp=tmp+' '+i
                else:
                    cell=cell+tmp+'\n    '  
                    tmp=' '+i
            cell+=tmp
            
        return cell

##  Print the generated MCNP input deck to file 
# @param eta [ETA parameters object] An object that contains all of the constraints required to initialize the geometry
# @param tallySpectrum [Numpy array] Contains the energy structure to be used for the tally.
# @param geom [MCNP_Geometry object] The geometry for running the MCNP radiation trasport code. Contains the surfaces, cells, and material information
# @param settings [MCNP_Settings object] An object representing the settings for running the MCNP radiation trasport code. Contains the source, physics, 
#        and tally information.
# @param mats [dictionary of material objects] A materials library containing all relevant nulcear data required to run radiation transport codes. 
# @param num int The current parent number being generated
# @param adv_print boolean (optional) An optional indicator to determine whether to print weight window and source bias information in the input file from 
#        ADVANTG outputs. 
def Print_MCNP_Input(eta,tallySpectrum, geom,settings,mats,num,adv_print=True):     
    path=os.path.abspath(os.path.join(os.path.abspath(os.getcwd()), os.pardir))+"/Results/Population/{}".format(num)
    # Delete previous input file if present
    if os.path.exists(path):
        if os.path.isfile("{}/ETA.inp".format(path)):
            os.remove("{}/ETA.inp".format(path))
    else:
        os.mkdir("{}".format(path))

    # Create and open input file 
    try:
        with open("{}/ETA.inp".format(path), "w") as inp_file:  

            # Print the header      
            inp_file.write("ETA design for Parent #{}\n".format(num))

            # Print Cell Cards
            inp_file.write("c ****************************************************************************\n")
            inp_file.write("c  Cell Cards  \n")
            inp_file.write("c ****************************************************************************\n")
            for c in geom.cells:
                inp_file.write("{}".format(str(c)))

            # Print Surface Cards
            inp_file.write("\n")
            inp_file.write("c ****************************************************************************\n")
            inp_file.write("c  Surface Cards  \n")
            inp_file.write("c ****************************************************************************\n")
            for s in geom.surfaces:  
                inp_file.write("{}".format(str(s)))
           
            # Print Data Cards
            inp_file.write("\n")
            inp_file.write("c ****************************************************************************\n")
            inp_file.write("c  Data Cards  \n")
            inp_file.write("c ****************************************************************************\n")
            inp_file.write("c  Physics  \n")
            inp_file.write("{}".format(settings.phys))
            inp_file.write("NPS {}\n".format(settings.nps))
            inp_file.write("RAND GEN=2 STRIDE=1529\n".format(settings.nps))

            # Print Material Cards
            inp_file.write("c ****************************************************************************\n")
            inp_file.write("c  Materials  \n")
            inp_file.write("c ****************************************************************************\n")
            i=1
            for key in geom.matls:
                str1=mats[key].mcnp().split('\n')
                str1[2]="m{}".format(i)
                str2=[]
                for s in range(0,len(str1)):
                    if str1[s][:9]!="     8018" and str1[s][:10]!="     73180" and str1[s][:10]!="     74180":
                                str2.append(str1[s])
                inp_file.write("{}".format('\n'.join("{}".format(i) for i in str2)))
                i+=1

            # Calculate cos(theta)
            theta=cos(atan(eta.r_f/eta.tcc_dist))-0.01
            if cos(atan(eta.r_o/(eta.tcc_dist+(eta.r_o-eta.r_f)*tan(radians(eta.theta)))))-0.01 < theta:
                theta=cos(atan(eta.r_o/(eta.tcc_dist+(eta.r_o-eta.r_f)*tan(radians(eta.theta)))))-0.01
                
            # Print Source Cards
            inp_file.write("c ****************************************************************************\n")
            inp_file.write("c  Source  \n")
            inp_file.write("c ****************************************************************************\n")
            inp_file.write("SDEF PAR=n ERG=d2 POS=0 0 0 VEC=0 0 1 \n")  
            inp_file.write("#   SI2           SP2      $ Source Spectrum\n")  
            inp_file.write("     1.00000E-012   0.00000E+00\n")
            for e,p in settings.source:
                inp_file.write("     {:6e}  {:6e}\n".format(e,p))
            
            # If ADVANTG files exist, read and print ADVANTG edits
            if os.path.exists(path+"/inp_edits.txt") \
               and os.path.exists(path+"/wwinp") and adv_print==True:
                adv=''
                try:
                    with open(path+"/inp_edits.txt", "r") as f:
            
                        # Read the output file line by line
                        for line in f:        
                            adv+=line
                    # Close the file
                    f.close()

                except IOError as e:
                    module_logger.error("I/O error({0}): {1}".format(e.errno, e.strerror))  
                    module_logger.error("File not found was: {0}".format(os.path.abspath(os.getcwd())+"/inp_edits.txt"))  

                # Test that the file closed
                assert f.closed==True, "File ({}) did not close properly.".format(os.path.abspath(os.getcwd())+ "/inp_edits.txt")
                
                # Print ADVANTG edits
                inp_file.write("c ****************************************************************************\n")
                inp_file.write("c Edits by ADVANTG\n") 
                inp_file.write("{}".format(adv))
                
            # If only one exists, output an error    
            elif adv_print==True and (os.path.exists("inp_edits.txt") and os.path.exists("wwinp")==False) or \
                 adv_print==True and (os.path.exists("inp_edits.txt")==False and os.path.exists("wwinp")):
                module_logger.error("ADVANTG input edits exist, but there is no corresponding wwinp file.")
                sys.exit
                
                
            # Print Tally Cards
            inp_file.write("c ****************************************************************************\n")
            inp_file.write("c  Tallies  \n")
            inp_file.write("{}".format(settings.tally))
            inp_file.write("E0  \n")
            for e,p in tallySpectrum:
                inp_file.write("      {:6e}\n".format(e))

        # Close the file
        inp_file.close()
    
    except IOError as e:
        module_logger.error("I/O error({0}): {1}".format(e.errno, e.strerror)) 
        module_logger.error("File not found was: {0}".format(os.path.abspath(os.getcwd())+"/ETA.inp"))

    # Test that the file closed
    assert inp_file.closed==True, "File did not close properly."
    
## Read the generated MCNP output and return the tally results
#  @param path String The path, including filename, to the MCNP output file to be read
#  @param tnum String The number of the tally to be read
#  @return tally array Array of tally results  
def Read_Tally_Output(path, tnum):
    
    assert isinstance(path, str)==True, 'Path must be of type str.'
    assert isinstance(tnum, str)==True, 'tnum must be of type str.'
    
    # Initialize the tally
    tally=[]
    t=False
    
    # Create and open input file 
    try:
        with open(path, "r") as f:
            
            # Read the output file line by line
            for line in f:
            
                # Find key word for start of flux array
                split_list=line.strip().split()
                if len(split_list)>=3:
                    if split_list[0].strip()=="1tally" and split_list[1].strip()==tnum.strip() and split_list[2].strip()=="nps":
                        t=True
                        # Advance 12 lines
                        for i in range(0,11):
                            split_list=f.next().split()
                
                if t==True:
                    # Exit at end of Tally
                    if split_list[0].strip()== "total":
                        t=False
                    else:
                        # Store Tally
                        tally.append([float(split_list[0].strip()),float(split_list[1].strip())])

        # Close the file
        f.close()
    
    except IOError as e:
        module_logger.error("I/O error({0}): {1}".format(e.errno, e.strerror))    
        module_logger.error("File not found was: {0}".format(path))

    # Test that the file closed
    assert f.closed==True, "File ({}) did not close properly.".format(path)
    
    return np.asarray(tally)   

## Read the generated MCNP output and return the tally results
    # @param path String The path, including filename, to the MCNP output file to be read
    # @param tnum String The number of the tally to be read.  Returns the entire binned tally.
    # @param rnum String The number of the tally to be read for the total reactions only. 
    # @return tally array Array of tally results for the tally specified by tnum [Ebins, tally, uncertainty]
    # @return rxs array Total number of reactions for the tally specified by rx_num [tally, uncertainty]
    # @return weight float The total weight of the system
def Read_MCNP_Output(path, tnum, rnum):

    assert isinstance(path, str)==True, 'Path must be of type str.'
    assert isinstance(tnum, str)==True, 'tnum must be of type str.'
    assert isinstance(rnum, str)==True, 'rnum must be of type str.'
    
    # Initialize the tally
    tally=[]
    t=False
    r=False
    w=False
    
    # Create and open input file 
    try:
        with open(path, "r") as f:
            
            # Read the output file line by line
            for line in f:
            
                # Find key word for start of flux array
                split_list=line.strip().split()
                if len(split_list)>=3:
                    if split_list[0].strip()=="1tally" and split_list[1].strip()==tnum.strip() and split_list[2].strip()=="nps":
                        t=True
                        # Advance 11 lines
                        for i in range(0,11):
                            split_list=f.next().split()
                    if split_list[0].strip()=="1tally" and split_list[1].strip()==rnum.strip() and split_list[2].strip()=="nps":
                        r=True
                        # Advance 12 lines
                        for i in range(0,12):
                            split_list=f.next().split()
                    if split_list[0].strip()=="cell" and split_list[1].strip()=="mat" and split_list[2].strip()=="density":
                        w=True
                        # Advance 2 lines
                        for i in range(0,2):
                            split_list=f.next().split()
                
                # Store data if keyword located
                if t==True:
                    # Exit at end of Tally
                    if split_list[0].strip()== "total":
                        t=False
                    else:
                        # Store Tally
                        tally.append([float(split_list[0].strip()),float(split_list[1].strip())])
                        
                if r==True:
                    # Store the total and exit at end of Tally
                    if split_list[0].strip()== "total":
                        rxs=[float(split_list[1].strip()),float(split_list[2].strip())]
                        r=False
                        
                if w==True:
                    # Store the total and exit at end of Tally
                    if len(split_list) > 0:
                        if split_list[0].strip()== "total":
                            weight=float(split_list[2].strip())
                            w=False

        # Close the file
        f.close()
   	 
   	 # Test that the file closed
    	assert f.closed==True, "File ({}) did not close properly.".format(path)    

    except IOError as e:
        module_logger.error("I/O error({0}): {1}".format(e.errno, e.strerror))
        module_logger.error("File not found was: {0}".format(path))

    return np.asarray(tally), np.asarray(rxs), weight
