#######################################################################################################
#
# Module : ETA_Utilities.py
#
# Contains : Functions and methods required to set ETA contraints and objective functions.  
#            This module will:
#            1) Read ETA contraints from ETA constraints file
#            2) Store ETA related design parameters in ETA class
#            3) Evaluate fitness of objective function using U-optimility
#
# Author : James Bevins
#
# Last Modified: 23April17
#
#######################################################################################################

import logging
module_logger = logging.getLogger('Coeus.ETA_Utilities')

import numpy as np
import Utilities as util

from math import pi

#-------------------------------------------------------------------------------------------------------------#      
class ETA_Parameters:
     
    ## Creates an ETA object that stores the ETA design parameters, constraints, and objective function
    def __init__(self, min_fiss=5E8, max_weight=125.0,src=5E15,\
                 tcc_dist=15.24, debris_shield_thickness= 0.3, wall_thickness=0.5, snout_dist=52.14, cover_thickness=1.0, mount_thickness=2.4, 
                 face_radius=5.48, eta_or=9.39, cone_angle=70.22,\
                 debris_shield_mat="Al", struct_mat="Al", fill_mat="Air (dry near sea level)",fissile_mat='Pb',\
                 nas_th=0.014,nas_r=2.69,nas_mat='Al',\
                 nas_foil_th=[0.1, 0.1, 0.1, 0.1, 0.01],nas_foil_r=2.5,nas_foil_mat=['Zr', 'Zn', 'In', 'Al', 'Ta'],toad_loc='In',\
                 toad_mat='Al',toad_foil_th=[0.0254,0.0127],toad_foil_r=1.252,toad_foil_mat=['Au','Pb'],\
                 holder_mat='Al', holder_fill_mat='Fe', holder_thickness=2.0, \
                 max_vert=3, max_horiz=7):
        assert len(nas_foil_th)==len(nas_foil_mat), 'The number of thicknesses and materials for the NAS must be equal.'
        assert len(toad_foil_th)==len(toad_foil_mat), 'The number of thicknesses and materials for the TOAD must be equal.'
        
        ## float The maximum weight for the ETA assembly entered in kg
        # [default=125.0 [kg]]
        self.max_weight=max_weight
        ## float NIF source neutron in 4pi
        # [default=5E15]
        self.src=src
        
        ## float The distance from the front face of the ETA assembly from target chamber center (TCC) entered in cm
        # [default=15.24 [cm]]
        self.tcc_dist=tcc_dist
        ## float The thickness of the debris cover and conical section of the ETA assembly entered in cm
        # [default=0.3 [cm]]
        self.t_ds=debris_shield_thickness
        ## float The thickness of the walls of the ETA assembly entered in cm
        # [default=0.25 [cm]]
        self.t_w=wall_thickness
        ## float The distance to where the snout mounts in cm. Measured from target chamber center (TCC)
        # [default=52.14 [cm]]
        self.snout_dist=snout_dist
        ## float The thickness of the back cover for the ETA in cm
        # [default=1.0 [cm]]
        self.t_c=cover_thickness
        ## float The thickness of the mount connecting the ETA Nose Cone Assembly to the snout in cm
        # [default=2.4 [cm]]
        self.t_m=mount_thickness
        ## float The opening radius of the ETA assembly entered in cm. Measured from ETA centerline
        # [default=5.48 [cm]]
        self.r_f=face_radius
        ## float The maximum outer radius of the ETA assembly structure entered in cm. Measured from ETA centerline
        # [default=9.39 [cm]]
        self.r_o=eta_or
        ## float The opening angle of the cone measured in degrees. Measured from ETA face plane.
        # [default=8.89 [cm]]
        self.theta=cone_angle
        
        ## string The material used for the conical debris cover.  Must be a naturally occuring element or specified in the materials compendium.
        # [default="Al"]
        self.ds_mat=debris_shield_mat
        ## string The material used for the ETA structure.  Must be a naturally occuring element or specified in the materials compendium.
        # [default="Al"]
        self.struct_mat=struct_mat
        ## string The material used to fill all nonspecified volume in the in the ETA.  Must be a naturally occuring element or specified 
        # in the materials compendium.
        # [default="Air (dry near sea level)"]
        self.fill_mat=fill_mat
        ## string The fissile material used in the ETA.  
        # Must be a naturally occuring element or specified in the materials compendium.
        # [default='U']
        self.fissile_mat=fissile_mat
        
        ## float The thickness of the neutron activation spectrometer in cm
        # [default= 0.14 cm]
        self.t_nas=nas_th
        ## float The radius of the neutron activation spectrometer in cm.
        # [default= 2.69 [cm]]
        self.r_nas=nas_r
        ## string The material used for the NAS structure.  Must be a naturally occuring element or specified in the materials compendium.
        # [default="Al"]
        self.nas_mat=nas_mat
        
        ## list of floats The thickness of the neutron activation spectrometer foils in cm
        # [default= [0.1, 0.1, 0.1, 0.1, 0.01] cm]
        self.t_nas_f=nas_foil_th
        ## float The radius of the neutron activation spectrometer foils entered in cm.
        # [default= 2.5 [cm]]
        self.r_nas_f=nas_foil_r
        ## list of strings The material used for the neutron activation spectrometer foils in the in the ETA.  
        # Must be a naturally occuring element or specified in the materials compendium.
        # [default=['Zr', 'Zn', 'In', 'Al', 'Ta']]
        self.nas_mat_f=nas_foil_mat
        ## str String indicating the material that the TOAD assembly follows after in the stackup.
        # [default=['In']
        self.toad_loc=toad_loc
        
        ## string The material used for the TOAD structure.  Must be a naturally occuring element or specified in the materials compendium.
        # [default="Al"]
        self.toad_mat=toad_mat
        ## list of floats The thickness of the TOAD foils in cm
        # [default= [,0.0127] cm]
        self.t_toad=toad_foil_th
        ## float The radius of foils entered in cm.
        # [default= 1.252 [cm]]
        self.r_toad=toad_foil_r
        ## list of strings The material used for the neutron activation spectrometer foils in the in the ETA.  
        # Must be a naturally occuring element or specified in the materials compendium.
        # [default=['Au','U']]
        self.toad_mat_f=toad_foil_mat
        
        ## string The material used for the holder structure.  Must be a naturally occuring element or specified in the materials compendium.
        # [default="Al"]
        self.holder_mat=holder_mat
        ## string The material used to fill the space int he holder structrure.  Must be a naturally occuring element or specified in the materials compendium.
        # [default="Fe"]
        self.h_fill_mat=holder_fill_mat
        ## float The thickness of the holder for the NAS insertion assemblyentered in cm.
        # [default= 2 [cm]]
        self.t_h=holder_thickness
        
        ## int The maximum number of vertical macrobodies or components in the ETA geometry. This cn be reduced to increase run spead or
        # increased to obtain a better result.
        # [default=3]
        self.max_vert=max_vert
        ## int The maximum number of horizontal macrobodies or components in the ETA geometry. This cn be reduced to increase run spead
        # or increased to obtain a better result.
        # [default=7]
        self.max_horiz=max_horiz
        
        # Calculate the number of fissions based on the source strength and volume
        ind=next((i for i in enumerate(self.toad_mat_f) if i == self.fissile_mat), -1)
        ##  float The minimum number of fissions in fissile foil
        # [default=5e8]
        self.min_fiss=min_fiss/(self.src*self.r_toad**2*pi*self.t_toad[ind])
    
    def __repr__(self):
        return "ETA_Params({0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}, {11}, {12}, {13}, {14}, {15}, {16}, {17}, {18}, {19}, {20}, {21}, {22}, {23}, {24}, {25}, {26}, {27}, {28}, {29}, {30}, {31})".format(\
                self.min_fiss, self.max_weight, self.src,\
                self.tcc_dist, self.t_ds, self.t_w, self.snout_dist, self.t_c, self.t_m, self.r_f, self.r_o, self.theta,\
                self.ds_mat, self.struct_mat, self.fill_mat,self.fissile_mat,\
                self.t_nas, self.r_nas, self.nas_mat,\
                self.t_nas_f, self.r_nas_f, self.nas_mat_f, self.toad_loc,\
                self.toad_mat, self.t_toad, self.r_toad, self.toad_mat_f,\
                self.holder_mat, self.h_fill_mat, self.t_h,\
                self.max_vert, self.max_horiz)
    
    def __str__(self):
        header = ["\nETA design constraints and objective function:"]
        header += ["Minimum number of fissions = {} fissions".format(self.min_fiss)]
        header += ["Maximum ETA weight = {} kg".format(self.max_weight)]
        header += ["Source Neutrons in 4 pi = {} neutrons".format(self.src)]
        
        header += ["ETA distance from TCC = {} cm".format(self.tcc_dist)]
        header += ["Debris Shield thickness = {} cm".format(self.t_ds)]
        header += ["ETA structural thickness = {} cm".format(self.t_w)]
        header += ["Snout distance from TCC = {} cm".format(self.snout_dist)]
        header += ["ETA back cover thickness = {} cm".format(self.t_c)]
        header += ["ETA to snout mount thickness = {} cm".format(self.t_m)]
        header += ["ETA face radius = {} cm".format(self.r_f)]
        header += ["ETA cylinder outer radius= {} cm".format(self.r_o)]
        header += ["ETA cone opening angle = {} degrees".format(self.theta)]
        
        header += ["Debris Shield Material = {}".format(self.ds_mat)]
        header += ["ETA Structural Material = {}".format(self.struct_mat)]
        header += ["ETA Void Fill Material = {}".format(self.fill_mat)]
        header += ["Fissile Material = {}".format(self.fissile_mat)]
        
        header += ["NAS Thickness = {} cm".format(self.t_nas)]
        header += ["NAS Radius = {} cm".format(self.r_nas)]
        header += ["NAS Material = {}".format(self.nas_mat)]
        
        header += ["NAS Activation Foils = {}".format(self.nas_mat_f)]
        header += ["NAS Activation Foil Thickness = {} cm".format(self.t_nas_f)]
        header += ["NAS Activation Foil Radius = {} cm".format(self.r_nas_f)]
        header += ["TOAD Follows Material = {}".format(self.toad_loc)]
        
        header += ["TOAD Material = {}".format(self.toad_mat)]
        header += ["TOAD Activation Foils = {}".format(self.toad_mat_f)]
        header += ["TOAD Activation Foil Thickness = {} cm".format(self.t_toad)]
        header += ["TOAD Activation Foil Radius = {} cm".format(self.r_toad)]
        
        header += ["Holder Material = {}".format(self.holder_mat)]
        header += ["Holder Fill Material = {}".format(self.h_fill_mat)]
        header += ["Holder wall thickness = {}".format(self.t_h)]
        
        header += ["Max vertical components = {}".format(self.max_vert)]
        header += ["Max horizontal components = {}".format(self.max_horiz)]
        
        return "\n".join(header)+"\n"

    ## Parses an ETA constraints csv file. 
    #     The key word options are:
    #         Minimum Fissions
    #         ETA Max Weight
    #         Source Strength
            
    #         TCC to ETA Distance
    #         Debris Shield Thickness
    #         ETA Wall Thickness
    #         Snout Distance
    #         ETA Back Cover Thickness
    #         ETA to Snout Mount Thickness
    #         ETA Face Radius
    #         ETA Cone Outer Radius
    #         ETA Cone Opening Angle
            
    #         Debris Shield Material
    #         ETA Structural Material
    #         ETA Void Fill Material
    #         Fissile Foil

    #         NAS Thickness
    #         NAS Radius
    #         NAS Material
            
    #         NAS Activation Foils
    #         NAS Activation Foil Thickness
    #         NAS Activation Foil Radius
    #         TOAD Follows Material
            
    #         TOAD Material
    #         TOAD Activation Foils
    #         TOAD Activation Foil Thickness
    #         TOAD Activation Foil Radius
            
    #         Holder Material
    #         Holder Fill Material
    #         Holder Wall Thickness

    #         Max Vertical Components
    #         Max Horizontal Components
    def read_constraints(self, filename):

        # Open file
        try: 
            self.f = open(filename, 'r') 
            
            # Read the file line by line and store the values in the ETA_Params object
            for line in self.f:
                split_list=line.split(',')
                for case in util.Switch(split_list[0].strip().lower()):
                    if case('Minimum Fissions'.lower()): 
                        self.min_fiss=float(split_list[1].strip())
                        break
                    if case('ETA Max Weight'.lower()): 
                        self.max_weight=float(split_list[1].strip())
                        break
                    if case('Source Strength'.lower()): 
                        self.src=float(split_list[1].strip())
                        break
                        
                    if case('TCC to ETA Distance'.lower()): 
                        self.tcc_dist=float(split_list[1].strip())
                        break
                    if case('Debris Shield Thickness'.lower()): 
                        self.t_ds=float(split_list[1].strip())
                        break
                    if case('ETA Wall Thickness'.lower()): 
                        self.t_w=float(split_list[1].strip())
                        break
                    if case('Snout Distance'.lower()): 
                        self.snout_dist=float(split_list[1].strip())
                        break
                    if case('ETA Back Cover Thickness'.lower()): 
                        self.t_c=float(split_list[1].strip())
                        break
                    if case('ETA to Snout Mount Thickness'.lower()): 
                        self.t_m=float(split_list[1].strip())
                        break
                    if case('ETA Face Radius'.lower()): 
                        self.r_f=float(split_list[1].strip())
                        break
                    if case('ETA Cone Outer Radius'.lower()): 
                        self.r_o=float(split_list[1].strip())
                        break
                    if case('ETA Cone Opening Angle'.lower()): 
                        self.theta=float(split_list[1].strip())
                        break
                        
                    if case('Debris Shield Material'.lower()): 
                        self.ds_mat=split_list[1].strip()
                        break
                    if case('ETA Structural Material'.lower()): 
                        self.struct_mat=split_list[1].strip()
                        break
                    if case('ETA Void Fill Material'.lower()): 
                        self.fill_mat=split_list[1].strip()
                        break
                    if case('Fissile Foil'.lower()): 
                        self.fissile_mat=split_list[1].strip()
                        break
                        
                    if case('NAS Thickness'.lower()): 
                        self.t_nas=float(split_list[1].strip())
                        break
                    if case('NAS Radius'.lower()): 
                        self.r_nas=float(split_list[1].strip())
                        break
                    if case('NAS Material'.lower()): 
                        self.nas_mat=split_list[1].strip()
                        break
                        
                    if case('NAS Activation Foils'.lower()): 
                        self.nas_mat_f=[]
                        for i in range(1, len(split_list)):
                            self.nas_mat_f.append(split_list[i].strip())
                        break
                    if case('NAS Activation Foil Thickness'.lower()): 
                        self.t_nas_f=[]
                        for i in range(1, len(split_list)):
                            self.t_nas_f.append(float(split_list[i].strip()))
                        break
                    if case('NAS Activation Foil Radius'.lower()): 
                        self.r_nas_f=float(split_list[1].strip())
                        break
                    if case('TOAD Follows Material'.lower()): 
                        self.toad_loc=split_list[1].strip()
                        break
                    
                    if case('TOAD Material'.lower()): 
                        self.toad_mat=split_list[1].strip()
                        break    
                    if case('TOAD Activation Foils'.lower()): 
                        self.toad_mat_f=[]
                        for i in range(1, len(split_list)):
                            self.toad_mat_f.append(split_list[i].strip())
                        break
                    if case('TOAD Activation Foil Thickness'.lower()): 
                        self.t_toad=[]
                        for i in range(1, len(split_list)):
                            self.t_toad.append(float(split_list[i].strip()))
                        break
                    if case('TOAD Activation Foil Radius'.lower()): 
                        self.r_toad=float(split_list[1].strip())
                        break
                    
                    if case('Holder Material'.lower()): 
                        self.holder_mat=split_list[1].strip()
                        break    
                    if case('Holder Fill Material'.lower()): 
                        self.h_fill_mat=split_list[1].strip()
                        break
                    if case('Holder Wall Thickness'.lower()): 
                        self.t_h=float(split_list[1].strip())
                        break
                        
                    if case('Max Vertical Components'.lower()): 
                        self.max_vert=int(split_list[1].strip())
                        break
                    if case('Max Horizontal Components'.lower()): 
                        self.max_horiz=int(split_list[1].strip())
                        break
                        
                    if case('/'): 
                        break
                    if case(): # default, could also just omit condition or 'if True'
                        module_logger.warning("\n A user input ({}) was found in the ETA constraints file that does not match the allowed input types. Minimum Fissions, ETA Max Weight,Source Strength\
                        TCC to ETA Distance, Debris Shield Thickness, ETA Wall Thickness, Snout Distance, ETA Back Cover Thickness, ETA to Snout Mount Thickness, \
                        ETA Face Radius, ETA Cone Inner Radius, ETA Cone Opening Angle, \
                        Debris Shield Material, ETA Structural Material, ETA Void Fill Material, Fissile Mat,\
                        NAS Thickness, NAS Radius, NAS Material\
                        NAS Activation Foils, NAS Activation Foil Thickness, NAS Activation Foil Radius, TOAD Follows Material,\
                        TOAD Material, TOAD Activation Foils, TOAD Activation Foil Thickness, TOAD Activation Foil Radius, \
                        Holder Material,Holder Fill Material,Holder Wall Thickness,\
                        Max Vertical Components, Max Horizontal Components".format(split_list[0].strip()))
        
            # Close the file
            self.f.close()
            
        except IOError as e:
            module_logger.error("I/O error({0}): {1}".format(e.errno, e.strerror)) 
            module_logger.error("File not found was: {0}".format(filename)) 
       
        # Test that the file closed
        assert self.f.closed==True, "File did not close properly."
        
        # Calculate the number of fissions based on the source strength and volume
        for i in range(0,len(self.toad_mat_f)):
            if self.toad_mat_f[i]==self.fissile_mat:
                ind=i
        self.min_fiss=self.min_fiss/(self.src*self.r_toad**2*pi*self.t_toad[ind])
               
     