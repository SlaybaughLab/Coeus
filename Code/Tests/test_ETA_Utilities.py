#######################################################################################################
#
# Module : test_ETA_Utilities.py
#
# Contains : Routines to test user inputs and defaults for ETA_Utilities module
#
# Author : James Bevins
#
# Last Modified: 17Aug16
#
#######################################################################################################
""" ETA-Utilities tests """

import os.path
import nose

import numpy as np

from ETA_Utilities import ETA_Parameters

from unittest import TestCase
from nose.plugins.skip import SkipTest
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_in
    
#-------------------------------------------------------------------------------------------------------------#  
print os.getcwd()
eta_path=os.getcwd()+"/Tests/files_test_Coeus/test_eta_constraints.csv"    
obj_path=os.getcwd()+"/Tests/files_test_Coeus/test_obj_spectrum.csv"

spectrum = np.array([[1.0E-8,1.0E-5], [1.0E-7,2.0E-5], [1.0E-6,1.0E-1], [1.24E-2,1.2546568188E+2], 
            [2.345E+1,2.53E+12]])

test_parameters_repr="ETA_Params(normalized differential, [], 1.59896054911e-06, 125.0, 5e+15, 15.24, 0.3, 0.5, 52.14, 1.0, 2.4, 5.48, 9.39, 70.22, Al, Al, Air (dry near sea level), Pb, 0.014, 2.69, Al, [0.1, 0.1, 0.1, 0.1, 0.01], 2.5, ['Zr', 'Zn', 'In', 'Al', 'Ta'], In, Al, [0.0254, 0.0127], 1.252, ['Au', 'Pb'], Al, Fe, 2.0, 3, 7)"

test_parameters_str="\nETA design constraints and objective function:\n\
Minimum number of fissions = 1.59896054911e-06 fissions\n\
Maximum ETA weight = 125.0 kg\n\
Source Neutrons in 4 pi = 5e+15 neutrons\n\
ETA distance from TCC = 15.24 cm\n\
Debris Shield thickness = 0.3 cm\n\
ETA structural thickness = 0.5 cm\n\
Snout distance from TCC = 52.14 cm\n\
ETA back cover thickness = 1.0 cm\n\
ETA to snout mount thickness = 2.4 cm\n\
ETA face radius = 5.48 cm\n\
ETA cylinder outer radius= 9.39 cm\n\
ETA cone opening angle = 70.22 degrees\n\
Debris Shield Material = Al\n\
ETA Structural Material = Al\n\
ETA Void Fill Material = Air (dry near sea level)\n\
Fissile Material = Pb\n\
NAS Thickness = 0.014 cm\n\
NAS Radius = 2.69 cm\n\
NAS Material = Al\n\
NAS Activation Foils = ['Zr', 'Zn', 'In', 'Al', 'Ta']\n\
NAS Activation Foil Thickness = [0.1, 0.1, 0.1, 0.1, 0.01] cm\n\
NAS Activation Foil Radius = 2.5 cm\n\
TOAD Follows Material = In\n\
TOAD Material = Al\n\
TOAD Activation Foils = ['Au', 'Pb']\n\
TOAD Activation Foil Thickness = [0.0254, 0.0127] cm\n\
TOAD Activation Foil Radius = 1.252 cm\n\
Holder Material = Al\n\
Holder Fill Material = Fe\n\
Holder wall thickness = 2.0\n\
Max vertical components = 3\n\
Max horizontal components = 7\n\
Objective function type = normalized differential\n\
\n\
Objective function spectra:\n\
Energy    Flux\n"
#-------------------------------------------------------------------------------------------------------------#     
def test_eta_parameters():
    eta_params=ETA_Parameters()
    assert_equal(eta_params.spectrum_type,'normalized differential')
    np.testing.assert_equal(eta_params.spectrum,np.array([]))
    assert_equal(eta_params.min_fiss,1.598960549105636e-06)
    assert_equal(eta_params.max_weight,125.0)
    assert_equal(eta_params.src,5E15)
    assert_equal(eta_params.tcc_dist,15.24)
    assert_equal(eta_params.t_ds,0.3)
    assert_equal(eta_params.t_w,0.5)
    assert_equal(eta_params.snout_dist,52.14)
    assert_equal(eta_params.t_c,1.0)
    assert_equal(eta_params.t_m,2.4)
    assert_equal(eta_params.r_f,5.48)
    assert_equal(eta_params.r_o,9.39) 
    assert_equal(eta_params.theta,70.22)
    assert_equal(eta_params.ds_mat,"Al") 
    assert_equal(eta_params.struct_mat,"Al") 
    assert_equal(eta_params.fill_mat,"Air (dry near sea level)")
    assert_equal(eta_params.fissile_mat,"Pb")
    assert_equal(eta_params.t_nas,0.014)
    assert_equal(eta_params.r_nas,2.69) 
    assert_equal(eta_params.nas_mat,"Al")
    assert_equal(eta_params.t_nas_f,[0.1, 0.1, 0.1, 0.1, 0.01])
    assert_equal(eta_params.r_nas_f,2.5) 
    assert_equal(eta_params.nas_mat_f,['Zr', 'Zn', 'In', 'Al', 'Ta'])
    assert_equal(eta_params.toad_loc,"In")
    assert_equal(eta_params.toad_mat,"Al")
    assert_equal(eta_params.t_toad,[0.0254, 0.0127])
    assert_equal(eta_params.r_toad,1.252) 
    assert_equal(eta_params.toad_mat_f,["Au","Pb"])
    assert_equal(eta_params.holder_mat,"Al")
    assert_equal(eta_params.h_fill_mat,"Fe")
    assert_equal(eta_params.t_h,2.0)
    assert_equal(eta_params.max_vert,3) 
    assert_equal(eta_params.max_horiz,7)  
    
def test_eta_parameters_repr():
    eta_params=ETA_Parameters()
    assert_equal(repr(eta_params),test_parameters_repr)
    
def test_eta__parameters_str():
    eta_params=ETA_Parameters()
    assert_equal(str(eta_params),test_parameters_str)
    
def test_eta_parameters_read_obj():
    eta_params=ETA_Parameters()
    ETA_Parameters.read_obj(eta_params,obj_path)
    assert_equal(eta_params.spectrum_type,'normalized differential')
    np.testing.assert_equal(eta_params.spectrum,spectrum)
    
def test_eta_parameters_read_constraints():
    eta_params=ETA_Parameters()
    ETA_Parameters.read_constraints(eta_params,eta_path)
    assert_almost_equal(eta_params.min_fiss,1.10524266E-5)
    assert_equal(eta_params.max_weight,100.263)
    assert_equal(eta_params.src,1E15)
    assert_equal(eta_params.tcc_dist,22)
    assert_equal(eta_params.t_ds,0.1)
    assert_equal(eta_params.t_w,0.5)
    assert_equal(eta_params.snout_dist,59.61)
    assert_equal(eta_params.t_c,0.8)
    assert_equal(eta_params.t_m,1.68)
    assert_equal(eta_params.r_f,3.64)
    assert_equal(eta_params.r_o,7.59) 
    assert_equal(eta_params.theta,7.51)
    assert_equal(eta_params.ds_mat,"Al") 
    assert_equal(eta_params.struct_mat,"Aluminum") 
    assert_equal(eta_params.fill_mat,"C") 
    assert_equal(eta_params.fissile_mat,"U")
    assert_equal(eta_params.t_nas,0.012)
    assert_equal(eta_params.r_nas,2.5) 
    assert_equal(eta_params.nas_mat,"Fe")
    assert_equal(eta_params.t_nas_f,[0.2, 0.1, 0.4])
    assert_equal(eta_params.r_nas_f,2.26) 
    assert_equal(eta_params.nas_mat_f,['Be', 'B', 'C'])
    assert_equal(eta_params.toad_loc,"B")
    assert_equal(eta_params.toad_mat,"Fe")
    assert_equal(eta_params.t_toad,[0.01, 0.02])
    assert_equal(eta_params.r_toad,1.2) 
    assert_equal(eta_params.toad_mat_f,["U","Au"])
    assert_equal(eta_params.holder_mat,"Al")
    assert_equal(eta_params.h_fill_mat,"Pb")
    assert_equal(eta_params.t_h,3.0)
    assert_equal(eta_params.max_vert,2) 
    assert_equal(eta_params.max_horiz,5) 