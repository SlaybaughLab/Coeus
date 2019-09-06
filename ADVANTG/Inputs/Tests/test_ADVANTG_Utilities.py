#######################################################################################################
#
# Module : test_ADVANTG_Utilities.py plus new test for User defined inputs
#
# Contains : Routines to test user inputs and defaults for ADVANTG_Utilities/ADVANTG_input module
#
# Author : James Bevins + SB
#
# Last Modified: 17Aug16/ August19
#
#######################################################################################################

""" ADVANTG-Utilities tests """

import os.path
import numpy as np

from ADVANTG_Utilities import ADVANTG_Settings, Print_ADVANTG_Input
from ETA_Utilities import ETA_Parameters
from NuclearData import Build_Matlib
from MCNP_Utilities import MCNP_Geometry, MCNP_Settings
from Gnowee_Utilities import Parent,Gnowee_Settings


from unittest import TestCase
import nose
from nose.plugins.skip import SkipTest
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_in
    
#-------------------------------------------------------------------------------------------------------------#         
constraint_path=os.getcwd()+"/Tests/files_test_Coeus/test_eta_constraints.csv"
mat_path=os.getcwd()+"/Tests/files_test_Coeus/test_eta_materials_compendium.csv"
set_path=os.getcwd()+"/Tests/files_test_Coeus/test_advantg_settings.csv"

test_settings_repr='ADVANTG Settings(dplus, cadis, mcnp silo, 24, True, 0.01, 1, 0.5, 0.5, 0.5, 0.25, 0.25, 0.05, 1.0)'

test_settings_str="\nADVANTG Program Settings:\n\
Multi-Group Library = dplus\n\
Solution Method = cadis\n\
Outputs = mcnp silo\n\
Adjoint Tally Number = 24\n\
Force Point Source = True\n\
Material Mix Tolerance = 0.01\n\
Scattering Order = 1\n\
ETA X Spacing Interval = 0.5\n\
ETA Y Spacing Interval = 0.5\n\
ETA Z Spacing Interval = 0.5\n\
Foil X Spacing Interval = 0.25\n\
Foil Y Spacing Interval = 0.25\n\
Foil Z Spacing Interval = 0.05\n\
External Spacing Interval = 1.0\n"
#-------------------------------------------------------------------------------------------------------------#     
def test_ADVANTG_setings():
    ADVANTG_set=ADVANTG_Settings()
    assert_equal(ADVANTG_set.lib,"dplus")
    assert_equal(ADVANTG_set.method,"cadis")
    assert_equal(ADVANTG_set.outputs,"mcnp silo")
    assert_equal(ADVANTG_set.tnum,24)
    assert_equal(ADVANTG_set.pt_src,"True")
    assert_equal(ADVANTG_set.mix_tol,0.01)
    assert_equal(ADVANTG_set.pn_order,1)
    assert_equal(ADVANTG_set.eta_x,0.5)
    assert_equal(ADVANTG_set.eta_y,0.5)
    assert_equal(ADVANTG_set.eta_z,0.5)
    assert_equal(ADVANTG_set.foil_x,0.25)
    assert_equal(ADVANTG_set.foil_y,0.25)
    assert_equal(ADVANTG_set.foil_z,0.05)
    assert_equal(ADVANTG_set.ext,1)
    
def test_ADVANTG_settings_repr():
    ADVANTG_set=ADVANTG_Settings()
    assert_equal(repr(ADVANTG_set),test_settings_repr)
    
def test_ADVANTG_settings_str():
    ADVANTG_set=ADVANTG_Settings()
    assert_equal(str(ADVANTG_set),test_settings_str)
    
def test_ADVANTG_settings_read_settings():
    ADVANTG_set=ADVANTG_Settings()
    ADVANTG_Settings.read_settings(ADVANTG_set,set_path)
    assert_equal(ADVANTG_set.lib,"DPLUS")
    assert_equal(ADVANTG_set.method,"CADIS")
    assert_equal(ADVANTG_set.outputs,"MCNP SILO")
    assert_equal(ADVANTG_set.tnum,14)
    assert_equal(ADVANTG_set.pt_src,"True")
    assert_equal(ADVANTG_set.mix_tol,0.1)
    assert_equal(ADVANTG_set.pn_order,1)
    assert_equal(ADVANTG_set.eta_x,2.0)
    assert_equal(ADVANTG_set.eta_y,2.0)
    assert_equal(ADVANTG_set.eta_z,2.0)
    assert_equal(ADVANTG_set.foil_x,1)
    assert_equal(ADVANTG_set.foil_y,1)
    assert_equal(ADVANTG_set.foil_z,0.2)
    assert_equal(ADVANTG_set.ext,2.5) 
    
def test_Print_ADVANTG_Input():
    adv_set=ADVANTG_Settings()
    geom=MCNP_Geometry()
    eta_params=ETA_Parameters()
    mat_lib=Build_Matlib(mat_path)
    mcnp=MCNP_Settings()
    gs=Gnowee_Settings()
    geom.init_geom(eta_params, mat_lib)
    p=Parent(0, eta_params, geom, gs, mcnp, mat_lib,[''], 0)
    for i in range(0,2):
        Print_ADVANTG_Input(eta_params,p.geom,adv_set,i) 
