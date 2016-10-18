#######################################################################################################
#
# Module : test_Metaheuristics.py
#
# Contains : Routines to test user inputs and defaults for Metaheuristics module
#
# Author : James Bevins
#
# Last Modified: 17Aug16
#
#######################################################################################################
""" Metaheuristic tests """

import os

from Gnowee_Utilities import Gnowee_Settings, Parent 
from NuclearData import Calc_Moderating_Ratio, Build_Matlib
from MCNP_Utilities import MCNP_Geometry, MCNP_Settings
from Metaheuristics import Mat_Levy_Flights
from ETA_Utilities import ETA_Parameters

from unittest import TestCase
import nose
from nose.plugins.skip import SkipTest
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_in
    
#-------------------------------------------------------------------------------------------------------------#         
set_path=os.getcwd()+"/Tests/files_test_Coeus/test_gnowee_settings.csv"
mat_path=os.getcwd()+"/Tests/files_test_Coeus/test_eta_materials_compendium.csv"
#-------------------------------------------------------------------------------------------------------------#     
def test_Mat_Levy_Flights():
    geom=MCNP_Geometry()
    mcnp=MCNP_Settings()
    eta=ETA_Parameters()
    gs=Gnowee_Settings()
    gs.p=2
    mats=Build_Matlib(mat_path)
    geom.init_geom(eta, mats)
    pop=[]
    pop.append(Parent(0,eta,geom,gs,mcnp,mats,[''],0))
    pop.append(Parent(1,eta,geom,gs,mcnp,mats,[''],0))
    mr=Calc_Moderating_Ratio(mats)
    Mat_Levy_Flights(pop,mats,mr,gs,eta.fissile_mat)