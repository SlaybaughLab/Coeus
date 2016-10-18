#######################################################################################################
#
# Module : test_Gnowee_Utilities.py
#
# Contains : Routines to test user inputs and defaults for Gnowee_Utilities module
#
# Author : James Bevins
#
# Last Modified: 17Aug16
#
#######################################################################################################
""" Gnowee-Utilities tests """

import numpy as np

import os.path

from ETA_Utilities import ETA_Parameters
from NuclearData import Build_Matlib
from Gnowee_Utilities import Gnowee_Settings, Parent, Calc_Fitness, Timeline, Pop_Update, Calc_Fitness, Rejection_Bounds
from MCNP_Utilities import MCNP_Geometry, MCNP_Settings, Print_MCNP_Input
from NuclearData import Build_Matlib, Calc_Moderating_Ratio
from Utilities import Run_Transport_Threads

from unittest import TestCase
import nose
from nose.plugins.skip import SkipTest
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_in
    
#-------------------------------------------------------------------------------------------------------------#         
test_repr='Gnowee Settings(25, lhc, 0.25, 0.2, 0.4, 10000, 100000, 1e-06, 200, 0.01, 0.01, 1.5, 1.0, 1, 10.0)'

test_str='\nGnowee Optimization Settings:\nPopulation Size = 25\nInitial sampling method = lhc\nDiscovery Fraction = 0.2\nElite fraction = 0.2\nLevy fraction = 0.4\nMaximum number of genrations = 10000\nMaximum number of function evaluations = 100000\nStall convergence tolerance = 1e-06\nStall iteration limit = 200\nOptimal fitness = 0.01\nOptimal convergence tolerance = 0.01\nLevy exponent = 1.5\nLevy scale unit = 1.0\nLevy independent variables = 1\nStep size scaling factor = 10.0\n'

set_path=os.getcwd()+"/Tests/files_test_Coeus/test_gnowee_settings.csv"
mcnp_path=os.getcwd()+"/Tests/files_test_Coeus/test_mcnp_settings.csv"
source_path=os.getcwd()+"/Tests/files_test_Coeus/test_source.csv"
eta_path=os.getcwd()+"/Tests/files_test_Coeus/test_eta_constraints.csv"
mat_path=os.getcwd()+"/Tests/files_test_Coeus/test_eta_materials_compendium.csv"
obj_path=os.getcwd()+"/Tests/files_test_Coeus/test_obj_spectrum.csv"

#-------------------------------------------------------------------------------------------------------------#     
def test_gnowee_settings():
    g_set=Gnowee_Settings()
    assert_equal(g_set.p,25)
    assert_equal(g_set.s,'lhc')
    assert_equal(g_set.fd,0.25)
    assert_equal(g_set.fe,0.20)
    assert_equal(g_set.fl,0.40)
    assert_equal(g_set.gm,10000)
    assert_equal(g_set.em,100000)
    assert_equal(g_set.ct,1e-6)
    assert_equal(g_set.sl,200)
    assert_equal(g_set.of,0.01)
    assert_equal(g_set.ot,1e-2)
    assert_equal(g_set.a,1.5)
    assert_equal(g_set.g,1.0)
    assert_equal(g_set.n,1)
    assert_equal(g_set.sf,10.0)
    assert_equal(repr(g_set),test_repr)
    assert_equal(str(g_set),test_str)
    
def test_gnowee_settings_read_settings():
    g_set=Gnowee_Settings()
    Gnowee_Settings.read_settings(g_set,set_path)
    assert_equal(g_set.p,18)
    assert_equal(g_set.s,'random')
    assert_equal(g_set.fd,0.59)
    assert_equal(g_set.fe,0.26)
    assert_equal(g_set.fl,1.0)
    assert_equal(g_set.gm,2156)
    assert_equal(g_set.em,2620)
    assert_equal(g_set.ct,0.26)
    assert_equal(g_set.sl,4665)
    assert_equal(g_set.of,236.56)
    assert_equal(g_set.ot,0.2)
    assert_equal(g_set.a,0.5)
    assert_equal(g_set.g,1.2)
    assert_equal(g_set.n,10)
    assert_equal(g_set.sf,120.3)
    
def test_parent():
    geom=MCNP_Geometry()
    mcnp=MCNP_Settings()
    eta=ETA_Parameters()
    gs=Gnowee_Settings()
    mats=Build_Matlib(mat_path)
    geom.init_geom(eta, mats)
    p=Parent(0,eta,geom,gs,mcnp,mats,[''],gs.p-1)
    assert_equal(p.ident,0)
    assert_equal(p.fit,1E15)
    assert_equal(p.rset.nps,1E6)
    assert_equal(p.fixed_mats,8)
    
def test_timeline():
    pop=[]
    geom=MCNP_Geometry()
    mcnp=MCNP_Settings()
    eta=ETA_Parameters()
    gs=Gnowee_Settings()
    mats=Build_Matlib(mat_path)
    geom.init_geom(eta, mats)
    pop.append(Parent(0,eta,geom,gs,mcnp,mats,[''],gs.p-1,fitness=15.0))    
    pop.append(Parent(1,eta,geom,gs,mcnp,mats,[''],gs.p-1,fitness=10.0))
    
    history=Timeline()
    history.update(pop, 1, 20)
    assert_equal(history.tline[-1].g,0)
    assert_equal(history.tline[-1].e,20)
    assert_equal(history.tline[-1].f,10.0)
    assert_equal(history.tline[-1].n,1E6)
    assert_equal(history.tline[-1].i,1)
    assert_equal(repr(history),"0  20  1.000000e+01  1.00e+06  1\n")
    assert_equal(str(history),"\nGenerations  Evaluations  Fitness  NPS  ID\n    0     20   10.0000  1.00e+06  1")
    
def test_pop_update():
    geom=MCNP_Geometry()
    mcnp=MCNP_Settings(mcnp_path)
    mcnp.read_source(source_path)
    eta=ETA_Parameters(eta_path)
    eta.read_obj(obj_path)
    gs=Gnowee_Settings()
    mats=Build_Matlib(mat_path)
    geom.init_geom(eta, mats)
    geom.fin_geom(eta, mats)
    old=[]
    old.append(Parent(0,eta,geom,gs,mcnp,mats,[''],0))
    new=[]
    new.append(Parent(0,eta,geom,gs,mcnp,mats,[''],0))
    new[0].fit=1E14
    Pop_Update(old,new,mcnp.nps,eta,mats,Run_Transport_Threads)
    assert_equal(old[0].ident,0)
    assert_equal(old[0].fit,1E14)
    new[0].fit=0.45
    Pop_Update(old,new,mcnp.nps)
    assert_equal(old[0].fit,0.45)
    assert_equal(old[0].rset.nps,10*mcnp.nps)
    new[0].fit=0.105
    Pop_Update(old,new,mcnp.nps)
    assert_equal(old[0].fit,0.105)
    assert_equal(old[0].rset.nps,100*mcnp.nps)
    new[0].fit=0.437254e-01
    Pop_Update(old,new,mcnp.nps)
    assert_equal(old[0].fit,0.0437254)
    assert_equal(old[0].rset.nps,100*mcnp.nps)
    
def test_Rejection_Bounds():
    cur=np.array([0.25,0.35,0.45])
    new=np.array([-1.25,2.35,3.45])
    lb=np.array([0,0,0])
    ub=np.array([3,3,3])
    new=Rejection_Bounds(cur,new,lb,ub,change_count=0)
    np.testing.assert_equal(new,np.array([0.25,2.35,0.45]))
    cur=[0.25,0.35,0.45]
    new=[-1.25,2.35,3.45]
    lb=[0,0,0]
    ub=[3,3,3]
    new=Rejection_Bounds(cur,new,lb,ub,change_count=0)
    assert_equal(new,[0.25,2.35,0.45])