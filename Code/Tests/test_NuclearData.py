#######################################################################################################
#
# Module : test_NuclearData.py
#
# Contains : Routines to test user inputs and defaults for NuclearData module
#
# Author : James Bevins
#
# Last Modified: 05May16
#
#######################################################################################################
""" NuclearData tests """

from pyne.dbgen.materials_library import make_matslib, make_elements

from NuclearData import Build_Matlib, Set_Density, Strip_Undesireables, Calc_Moderating_Ratio, Moderating_Ratio

from unittest import TestCase
import nose
import os
from nose.plugins.skip import SkipTest
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_in
    
#-------------------------------------------------------------------------------------------------------------#         
mat_path=os.getcwd()+"/Tests/files_test_Coeus/test_eta_materials_compendium.csv"

#-------------------------------------------------------------------------------------------------------------#     
def test_set_density1():
    mats=make_elements()
    mats=Set_Density(mats)
    assert_equal(mats['H'].density,0.0000899)
    assert_equal(mats['Cd'].density,8.65)
    assert_equal(mats['Pb'].density,11.34)  
    
def test_set_density():
    mats=make_matslib(mat_path)
    mats=Set_Density(mats)
    assert_equal(mats['H'].density,0.0000899)
    assert_equal(mats['Cd'].density,8.65)
    assert_equal(mats['Pb'].density,11.34)
    
def test_strip_undesirables1():
    mats=make_elements()
    l=len(mats)
    mats=Strip_Undesireables(mats,False,False,False)
    assert_equal(len(mats),l-1)
    mats=Strip_Undesireables(mats,True,False,False)
    assert_equal(len(mats),l-11)
    mats=Strip_Undesireables(mats,True,True,False)
    assert_equal(len(mats),l-14)
    mats=Strip_Undesireables(mats,True,True,True)
    assert_equal(len(mats),l-33)
    
def test_strip_undesirables2():
    mats=make_elements()
    l=len(mats)
    mats=Strip_Undesireables(mats,True,True,True)
    mats=Set_Density(mats)
    assert_equal(len(mats),l-33)
    
def test_build_matlib():
    mats=Build_Matlib(mat_path)
    assert_equal(mats['Cd'].density,8.65)
    assert_equal(mats['Pb'].density,11.34)
    assert_equal(len(mats),53)
    
def test_Moderating_Ratio():
    test=Moderating_Ratio("Pb",10.25,0.126)
    assert_equal(test.name,"Pb")
    assert_equal(test.mr_1MeV,10.25)
    assert_equal(test.mr_14MeV,0.126)
    
def test_Calc_Moderating_Ratio():
    mats=Build_Matlib(mat_path)
    test=Calc_Moderating_Ratio(mats)
    ind=next((i for i, item in enumerate(test) if item.name == 'Pb'), -1)
    assert_equal(test[ind].name,"Pb")
    assert_almost_equal(test[ind].mr_1MeV,7.8756626964889094)
    assert_almost_equal(test[ind].mr_14MeV,0.010808796195508964)