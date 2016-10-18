#######################################################################################################
#
# Module : test_MCNP_Utilities.py
#
# Contains : Routines to test user inputs and defaults for MCNP_Utilities module
#
# Author : James Bevins
#
# Last Modified: 17Aug16
#
#######################################################################################################
""" MCNP-Utilities tests """

import os.path
import numpy as np

from MCNP_Utilities import MCNP_Settings, MCNP_Surface, MCNP_Cell, MCNP_Geometry, Print_MCNP_Input, Read_Tally_Output
from MCNP_Utilities import Read_MCNP_Output
from ETA_Utilities import ETA_Parameters
from NuclearData import Build_Matlib

from unittest import TestCase
import nose
from nose.plugins.skip import SkipTest
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_in
    
#-------------------------------------------------------------------------------------------------------------#         
source = [[1.0E-8,0.9], [1.0E-7,0.5], [1.0E-6,1.0E-1], [1.24E-2,10], 
            [2.345E+1,26.26]]

set_path=os.getcwd()+"/Tests/files_test_Coeus/test_mcnp_settings.csv"
src_path=os.getcwd()+"/Tests/files_test_Coeus/test_source.csv"
constraint_path=os.getcwd()+"/Tests/files_test_Coeus/test_eta_constraints.csv"
materials_library_path=os.getcwd()+"/Tests/files_test_Coeus/test_eta_materials_compendium.csv"

test_surf_repr='MCNP Surface(600, TRC, vx=1.0, vy=2.0, vz=3.0, hx=2.5, hy=23.6, hz=23.56, r1=3.4, r2=1.0, c=test)'

test_surf_str='700  px    2.00000  $test\n'
test_cell_repr='MCNP Cell:(1, mat=10, units=atom, density=0.0422, booleam geom=500 -501, n imp=1, p_imp=0, comment=)'

test_cell_str1='1  10  -4.22000e-02  500 -501  imp:n=1 imp:p=0  $\n'
test_cell_str2='1  10  -4.22000e-02  500 -501  imp:n=1 imp:p=0  $\n'
test_cell_str3='1  10            500 -501  imp:n=1 imp:p=0 $\n'
test_cell_str4='1  10  -4.22300e-02  (500 -501):(502 -503):(504 -505):(506 -507):(508\n     -509):(509 -510)  imp:n=1 imp:p=0  $\n'

mat_card="C name: Air (dry near sea level)\nC density = 0.0\nm?\n     6012 -1.2256e-04\n     6013 -1.4365e-06\n     7014 -7.5527e-01\n     8016 -2.3178e-01\n     18036 -3.8527e-05\n     18038 -7.6673e-06\n     18040 -1.2781e-02\n"

test_geom_str1='MCNP geometry instance properties:\nMCNP Surfaces:\n509  TRC    1.00000   2.00000   3.00000    2.50000  23.60000  23.56000   \n     3.40000    1.00000  $one\n\n504  px    2.00000  $two\n\n505  Py   -2.00000  $three\n\nMCNP Cells:\n1  11  -4.22000e-02  500 -501  imp:n=1 imp:p=0  $\n\n2  12  -4.22000e-02  500 -501  imp:n=1 imp:p=0  $\n\n3  13  -4.22000e-02  500 -501  imp:n=1 imp:p=0  $\n\nMCNP Materials:\nAir (dry near sea level)\nAl\n' 

base_geom='MCNP geometry instance properties:\nMCNP Surfaces:\n500  TRC    0.00000   0.00000  16.12650    0.00000   0.00000  14.35147   \n     0.00001    5.16119  $inner debris cover\n\n501  TRC    0.00000   0.00000  15.24000    0.00000   0.00000  15.23797   \n     0.00001    5.48000  $outer debris cover\n\n502  TRC    0.00000   0.00000  30.77797    0.00000   0.00000  10.57235   \n     5.26908    9.07119  $inner cone\n\n503  TRC    0.00000   0.00000  30.47797    0.00000   0.00000  10.87235   \n     5.48000    9.39000  $outer cone\n\n504  RCC    0.00000   0.00000  41.35032    0.00000   0.00000   9.78968   \n     8.89000  $inner cylinder\n\n505  RCC    0.00000   0.00000  41.35032    0.00000   0.00000   9.78968   \n     9.39000  $outer cylinder\n\n506  RCC    0.00000   0.00000  51.14000    0.00000   0.00000   1.00000   \n     9.39000  $cover\n\n507  RCC    0.00000   0.00000  52.14000    0.00000   0.00000   2.40000   \n     5.63400  $adapter\n\nMCNP Cells:\n1   1  -2.70000e+00  500 -501  imp:n=1 imp:p=0  $\n\n2   1  -2.70000e+00  502 -503  imp:n=1 imp:p=0  $\n\n3   1  -2.70000e+00  504 -505  imp:n=1 imp:p=0  $\n\n4   1  -2.70000e+00  -506  imp:n=1 imp:p=0  $\n\n5   1  -2.70000e+00  -507  imp:n=1 imp:p=0  $\n\nMCNP Materials:\nAl\nZr\nZn\nIn\nTa\nAu\nPb\nFe\n'
#-------------------------------------------------------------------------------------------------------------#     
def test_mcnp_surface1():
    surf=MCNP_Surface(500,"SO",r=2.534,comment="test")
    assert_equal(surf.name,500)
    assert_equal(surf.s_type,"SO")
    assert_equal(surf.r,2.534) 
    assert_equal(surf.c,"test") 
    
def test_mcnp_surface2():
    surf=MCNP_Surface(501,"CX",r=2.534,comment="test")
    assert_equal(surf.name,501)
    assert_equal(surf.s_type,"CX")
    assert_equal(surf.r,2.534) 
    assert_equal(surf.c,"test") 
    
def test_mcnp_surface3():
    surf=MCNP_Surface(502,"CY",r=2.534,comment="test")
    assert_equal(surf.name,502)
    assert_equal(surf.s_type,"CY")
    assert_equal(surf.r,2.534) 
    assert_equal(surf.c,"test") 
    
def test_mcnp_surface4():
    surf=MCNP_Surface(503,"CZ",r=2.534,comment="test")
    assert_equal(surf.name,503)
    assert_equal(surf.s_type,"CZ")
    assert_equal(surf.r,2.534) 
    assert_equal(surf.c,"test") 
    
def test_mcnp_surface5():
    surf=MCNP_Surface(504,"px",d=2.0,comment="test")
    assert_equal(surf.name,504)
    assert_equal(surf.s_type,"px")
    assert_equal(surf.d,2)
    assert_equal(surf.c,"test") 
    
def test_mcnp_surface6():
    surf=MCNP_Surface(505,"Py",d=-2.0,comment="test")
    assert_equal(surf.name,505)
    assert_equal(surf.s_type,"Py")
    assert_equal(surf.d,-2)
    assert_equal(surf.c,"test") 
    
def test_mcnp_surface7():
    surf=MCNP_Surface(506,"PZ",d=2.369,comment="test")
    assert_equal(surf.name,506)
    assert_equal(surf.s_type,"PZ")
    assert_equal(surf.d,2.369)
    assert_equal(surf.c,"test") 
    
def test_mcnp_surface8():
    surf=MCNP_Surface(507,"RCC",vx=1.0,vy=2.0,vz=3.0,hx=2.5,hy=23.6,hz=23.56,r=3.4,comment="test")
    assert_equal(surf.name,507)
    assert_equal(surf.s_type,"RCC")
    assert_equal(surf.vx,1.0) 
    assert_equal(surf.vy,2.0) 
    assert_equal(surf.vz,3.0) 
    assert_equal(surf.hx,2.5) 
    assert_equal(surf.hy,23.6) 
    assert_equal(surf.hz,23.56) 
    assert_equal(surf.r,3.4) 
    assert_equal(surf.c,"test") 
    
def test_mcnp_surface9():
    surf=MCNP_Surface(508,"RPP",x_min=1.0,x_max=2.0,y_min=3.0,y_max=4.0,z_min=5.0,z_max=6.0,comment="test")
    assert_equal(surf.name,508)
    assert_equal(surf.s_type,"RPP")
    assert_equal(surf.x_min,1.0) 
    assert_equal(surf.x_max,2.0) 
    assert_equal(surf.y_min,3.0) 
    assert_equal(surf.y_max,4.0) 
    assert_equal(surf.z_min,5.0) 
    assert_equal(surf.z_max,6.0) 
    assert_equal(surf.c,"test") 
    
def test_mcnp_surface10():
    surf=MCNP_Surface(509,"TRC",vx=1.0,vy=2.0,vz=3.0,hx=2.5,hy=23.6,hz=23.56,r1=3.4,r2=1.0,comment="test")
    assert_equal(surf.name,509)
    assert_equal(surf.s_type,"TRC")
    assert_equal(surf.vx,1.0) 
    assert_equal(surf.vy,2.0) 
    assert_equal(surf.vz,3.0) 
    assert_equal(surf.hx,2.5) 
    assert_equal(surf.hy,23.6) 
    assert_equal(surf.hz,23.56) 
    assert_equal(surf.r1,3.4) 
    assert_equal(surf.r2,1.0) 
    assert_equal(surf.c,"test") 
    
def test_mcnp_surface_repr():
    surf=MCNP_Surface(600,"TRC",vx=1.0,vy=2.0,vz=3.0,hx=2.5,hy=23.6,hz=23.56,r1=3.4,r2=1.0,comment="test")
    assert_equal(repr(surf),test_surf_repr)
    
def test_mcnp_surface_str():
    surf=MCNP_Surface(700,"px",d=2.0,comment="test")
    assert_equal(str(surf),test_surf_str)
    
def test_mcnp_cell1():
    cell=MCNP_Cell(1,10,"atom",0.0422, "500 -501", (1,0))
    assert_equal(cell.name,1)
    assert_equal(cell.m,10)
    assert_equal(cell.units,"atom")
    assert_equal(cell.d,0.0422)
    assert_equal(cell.geom,"500 -501")
    assert_equal(cell.imp,(1,0))
    
def test_mcnp_cell_repr():
    cell=MCNP_Cell(1,10,"atom",0.0422, "500 -501", (1,0))
    assert_equal(repr(cell),test_cell_repr)
    
def test_mcnp_cell_str1():
    cell=MCNP_Cell(1,10,"atom",0.0422, "500 -501", (1,0))
    assert_equal(str(cell),test_cell_str1)
    
def test_mcnp_cell_str2():
    cell=MCNP_Cell(1,10,"mass",0.0422, "500 -501", (1,0))
    assert_equal(str(cell),test_cell_str2)
    
def test_mcnp_cell_str3():
    cell=MCNP_Cell(1,10,"void",0.0422, "500 -501", (1,0))
    assert_equal(str(cell),test_cell_str3)
    
def test_mcnp_cell_str4():
    cell=MCNP_Cell(1,10,"mass",0.04223, "(500 -501):(502 -503):(504 -505):(506 -507):(508 -509):(509 -510)", (1,0))
    assert_equal(str(cell),test_cell_str4)
    
def test_mcnp_addMat():
    eta_params=ETA_Parameters()
    mat_lib=Build_Matlib(materials_library_path)
    geom=MCNP_Geometry()
    geom.add_matls(mat_lib, [eta_params.fill_mat,eta_params.struct_mat])
    assert_equal(geom.matls[0],eta_params.fill_mat)
    assert_equal(geom.matls[1],eta_params.struct_mat)
    assert_equal(mat_lib[geom.matls[0]].mcnp(),mat_card)
    geom.add_matls(mat_lib, [eta_params.fill_mat,eta_params.struct_mat])
    assert_equal(len(geom.matls),4)
    
def test_mcnp_addSurf():
    geom=MCNP_Geometry()
    geom.add_surf(MCNP_Surface(509,"TRC",vx=1.0,vy=2.0,vz=3.0,hx=2.5,hy=23.6,hz=23.56,r1=3.4,r2=1.0))
    assert_equal(geom.surfaces[0].name,509)
    assert_equal(geom.surfaces[0].s_type,"TRC")
    assert_equal(len(geom.surfaces),1)
    geom.add_surf([MCNP_Surface(504,"px",d=2.0),MCNP_Surface(505,"Py",d=-2.0)])
    assert_equal(geom.surfaces[1].name,504)
    assert_equal(geom.surfaces[1].s_type,"px")
    assert_equal(geom.surfaces[2].name,505)
    assert_equal(geom.surfaces[2].s_type,"Py")
    assert_equal(len(geom.surfaces),3)
    geom.add_surf(MCNP_Surface(504,"px",d=2.0))
    assert_equal(len(geom.surfaces),3)
    
def test_mcnp_addCell():
    geom=MCNP_Geometry()
    geom.add_cell(MCNP_Cell(1,11,"atom",0.0422, "500 -501", (1,0)))
    assert_equal(geom.cells[0].name,1)
    assert_equal(geom.cells[0].m,11)
    assert_equal(len(geom.cells),1)
    geom.add_cell([MCNP_Cell(2,12,"atom",0.0422, "500 -501", (1,0)),MCNP_Cell(3,13,"atom",0.0422, "500 -501", (1,0))])
    assert_equal(geom.cells[1].name,2)
    assert_equal(geom.cells[1].m,12)
    assert_equal(geom.cells[2].name,3)
    assert_equal(geom.cells[2].m,13)
    assert_equal(len(geom.cells),3)
    geom.add_cell(MCNP_Cell(1,10,"atom",0.0422, "500 -501", (1,0)))
    assert_equal(len(geom.cells),3)
    
def test_mcnp_geometry():
    geom=MCNP_Geometry()
    geom.add_cell([MCNP_Cell(1,11,"atom",0.0422, "500 -501", (1,0)), MCNP_Cell(2,12,"atom",0.0422, "500 -501",\
                  (1,0)),MCNP_Cell(3,13,"atom",0.0422, "500 -501", (1,0))])
    assert_equal(len(geom.cells),3)
    geom.add_surf([MCNP_Surface(509,"TRC",vx=1.0,vy=2.0,vz=3.0,hx=2.5,hy=23.6,hz=23.56,r1=3.4,r2=1.0,comment="one"),\
                    MCNP_Surface(504,"px",d=2.0,comment="two"),MCNP_Surface(505,"Py",d=-2.0,comment="three")])
    assert_equal(len(geom.surfaces),3)
    eta_params=ETA_Parameters()
    mat_lib=Build_Matlib(materials_library_path)
    geom.add_matls(mat_lib, [eta_params.fill_mat,eta_params.struct_mat])
    assert_equal(len(geom.matls),2)
    assert_equal(repr(geom),"MCNP geometry instance(There are 3 cells, 3 surfaces, and 2 materials used.)")
    assert_equal(str(geom),test_geom_str1)

def test_init_geometry():    
    geom=MCNP_Geometry()
    mat_lib=Build_Matlib(materials_library_path)
    eta_params=ETA_Parameters()
    geom.init_geom(eta_params, mat_lib)
    assert_equal(str(geom),base_geom)
    
def test_mcnp_settings_read_source():
    mcnp_set=MCNP_Settings()
    MCNP_Settings.read_source(mcnp_set,src_path)
    assert_equal(mcnp_set.source,source)
    
def test_Print_MCNP_Input():
    mcnp_set=MCNP_Settings()
    MCNP_Settings.read_source(mcnp_set,src_path)
    geom=MCNP_Geometry()
    mat_lib=Build_Matlib(materials_library_path)
    eta_params=ETA_Parameters()
    for i in range(0,25):
        Print_MCNP_Input(eta_params,geom,mcnp_set,mat_lib,i,adv_print=False)
    
def test_mcnp_settings_read_settings():
    mcnp_set=MCNP_Settings()
    MCNP_Settings.read_settings(mcnp_set,set_path)
    assert_equal(mcnp_set.phys,"Print\nMODE n\n")
    assert_equal(mcnp_set.nps,1E4)
    assert_equal(mcnp_set.tally,"")
    mcnp_set.set_tallies(3,2)
    assert_equal(mcnp_set.tally,'FC14 Fission Reaction Rate (Fissions per cm^3 per src particle)\nF14:n 3\nFM14  (-1 2 -6)     $Flux * atom density of material 2 * sigma f\nFC24 Uranium Flux Spectra (Number per cm^2 per src neutron)\nF24:n 3\n')
    
def test_Read_Tally_Output():
    tally=Read_Tally_Output(os.getcwd()+'/Tests/files_test_Coeus/ETA.out','24')  
    tgt=np.array([[1.0000E-08,0.00000E+00],[1.0000E-07,0.00000E+00],[1.0000E-06,0.00000E+00],[1.2400E-02,0.00000E+00],\
                 [2.3450E+01,8.27144E-06]])   
    np.testing.assert_equal(tally,tgt)
    
def test_Read_MCNP_Output():
    (t,r,w)=Read_MCNP_Output(os.getcwd()+'/Tests/files_test_Coeus/ETA.out','24','14')  
    tally=np.array([[1.0000E-08,0.00000E+00],[1.0000E-07,0.00000E+00],[1.0000E-06,0.00000E+00],[1.2400E-02,0.00000E+00],\
                 [2.3450E+01,8.27144E-06]])   
    rx=np.array([1.71263E-07, 0.3899])
    weight=4.36759E+03
    np.testing.assert_equal(t,tally)
    np.testing.assert_equal(r,rx)
    np.testing.assert_equal(w,weight)
    
def test_Read_MCNP_Output2():
    (t,r,w)=Read_MCNP_Output(os.path.abspath(os.getcwd())+'/Tests/files_test_Coeus/Population/0/tmp/ETA.out','24','14')  
    tally=np.array([[4.1399E-07,0.00000E+00],[1.1253E-06,0.00000E+00],[3.0590E-06,0.00000E+00],[1.0677E-05,0.00000E+00],[2.9023E-05,0.00000E+00],[1.0130E-04,0.00000E+00],[2.7536E-04,0.00000E+00],[5.8295E-04,2.09920E-10],[1.2341E-03,3.20839E-10],[3.3546E-03,3.21444E-10],[1.0333E-02,1.14966E-08],[2.1875E-02,4.11982E-08],[2.4788E-02,2.37702E-08],[3.4307E-02,6.19518E-08],[5.2475E-02,1.54382E-07],[1.1109E-01,7.82732E-07],[1.5764E-01,7.08239E-07],[2.4724E-01,1.31474E-06],[3.6883E-01,1.78305E-06],[5.5023E-01,2.10440E-06],[6.3928E-01,8.91525E-07],[7.4274E-01,9.33451E-07],[8.2085E-01,5.89541E-07],[9.6164E-01,9.74737E-07],[1.1080E+00,8.03300E-07],[1.4227E+00,1.39118E-06],[1.8268E+00,1.35067E-06],[2.3069E+00,9.98943E-07],[2.3852E+00,1.07852E-07],[3.0119E+00,7.99538E-07],[4.0657E+00,5.61116E-07],[4.7237E+00,2.38243E-07],[4.9659E+00,5.17845E-08],[6.3763E+00,2.43730E-07],[7.4082E+00,1.40065E-07],[8.1873E+00,7.72939E-08],[9.0484E+00,7.11370E-08],[1.0000E+01,9.62330E-08],[1.1052E+01,5.15155E-08],[1.2214E+01,1.05864E-07],[1.2523E+01,3.28099E-08],[1.3840E+01,5.13447E-07],[1.4191E+01,1.61314E-06],[1.4918E+01,3.43095E-08],[1.6905E+01,0.00000E+00],[1.9640E+01,0.00000E+00]])   
    rx=np.array([1.23990E-06,0.0131])
    weight=8.27728E+04
    np.testing.assert_equal(t,tally)
    np.testing.assert_equal(r,rx)
    np.testing.assert_equal(w,weight)