#######################################################################################################
#
# Module : test_Utilities.py
#
# Contains : Routines to test user inputs and defaults for Utilities module
#
# Author : James Bevins
#
# Last Modified: 17Aug16
#
#######################################################################################################
""" Utilities tests """

import os

import numpy as np

from Utilities import Cmd_Thread, Run_Transport_Threads, Run_Transport_PP, to_NormDiff, Uopt, Event, FuncThread 
from Utilities import FuncThreadWithReturn, LeastSquares, RelativeLeastSquares, Meta_Stats
from ETA_Utilities import ETA_Parameters
from NuclearData import Build_Matlib
from Gnowee_Utilities import Gnowee_Settings, Parent
from MCNP_Utilities import MCNP_Geometry, MCNP_Settings

from unittest import TestCase
import nose
from nose.plugins.skip import SkipTest
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_in
    
#-------------------------------------------------------------------------------------------------------------#         
test_mcnp='mcnp6 ../NSA_Proposal_ETA.inp NSA_Proposal_ETA.out'

#-------------------------------------------------------------------------------------------------------------#     
def test_cmd_thread():        
    test=Cmd_Thread(os.getcwd(),"mcnp6 i=ETA.inp")
    assert_equal(test.cwdir,os.getcwd())
    assert_equal(test.cmd,"mcnp6 i=ETA.inp")
    assert_equal(repr(test),"Thread instance({}, mcnp6 i=ETA.inp)".format(os.getcwd()))
    assert_equal(str(test),"\nThread Instance:\nThe current working directory is = {}\nThe cmd line input is = mcnp6 i=ETA.inp\n".format(os.getcwd()))  
    
def test_run_transport_pp():     
    # Initialize the colony
    eta_params=ETA_Parameters()
    gs=Gnowee_Settings()
    mcnp_set=MCNP_Settings(eta_params)
    mat_lib=Build_Matlib()
    base_eta=MCNP_Geometry()
    base_eta.init_geom(eta_params, mat_lib)
    colony=[]
    ids=[]
    for i in range(0,gs.p):
        colony.append(Parent(i, eta_params, base_eta, gs, mcnp_set, mat_lib, [''], i))
        ids.append(i)
    # Run the code
    Run_Transport_PP(ids,1,'mcnp6')
    
def test_run_transport_threads():     
    # Initialize the colony
    eta_params=ETA_Parameters()
    gs=Gnowee_Settings()
    mcnp_set=MCNP_Settings(eta_params)
    mat_lib=Build_Matlib()
    base_eta=MCNP_Geometry()
    base_eta.init_geom(eta_params, mat_lib)
    colony=[]
    ids=[]
    for i in range(0,gs.p):
        colony.append(Parent(i, eta_params, base_eta, gs, mcnp_set, mat_lib, [''], i))
        ids.append(i)
    # Run the code
    Run_Transport_Threads(ids,1,'mcnp6')
    
def test_to_normdiff(): 
    inp=np.array([[4.13990E-07,0.00000E+00],[1.12530E-06,0.00000E+00],[3.05900E-06,0.00000E+00],[1.06770E-05,0.00000E+00],[2.90230E-05,0.00000E+00],[1.01300E-04,0.00000E+00],[2.75360E-04,0.00000E+00],[5.82950E-04,2.09920E-10],[1.23410E-03,3.20839E-10],[3.35460E-03,3.21444E-10],[1.03330E-02,1.14966E-08],[2.18750E-02,4.11982E-08],[2.47880E-02,2.37702E-08],[3.43070E-02,6.19518E-08],[5.24750E-02,1.54382E-07],[1.11090E-01,7.82732E-07],[1.57640E-01,7.08239E-07],[2.47240E-01,1.31474E-06],[3.68830E-01,1.78305E-06],[5.50230E-01,2.10440E-06],[6.39280E-01,8.91525E-07],[7.42740E-01,9.33451E-07],[8.20850E-01,5.89541E-07],[9.61640E-01,9.74737E-07],[1.10800E+00,8.03300E-07],[1.42270E+00,1.39118E-06],[1.82680E+00,1.35067E-06],[2.30690E+00,9.98943E-07],[2.38520E+00,1.07852E-07],[3.01190E+00,7.99538E-07],[4.06570E+00,5.61116E-07],[4.72370E+00,2.38243E-07],[4.96590E+00,5.17845E-08],[6.37630E+00,2.43730E-07],[7.40820E+00,1.40065E-07],[8.18730E+00,7.72939E-08],[9.04840E+00,7.11370E-08],[1.00000E+01,9.62330E-08],[1.10520E+01,5.15155E-08],[1.22140E+01,1.05864E-07],[1.25230E+01,3.28099E-08],[1.38400E+01,5.13447E-07],[1.41910E+01,1.61314E-06],[1.49180E+01,3.43095E-08],[1.69050E+01,0.00000E+00],[1.96400E+01,0.00000E+00]])  
    tgt=np.array([0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0001,0.0007,0.0008,0.0031,0.0039,0.0169,0.0230,0.0504,0.0602,0.0877,0.0506,0.0547,0.0299,0.0444,0.0351,0.0620,0.0573,0.0578,0.0154,0.1007,0.0298,0.0170,0.0060,0.0205,0.0062,0.0052,0.0040,0.0034,0.0034,0.0033,0.0019,0.0117,0.0234,0.0624,0.0468,0.0005,])   
    tgt=np.round(tgt,decimals=2) 
    
    out=to_NormDiff(inp)
    out=np.round(out,decimals=2)
    np.testing.assert_equal(out,tgt)
    
def test_Uopt(): 
    candidate=np.array([8.0E+0,9.0E+02,8.0E+02,8.0E+02,7.0E+02,6.0E+02])
    design=np.array([9.0E+0,3.0E+02,4.0E+02,1.0E+03,7.0E+02,5.5E+02])
    expected_result=np.array([1.0E+0,6.0E+02,4.0E+02,2.0E+02,0.0,0.5E+02])
    result=Uopt(candidate,design)
    np.testing.assert_equal(result,np.sum(expected_result))
    
def test_LeastSquares(): 
    candidate=np.array([8.0E+0,9.0E+02,8.0E+02,8.0E+02,7.0E+02,6.0E+02])
    design=np.array([9.0E+0,3.0E+02,4.0E+02,1.0E+03,7.0E+02,5.5E+02])
    expected_result=np.array([1.00E+00,3.60E+05,1.60E+05,4.00E+04,0.00E+00,2.50E+03])
    result=LeastSquares(candidate,design)
    np.testing.assert_almost_equal(result,np.sum(expected_result))
    
def test_RelativeLeastSquares(): 
    candidate=np.array([2.10084034E-03,2.36344538E-01,2.10084034E-01,2.10084034E-01,1.83823529E-01,1.57563025E-01])
    design=np.array([3.04156810E-03,1.01385603E-01,1.35180804E-01,3.37952011E-01,2.36566408E-01,1.85873606E-01])
    expected_result=np.array([2.90958049E-04,1.79649905E-01,4.15036276E-02,4.83802998E-02,1.17591133E-02,4.31201072E-03])
    result=RelativeLeastSquares(candidate,design)
    np.testing.assert_almost_equal(result,np.sum(expected_result))
    
    candidate=np.array([0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00,2.849504E-07,1.151055E-06,1.713916E-06,4.650424E-06,2.685478E-05,2.012203E-04,7.592043E-04,6.746513E-04,2.667873E-03,4.510863E-03,1.784273E-02,2.211540E-02,5.474581E-02,6.644744E-02,7.973681E-02,4.231314E-02,4.725394E-02,2.657334E-02,3.946594E-02,3.201888E-02,5.873579E-02,5.147567E-02,4.711781E-02,1.230301E-02,8.064764E-02,2.258204E-02,1.372748E-02,4.738375E-03,1.544696E-02,5.021177E-03,4.028661E-03,3.287886E-03,3.179645E-03,3.600031E-03,4.670554E-03,2.959834E-03,1.885652E-02,3.708163E-02,9.843502E-02,7.375563E-02,9.867418E-04])
    design=np.array([4.066918E-13,2.825704E-11,7.962365E-10,7.373641E-09,2.404237E-08,1.185623E-07,4.099819E-07,1.204311E-06,4.002197E-06,2.135014E-05,1.563630E-04,6.279676E-04,5.722449E-04,2.080694E-03,2.725402E-03,1.179853E-02,1.589037E-02,3.674752E-02,4.712247E-02,7.091333E-02,4.365288E-02,4.831482E-02,2.674840E-02,4.052126E-02,3.298227E-02,6.017881E-02,5.755246E-02,6.032581E-02,1.847843E-02,1.261844E-01,4.444340E-02,2.998172E-02,1.113570E-02,3.622714E-02,9.127742E-03,5.905254E-03,3.046808E-03,1.812024E-03,1.438887E-03,1.387365E-03,6.852790E-04,6.874538E-03,1.843700E-02,5.523352E-02,6.085823E-02,9.803895E-03])
    expected_result=0.188372577544
    result=RelativeLeastSquares(candidate,design)
    np.testing.assert_almost_equal(result,expected_result)
    
def test_functhreadwithreturn():        
    t=[]
    result=[]
    
    candidate1=np.array([8.0E+0,9.0E+02,8.0E+02,8.0E+02,7.0E+02,6.0E+02])
    design1=np.array([9.0E+0,3.0E+02,4.0E+02,1.0E+03,7.0E+02,5.5E+02])
    expected_result1=np.array([1.0E+0,6.0E+02,4.0E+02,2.0E+02,0.0,0.5E+02])
    t.append(FuncThreadWithReturn(target=Uopt,args=(candidate1,design1,)))
    
    candidate2=np.array([5.0E+0,6.0E+02,5.0E+02,5.0E+02,4.0E+02,3.0E+02])
    design2=np.array([9.0E+0,3.0E+02,4.0E+02,1.0E+03,7.0E+02,5.5E+02])
    expected_result2=np.array([4.0E+0,3.0E+02,1.0E+02,5.0E+02,300.0,2.5E+02])
    t.append(FuncThreadWithReturn(target=Uopt,args=(candidate2,design2,)))
    
    t[0].start()
    t[1].start()
    result.append(t[0].join())
    result.append(t[1].join())
    
    np.testing.assert_equal(result[0],np.sum(expected_result1))
    np.testing.assert_equal(result[1],np.sum(expected_result2))
    
def test_functhread():        
    t=[]
    
    t.append(FuncThread(Event, 10,126,0.5678,1000.0,1))
    t.append(FuncThread(Event, 11,127,0.5679,1000.0,2))
    
    t[0].start()
    t[1].start()
    t[0].join()
    t[1].join()
    
    # Not best example because instances are lost to the ether...
    
def test_Event(): 
    test=Event(10,126,0.5678,1E6,5)
    assert_equal(test.g,10)
    assert_equal(test.e,126)
    assert_equal(test.f,0.5678)
    assert_equal(test.n,1E6)
    assert_equal(test.i,5)
    assert_equal(repr(test),"Event instance(10, 126, 0.5678, 1000000.0, 5)")
    assert_equal(str(test),"\nEvent Instance:\nGeneration # = 10\nNumber of Evaluations = 126\nFitness = 0.5678\nNPS = 1000000.0\nIdentity = 5\n")
    
def test_Meta_Stats():
    test=Meta_Stats()
    assert_equal(repr(test),'       0,0;               0,0;            0,0;             0,0;             0,0;            0,0;          0,0;       0,0;\n')
    assert_equal(str(test),'Mat_Levy_Flights  Cell_Levy_Flights  Mutate_Mats  Partial_Inversion  A_Operator  Study_Operator  Three_opt  Discard\n       0,0;               0,0;            0,0;             0,0;             0,0;            0,0;          0,0;       0,0;\n\n')
    test.update(mat_mut=(1,2),discard=(2,3))
    assert_equal(repr(test),'       0,0;               0,0;            1,2;             0,0;             0,0;            0,0;          0,0;       2,3;\n')
    test.write(header=True)