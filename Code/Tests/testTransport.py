"""!
@file testTransport.py
@package CoeusTesting

@defgroup testTransport testTransport

@brief Routines to test testTransport module.

@author James Bevins

@date 18Aug19
"""

import os

from nose.tools import assert_equal
from UserInputs import UserInputs
from Transport import Transport

#-----------------------------------------------------------------------------#
# Assumed inputs
PATH = os.getcwd() +'\\Tests\\files_test_Coeus\\'
INPUTFNAME = PATH + 'test_user_inputs.txt'

#-----------------------------------------------------------------------------#
def test_Transport_MCNP():
    """
    Test the default object creation for MCNP input.
    """
    # Read user input
    inputs = UserInputs(coeusInputPath=INPUTFNAME)
    inputs.read_inputs()

    # Read transport input (MCNP Example)
    trans = Transport(PATH+inputs.transInput)
    assert_equal(trans.transPath, PATH+"test_mcnp.inp")
    assert_equal(trans.code, "mcnp6")
    assert_equal(trans.sampVars, {'h1': ' 1:5> decimal',
                                  'h2': ' 1:5> decimal',
                                  'h3': ' 1:5> decimal',
                                  'h4': ' 1:5> decimal',
                                  'h5': ' 1:5> decimal',
                                  'h6': ' 1:5> decimal',
                                  'd1': ' 1:3> integer',
                                  'r1': ' 0.5:3> decimal',
                                  'mat1': ' 1,2,3,4,5,6,7,8,9,10,11,12> material'})
    assert_equal(trans.corrVars, {'dens1': ' -2.7,-7.8,-6.5,-8.9,-7.3,-2.7,-16.6,-19.3,-18.7,-11.3,-1.16500e-09,101> mat1',
                                  'nu1': ' 1,2,3,4,5,6,7,8,9,10,11,12> mat1'})
    assert_equal(trans.transInput[0:78], "C ****************************************************************************")
    assert_equal(trans.transInput[-12:], "2.010000e+01")
