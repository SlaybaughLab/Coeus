#######################################################################################################
#
# Module : test_UserInputs.py
#
# Contains : Routines to test UserInputs module
#
# Author : James Bevins
#
# Last Modified: 18Aug19
#
#######################################################################################################

import os 

import numpy as np

from UserInputs import UserInputs
from nose.tools import assert_equal
    
#-------------------------------------------------------------------------------------------------------------#         
# Assumed inputs
path = os.getcwd() +'\\Tests\\files_test_Coeus\\'
inputFname = path + 'test_user_inputs.txt'
print(inputFname)

#-------------------------------------------------------------------------------------------------------------#     
# Test the default object creation
def test_UserInputs():
    inputs = UserInputs(coeusInputPath=inputFname)
    assert_equal(inputs.coeusInput, path+"test_user_inputs.txt")
    if os.path.isfile(inputs.coeusInput):
        print("\nLoading Coeus input file located at: {}".format(
                                                           inputs.coeusInput))
    else:
        print("\nNo user supplier input file located a: {}".format(
                                                           inputs.coeusInput))
    assert_equal(inputs.transInput, None)
    assert_equal(inputs.advantgInput, None)
    assert_equal(inputs.code,"mcnp6")
   
# Test the file read
def test_read_inputs():
    inputs = UserInputs(coeusInputPath=inputFname)
    objFunc = inputs.read_inputs()

    # Test UserInputs object
    assert_equal(inputs.coeusInput, path+"test_user_inputs.txt")
    assert_equal(inputs.transInput, "test_mcnp.inp")
    assert_equal(inputs.advantgInput, "test_advantg.inp")
    assert_equal(inputs.code, "mcnp6")

    # Test Objectives object
    assert_equal(objFunc.func.__name__, "relative_least_squares")
    assert_equal(objFunc.funcTally,"24")
    assert_equal(objFunc.objType, "spectrum")
    assert_equal(objFunc.objForm, 0)
    np.testing.assert_equal(objFunc.objective, np.asarray(
                                    [[4.1399e-07, 4.6800e-15], 
                                     [1.1253e-06, 3.1300e-13],
                                     [3.0590e-06, 7.9900e-12],
                                     [1.0677e-05, 4.4200e-11],
                                     [2.9023e-05, 1.0700e-10],
                                     [1.0130e-04, 6.6300e-10],
                                     [2.7536e-04, 2.2000e-09],
                                     [5.8295e-04, 6.5900e-09],
                                     [1.2341e-03, 2.0900e-08],
                                     [3.3546e-03, 1.3600e-07],
                                     [1.0333e-02, 1.1500e-06],
                                     [2.1875e-02, 3.9500e-06],
                                     [2.4788e-02, 1.3500e-06],
                                     [3.4307e-02, 5.2900e-06],
                                     [5.2475e-02, 1.3400e-05],
                                     [1.1109e-01, 6.2000e-05],
                                     [1.5764e-01, 6.1000e-05],
                                     [2.4724e-01, 1.2500e-04],
                                     [3.6883e-01, 1.6400e-04],
                                     [5.5023e-01, 2.2000e-04],
                                     [6.3928e-01, 9.4200e-05],
                                     [7.4274e-01, 9.7300e-05],
                                     [8.2085e-01, 6.6100e-05],
                                     [9.6164e-01, 1.0600e-04],
                                     [1.1080e+00, 9.5000e-05],
                                     [1.4227e+00, 1.6600e-04],
                                     [1.8268e+00, 1.6200e-04],
                                     [2.3069e+00, 1.4300e-04],
                                     [2.3852e+00, 1.9600e-05],
                                     [3.0119e+00, 1.2900e-04],
                                     [4.0657e+00, 1.3300e-04],
                                     [4.7237e+00, 4.8200e-05],
                                     [4.9659e+00, 1.3200e-05],
                                     [6.3763e+00, 4.6400e-05],
                                     [7.4082e+00, 1.5100e-05],
                                     [8.1873e+00, 6.4900e-06],
                                     [9.0484e+00, 4.7500e-06],
                                     [1.0000e+01, 3.6900e-06],
                                     [1.1052e+01, 3.5500e-06],
                                     [1.2214e+01, 4.4300e-06],
                                     [1.2523e+01, 1.3400e-06],
                                     [1.3840e+01, 5.3200e-05],
                                     [1.4191e+01, 1.4400e-04],
                                     [1.4918e+01, 8.1900e-05],
                                     [1.6905e+01, 7.1500e-08],
                                     [1.9640e+01, 4.9100e-09]]))
    
   