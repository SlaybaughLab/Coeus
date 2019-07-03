# -*- coding: utf-8 -*-
"""!
Created on Sat Mar 24 08:15:15 2018

@author: Sandra Bogetic
"""
import os
import pandas as pd
import numpy as np
import random
import GnoweeHeuristics
from GnoweeHeuristics import GnoweeHeuristics
#from GnoweeHeuristics import disc_levy_flight
#from GnoweeHeuristics import cont_levy_flight

from Interface_test import importfile, dict2df, dict2dfd, changedf
from Interface_test import df2gnowee, changef2 

from unittest import TestCase
import nose
from nose.plugins.skip import SkipTest
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, 
    assert_almost_equal, assert_true, assert_false, assert_in
from numpy.testing import assert_array_equal

#---------------------------------------------------------------



df.isnull().values.any()

df_test=pd.DataFrame(columns=['values','data type','lower boundary','upper boundary','number type','pop1','pop2','pop3','pop4','pop5','pop6','pop7','pop8','pop9','pop10'])
df_test=pd.DataFrame(index=['r1', 'h2', 'h3', 'h1', 'h6', 'h4', 'h5', 'd1', 'mat1'])


test_value_mat1 = df.at['mat1', 'data type']
test_value_h1 = df.at['h1', 'data type']
test_value_d1 = df.at['d1', 'data type']
test_value_h2 = df.at['h2', 'data type']
test_value_h3 = df.at['h3', 'data type']
test_value_h4 = df.at['h4', 'data type']
test_value_h5 = df.at['h5', 'data type']
test_value_h6 = df.at['h6', 'data type']
test_value_r1 = df.at['r1', 'data type']


np.array_equal(df.columns,df_test.columns)
np.array_equal(df.index,df_test.index)
np.array_equal(df['lower boundary'].values,[0.5, 1. , 1. , 1. , 1. , 1. , 1. , 1. , 1.])
    
def test_variables():    
    assert_equal(test_value_mat1,"d")
    assert_equal(test_value_h1,"c")
    assert_equal(test_value_d1,"c")
    assert_equal(test_value_r1,"c")
    assert_equal(test_value_h2,"c")
    assert_equal(test_value_h3,"c")
    assert_equal(test_value_h4,"c")
    assert_equal(test_value_h5,"c")
    assert_equal(test_value_h6,"c")


assert_array_equal(array_d, array_c)
assert_array_equal(df.iloc[:,6],array.transpose()[:,2])

df_test=df

#def Read_Tally_Output(path, tnum):

#assert isinstance(path, str)==True

