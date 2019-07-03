# -*- coding: utf-8 -*-
"""!
Created on Sat Mar 24 08:15:15 2018

@author: Sandra Bogetic
"""
'''
Contains: This part of the code, is the beginning of a new COEUS interface system with Gnowee. 
          The idea here is to generalize to code so that any MCNP input file can be read as complex 
          as the user need to. Also ideally we are introducing the possibilities of having
          any type of ETA shape/contrainst/variables that can include also variables that are 
          external to the ETA itself. 
'''
import os
import pandas as pd
import numpy as np
import random
import GnoweeHeuristics
from GnoweeHeuristics import GnoweeHeuristics
#from GnoweeHeuristics import disc_levy_flight
#from GnoweeHeuristics import cont_levy_flight

"""
Contains: The file that is read here is an MCNP looking like file that contains to sections:
          - variables, those are identified as following: "# <>"  
          - The variables can be discrete(','), continuous(':'), binary(':'), 
            combinatorial.
          - The number type available are: integer and decimal (the decimal goes down to mm)
          - variables starting with '@' is if they are dependend to another variable
            that is defined at the end of the sentence description. E.g. if a material is 
            a discrete variable, it's important to make sure thaat the right density is associated
            with the right material
         
          @filename:choose an input file from inside the same folder
"""

import logging

module_logger = logging.getLogger('MCNP_interface')


import os
import sys
import copy as cp

from pandas.util.testing import assert_frame_equal

# choose an input file from inside the same folder

filename='MCNP_input.txt'

"""
Contains: this section imports the first part of the input file
          @dict : dictionary containing all the variables, the dictionary 
          with all the variables that are included in the first section of the file.
          The datafram identifies the type of variables that we have.
          Second part:
          
          @f2: text string, this is the MCNP section of the file with the variables within
          this section returns the dictionary with the variables, the MCNP input file
          with the variables inside

"""

def importfile(filename):
# Separation to two sections
    dict={};dict2={}
    
    try:
        # Open MCNP file
        
        txtfile=open(filename)
        txt=txtfile.read()
        i=txt.find('variables')
        txtfile.seek(i)
        
        print(txtfile.readline())
        
        for line in txtfile:
            
            if line[0]=='#':
            #.rstrip('\n') !='C':
                            line=line.rstrip('\n').replace('# <','').replace(' =','=')
                        
                            #val, unit=line.split('>')               
                            variable,value=line.split('=')
                            dict[variable]=value
            elif line[0]=='@':
                line=line.rstrip('\n').replace('@ <','').replace(' =','=')
                variable,value=line.split('=')
                dict2[variable]=value
            elif line[0]=='C': break  
            
        print(dict)
        
        j=txt.find('C  Cell Cards')
        txtfile.seek(j)
        print(txtfile.readline())       
        f2=txtfile.read()
        
        txtfile.close()
        
        return dict, dict2, f2
    
    #Close the file
    
        
    except IOError as e:
            module_logger.error( "I/O error({0}): {1}".format(e.errno, e.strerror)) 
            module_logger.error("File not found was: {0}".format(filename))  
            
             # Test that the file closed
    assert txtfile.closed==True, "File ({}) did not close properly.".format(filename)

"""
Contains: this section converts the dictionary into a dataframe containing columns for
          the original values, it's type, lower and upper boundary, an (empty) column for a randomly 
          selected value and unit
          @df : dataframe: contains columns where for each of the variables we have: name, type, boundaries, type of number
          (units also if one needs)
"""

def dict2df(dict):
    
    df=pd.DataFrame(index=list(dict.keys()),columns=['values','data type','lower boundary','upper boundary','number type'])
    lowb=[];upb=[];dtype=[];discvals=[];units=[]
    for vals in list(dict.values()):
        vals, unit=vals.split('>')
        unit=unit.replace('  ','')
        try:
            lb=float(vals); ub=float(vals); dv=float(vals); dt='integer'
        except:
            vals=vals.replace(' ','')
            if ':' in vals:
                lb,ub = vals.split(':'); dv=vals; dt='c'  
            elif ',' in vals:
                disc=vals.split(',')
                lb,ub=disc[0],disc[-1]; dt='d'
                dv=list(int(x) for x in  disc)   
            elif vals.isdigit() == True and vals[0] == '0':
                lb,ub=vals[0],vals[-1]; dv=list([int(lb), int(ub)]); dt='binary'
            elif vals.isdigit() == True:
                lb,ub=vals[0],vals[-1]; dv=list([int(lb), int(ub)]); dt='combinatory'
            else:
                lb='nan'; ub='nan'; dv=vals; dt='function'
                
        lowb.append(float(lb));upb.append(float(ub));dtype.append(dt)
        units.append(unit);discvals.append(dv)
                
    df['lower boundary']=lowb; df['upper boundary']=upb; df['data type']=dtype
    df['number type']=units ;df['values']=discvals
    
    return df
'''
Contains: @dict2dfd cretes a dataframe for those variables that are dependemnt
          on the actual variables
          E.g. If materials are variables in the MCNP file as m1,m2 and so on we need to make
          sure in MCNP we have the right density added
          Several dependent variable can be created
          
          Return @dfd : dependent data frame
'''

def dict2dfd(dict2):
    dfd=pd.DataFrame(index=list(dict2.keys()),columns=['values','data type','random value','variable'])
    dtype=[];discvals=[];units=[]
    for vals in list(dict2.values()):
        vals, unit=vals.split('>')
        unit=unit.replace('  ','')
        try:
            dv=float(vals); dt='integer'
        except:
            vals=vals.replace(' ','')
            if ':' in vals:
                dv=vals; dt='c'  
            elif ',' in vals:
                disc=vals.split(',')
                dt='d'
                dv=list(float(x) for x in  disc)   
            elif vals.isdigit() == True and vals[0] == '0':
                lb,ub=vals[0],vals[-1]; dv=list([int(lb), int(ub)]); dt='binary'
            elif vals.isdigit() == True:
                lb,ub=vals[0],vals[-1]; dv=list([int(lb), int(ub)]); dt='combinatory'
            else:
                lb='nan'; ub='nan'; dv=vals; dt='function'
        dtype.append(dt)
        units.append(unit);discvals.append(dv)
                
    dfd['data type']=dtype
    dfd['variable']=units ;dfd['values']=discvals
    return dfd


"""
Contains: Here @changedf assigns a random value within the possible boundaries to each variable,
          Gnowee needd to receive an array 'pop' from COEUS that is made up of several 
          list of pop variables, 25 are recommended for actual run, 10 is enough for testing
          
          @df : it is still the dataframe with added columns of random values for each variable 
          up to 10 or 25 populations
"""

#population=gh.population

def changedf(df, dfd):
    for p in range(1,11):
        for f in list(df.index):
            if df.loc[f].iloc[1]=='integer':
                df.at[f,'pop' + str(p)]=df.loc[f].iloc[2]
            elif df.loc[f].iloc[1]=='d':
                random.seed(9001)
                r=random.randint(0,len(df.loc[f].iloc[0])-1)
                df.at[f,'pop' + str(p)]=df.loc[f].iloc[0][r]
                for ff in list(dfd.index):
                    if f+'\r'==str(dfd.loc[ff].iloc[3]).replace(' ',''):
                        dfd.at[ff,'pop' + str(p)]=dfd.loc[ff].iloc[0][r] 
            elif df.loc[f].iloc[1]=='c':
                if df.loc[f].iloc[4]==' decimal\r':
                    random.seed(9001)
                    ru=round(random.uniform(df.loc[f].iloc[2],df.loc[f].iloc[3]),2)
                    df.at[f,'pop' + str(p)]=ru
                elif df.loc[f].iloc[4]==' integer\r':
                    random.seed(9001)
                    ru=random.randint(df.loc[f].iloc[2],df.loc[f].iloc[3])
                    df.at[f,'pop' + str(p)]=ru
            elif df.loc[f].iloc[1]=='binary':
                random.seed(9001)
                r=random.randint(0,1)
                df.at[f,'pop' + str(p)]=r
            elif df.loc[f].iloc[1]=='combinatory':
                random.seed(9001)
                r=random.randint(0,1)
                df.at[f,'pop' + str(p)]=df.loc[f].iloc[r]
                
               

        for f in list(df.index):
            if df.loc[f].iloc[1]=='function':
                fct=str(df.loc[f].iloc[0])
                for ff in list(df.index):
                    fct=fct.replace(ff,str(df.loc[ff].iloc[-1]))
                try:
                    fct=eval(fct)
                    df.at[f,'pop' + str(p)]=fct
                except:
                    print('there is an undefined variable in the function')
#       
                    assert_frame_equal(df, dfd, check_dtype=False)
    return df, dfd

"""
create input variables for gnowee Heuristics to call gnowee and create new values:
"""
def df2gnowee(df):
    vT=df['data type'].astype(str).values.tolist()
    lB=df['lower boundary'].loc[df['data type']=='c'].values.tolist()
    uB=df['upper boundary'].loc[df['data type']=='c'].values.tolist()
    vN=df.index.astype(str).values.tolist()
    dV=df['values'].loc[df['data type']=='d'].tolist()
    
    gh=GnoweeHeuristics(objective=None, constraints=[], lowerBounds=lB, upperBounds=uB, varType=vT, discreteVals=dV, pltTitle='', histTitle='', varNames=vN)
    
    pop_c=df.iloc[:,5:].loc[df['data type']=='c'].values.transpose().tolist()
    pop_d=df.iloc[:,5:].loc[df['data type']=='d'].values.transpose().tolist()
    
#    GnoweeHeuristics.cont_levy_flight(gh, pop_c)
#    GnoweeHeuristics.disc_levy_flight(gh, pop_d)
    
    
    return gh, pop_c, pop_d

'''
Contains: @changef2 takes @f2, the MCNP input file with the variables and subsitutes 
          each variable in the right location with the values from 'pop'
          
          @df : data frame for the variables
          @dfd : dependent variables
          @f2 : MCNP input file with the variables
          @pp : population number
          
Returns:  @f2new : MCNP files with the new values, here done first with the random values given to 
         Gnowee, then the same is done using the values Gnowee

'''

def changef2(df, f2, pp):
    f2new=f2
    for f in list(df.index):
        f2new=f2new.replace(f,str(df.at[f,pp]))
        
    for ff in list(dfd.index):
        f2new=f2new.replace(ff,str(dfd.at[ff,pp]))
       
    while f2new.find('(')>0:
            b=f2new.find('(')
            e=f2new.find(')')
            little_math=f2new[b:e+1].replace('<','').replace('>','')
            try:
                sol=round(eval(little_math),2)
                f2new=f2new.replace(f2new[b:e+1],str(sol))
            except: 
                print('there is an undefined variable between brackets - break')
                break
                

    return f2new


dict, dict2, f2 = importfile(filename)
df=dict2df(dict)
df=df.sort_values('data type')
dfd=dict2dfd(dict2)
df, dfd= changedf(df, dfd)



#CREATE GNOWEE VARIABLES

gh, pop_c, pop_d=df2gnowee(df)

print(gh)
#
array=df.iloc[: , 5:].values.transpose()

print(array)

# CALL LEVY FLIGHT OPTIONS


levy_disc=GnoweeHeuristics.disc_levy_flight(gh,array)
print(levy_disc)

array_d=pd.DataFrame(levy_disc[0])
array_d=np.transpose(array_d)
df.iloc[:,5:]=array_d.values

array_new=df.iloc[: , 5:].values.transpose()

levy_cont=GnoweeHeuristics.cont_levy_flight(gh,array_new)
print(levy_cont)

array_c=pd.DataFrame(levy_disc[0])
array_c=np.transpose(array_c)
array_c=array_c.round(2)
df.iloc[:,5:]=array_c.values

"""
    assert_frame_equal(df, dfp, check_names=True)
    assert_frame_equal(df, dfp, check_column_type=True)     
    assert_frame_equal(df, dfp, check_check_frame_type=True)
"""
# ***********************************************************

file = open("mcnp.txt", "w") 
file.write(f2) 
file.close()

assert file.closed==True, "File ({}) did not close properly.".format(filename)


dfp=df

'''
Contains: roduces as  many MCNP input files as we have population lists in pop
'''

for p in range(1,11):
    pp='pop'+str(p)
    f2new= changef2(df, f2, pp)
    file = open("mcnp_pop"+str(p)+".txt", "w") 
    file.write(f2new) 
    file.close() 
    
    dfp=dfp.drop(columns=['pop'+ str(p)])

assert file.closed==True, "File ({}) did not close properly.".format(filename)

'''
Genertate SBATCH file
'''


file = open("run_mcnp.sh", "w")
file.write("#!/bin/sh\n\n")
file.write("#SBATCH --time=02:30:00\n")
file.write("# Job name:\n")
file.write("#SBATCH --job-name=mc20\n")
file.write("# QoS:\n")
file.write("#SBATCH --partition=savio2\n")
file.write("# Account:\n")
file.write("#SBATCH --ntasks=20\n")
file.write("#SBATCH --output=arrayJob_%A_%a.out\n")
file.write("#SBATCH --error=arrayJob_%A_%a.err\n")
file.write("# Array:\n")
file.write("#SBATCH --array=1,2,3,4,5,7,8,9,10\n\n")
file.write("module load openmpi\n")
file.write("mpirun mcnp6.mpi i=mcnp_pop$SLURM_ARRAY_TASK_ID.txt o=mcnp_pop$SLURM_ARRAY_TASK_ID.out\n")
file.close()
    
# Run SBATCH with MCNP


cmd = 'sbatch run_mcnp.sh'
print cmd
os.system(cmd)


