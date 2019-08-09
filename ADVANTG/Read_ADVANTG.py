# -*- coding: utf-8 -*-
"""
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
import pandas as pd
import random

import Gnowee
from ObjectiveFunction import ObjectiveFunction
from GnoweeUtilities import ProblemParameters

from GnoweeHeuristics import GnoweeHeuristics
from GnoweeHeuristics import disc_levy_flight
from GnoweeHeuristics import cont_levy_flight

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
filename='runCadis.adv'

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
    dict={};dict2={}
    txtfile=open(filename)
    txt=txtfile.read()
    #    [f1, f2]=txt.split('***********************************')
    i=txt.find('variables')
    txtfile.seek(i)
    print(txtfile.readline())
    for line in txtfile:
        if line[0]=='#':#.rstrip('\n') !='C':
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
#print(f2)
    return dict, dict2, f2

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
                lb,ub = vals.split(':'); dv=vals; dt='continuous'
            elif ',' in vals:
                disc=vals.split(',')
                lb,ub=disc[0],disc[-1]; dt='discrete'
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

#df=df.sort_values('data type')

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
                dv=vals; dt='continuous'
            elif ',' in vals:
                disc=vals.split(',')
                dt='discrete'
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
population=gh.population

def changedf(df, dfd):
    for p in range(1,population+1):
        for f in list(df.index):
            if df.loc[f].iloc[1]=='integer':
                df.at[f,'pop' + str(p)]=df.loc[f].iloc[2]
            elif df.loc[f].iloc[1]=='discrete':
                r=random.randint(0,len(df.loc[f].iloc[0])-1)
                df.at[f,'pop' + str(p)]=df.loc[f].iloc[0][r]
                for ff in list(dfd.index):
                    if f+'\r'==str(dfd.loc[ff].iloc[3]).replace(' ',''):
                        dfd.at[ff,'pop' + str(p)]=dfd.loc[ff].iloc[0][r]
            elif df.loc[f].iloc[1]=='continuous':
                if df.loc[f].iloc[4]==' decimal\r':
                    ru=round(random.uniform(df.loc[f].iloc[2],df.loc[f].iloc[3]),2)
                    df.at[f,'pop' + str(p)]=ru
                elif df.loc[f].iloc[4]==' integer\r':
                    ru=random.randint(df.loc[f].iloc[2],df.loc[f].iloc[3])
                    df.at[f,'pop' + str(p)]=ru
    elif df.loc[f].iloc[1]=='binary':
        r=random.randint(0,1)
        df.at[f,'pop' + str(p)]=r
            elif df.loc[f].iloc[1]=='combinatory':
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
return df, dfd

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

def changef2(df, dfd, f2, pp):
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
dfd=dict2dfd(dict2)
df, dfd= changedf(df, dfd)



file = open("mcnp.txt", "w")
file.write(f2)
file.close()

'''
    Contains: @f2new produces as  many MCNP input files as we have population lists in pop
    '''

for p in range(1,population+1):
    pp='pop'+str(p)
    f2new= changef2(df, dfd, f2, pp)
    file = open("mcnp_pop"+str(p)+".txt", "w")
    file.write(f2new)
    file.close()

dfp=df.drop(columns=['pop1'])
dfp.to_csv('params.csv')

'''
    #\gh=GnoweeHeuristics(objective=None, constraints=[], lowerBounds=df['lower boundary'].values,
    upperBounds=df['upper boundary'].values,varType=df['data type'].astype(str).values.tolist(),
    pltTitle='', histTitle='', varNames=df.index.astype(str).values.tolist())
    print gh
    (timeline) = Gnowee.main(gh)
    print '\nThe result:\n', timeline[-1]
    '''
'''
    we generate here the array of the population lists
    array=df.iloc[: , 5:15].values
    '''

array=df.iloc[: , 5:15].values.transpose()

'''
    Look into calling Levy flight for continuous and discrete
    levy=disc_levy_flight(array)
    '''

'''
    After from Gnowee new list of arrays reinsert into dataframe:
    df.iloc[:,5:15]=array
    '''

'''
    gh=GnoweeHeuristics(objective=None, constraints=[], lowerBounds=dfc['lower boundary'].values,upperBounds=dfc['upper boundary'].values,varType=dfc['data type'].astype(str).values.tolist(),pltTitle='', histTitle='', varNames=dfc.index.astype(str).values.tolist())
    # working example
    #gnowee=ProblemParameters(objective=None, constraints=[],lowerBounds=[1,2],upperBounds=[3,4],varType=['c','c','d'],discreteVals=[['10','11','12','13']],pltTitle='', histTitle='', varNames=['a','b','c'])
    '''

