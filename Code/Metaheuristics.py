#######################################################################################################
#
# Module : Metaheuristics.py
#
# Contains : Gnowee Hybrid Metaheuristic Optimization algorithm methods.
#
# Author : James Bevins
#
# Last Modified: 20Aug16
#
#######################################################################################################

import logging
module_logger = logging.getLogger('Coeus.Metaheuristics')

import os
import sys
sys.path.insert(0,os.path.abspath(os.getcwd())+'/Sampling')

import copy as cp
import numpy as np
import SamplingMethods as sm

from random import random
from Gnowee_Utilities import Rejection_Bounds, Simple_Bounds, Parent
from math import sqrt, ceil, tan, radians
from Utilities import WeightedRandomGenerator
from scipy.stats import rankdata

#---------------------------------------------------------------------------------------#
def Mat_Levy_Flights(x,mats,mr,S,exclude):
    """
    Change cell materials based on Levy draw. The materials will be changed by either 
    a) using material library key list index numbers or b) moderating ratio (for both 1 and 14 MeV).  
    The choice will be based on a random number draw and be 33/33/33. 
    This provides for a decoupling to the moderation power to add a layer of more randomness, 
    while maintaining a physics based Levy flight process. 
   
    Parameters
    ==========
    x : list of parent objects
        The current parents representing system designs
    mats : dictionary of material objects
        A materials library containing all relevant nulcear data required to run radiation transport codes.
    mr : list of moderating ratio objects
        Contains the moderating ratios for the materials library used to guide the Levy flight
    S : Object    
        An object representing the settings for the optimization algorithm
    exclude : list    
        A list of materials to be excluded
   
    Optional
    ========   
   
    Returns
    =======
    tmp : list of parent objects
        The proposed parents representing new system designs
    """
    
    assert S.fl>=0 and S.fl <=1, 'The probability that a parent is used for global Levy search must exist on [0,1] & =%d' %S.fl
    
    tmp=[] # Local copy of parent that is modified
    keys=mats.keys()
    module_logger.debug("Keys: {}".format(keys))
        
    # Determine step size using Levy Flight
    for i in x:
        module_logger.debug("Parent materials: {}".format(i.geom.matls))
    step=sm.Levy(max(len(i.geom.matls) for i in x)-i.fixed_mats+1,len(x),alpha=S.a,gamma=S.g,n=S.n) #+1 b/c fill isn't counted
    module_logger.debug("The steps are: {}".format(step))
    module_logger.debug("{}, {}, {}, {}".format(S.a,S.g,S.n,1.0/S.sf))
    
    # Perform global search from fl parents
    for i in range(0,int(S.fl*S.p)):
        # Make a local copy
        tmp.append(cp.deepcopy(x[i]))
        
        # Select random number to determine permutation method.  
        # p<=0.33=Matl key index, 0.33<p<=0.66= 1 MeV Moderating Ratio, 0.66<p<=1.0= 14 MeV Moderating Ratio
        p=random()
        module_logger.debug("The decision variable p= {}".format(p))
        
        if p <= 0.33:
            #Calculate Levy flight based on material key index
            for j in range(tmp[-1].fixed_mats,len(tmp[-1].geom.matls)-1): #structural mats plus 1 void fill on end of list
                # Find current index of material
                ind=next((i for i, item in enumerate(keys) if item == tmp[-1].geom.matls[j]), -1)
                module_logger.debug("Step: {}, {}, {})".format(ind,int(step[i,j-tmp[-1].fixed_mats]),(ind+int(step[i,j-tmp[-1].fixed_mats]))%len(keys)))
                module_logger.debug("Old: {})".format(tmp[-1].geom.matls[j]))
                levy=(ind+int(step[i,j-tmp[-1].fixed_mats]))%len(keys)
                
                # Exclude Fissile Material from Geometry
                while keys[levy] in exclude:
                    resamp=sm.Levy(max(len(i.geom.matls) for i in x)-tmp[-1].fixed_mats+1,len(x),alpha=S.a,gamma=S.g,n=S.n) #+1 fill 
                    step[i,j-tmp[-1].fixed_mats]=resamp[i,j-tmp[-1].fixed_mats]
                    levy=(ind+int(step[i,j-tmp[-1].fixed_mats]))%len(keys)
                
                # Update material if a new material is selected
                if levy != ind:
                    tmp[-1].geom.matls[j]=keys[levy]    
                    module_logger.debug("New: {})".format(tmp[-1].geom.matls[j]))
                    for c in tmp[-1].geom.cells:
                        if c.m == j+1:
                            module_logger.debug("Levy: {})".format(levy))
                            module_logger.debug("keys[Levy]: {})".format(keys[levy]))
                            module_logger.debug("mats[keys[levy]]: {})".format(mats[keys[levy]]))
                            module_logger.debug("mats[keys[levy]].density: {})".format(mats[keys[levy]].density))
                            c.d=mats[keys[levy]].density

                module_logger.debug("New parent materials list: {})".format(tmp[-1].geom.matls))
        elif p > 0.33 and p <= 0.66:
            #Calculate Levy flight based on 1 MeV stopping ratio
            mr.sort(key=lambda x: x.mr_1MeV)
            for j in range(tmp[-1].fixed_mats,len(tmp[-1].geom.matls)-1): #structural mats plus 1 void fill on end of list
                # Find current index of material
                ind=next((i for i, item in enumerate(mr) if item.name == tmp[-1].geom.matls[j]), -1)                
                module_logger.debug("Step: {}, {}, {})".format(ind,int(step[i,j-tmp[-1].fixed_mats]),(ind+int(step[i,j-tmp[-1].fixed_mats]))%len(keys)))
                module_logger.debug("Old: {})".format(tmp[-1].geom.matls[j]))
                levy=(ind+int(step[i,j-tmp[-1].fixed_mats]))%len(keys)
                
                # Exclude Fissile Material from Geometry
                while mr[levy].name in exclude:
                    resamp=sm.Levy(max(len(i.geom.matls) for i in x)-tmp[-1].fixed_mats+1,len(x),alpha=S.a,gamma=S.g,n=S.n) #+1 fill 
                    step[i,j-tmp[-1].fixed_mats]=resamp[i,j-tmp[-1].fixed_mats]
                    levy=(ind+int(step[i,j-tmp[-1].fixed_mats]))%len(keys)
                
                # Update material if a new material is selected
                if levy != ind:
                    tmp[-1].geom.matls[j]=mr[levy].name   
                    module_logger.debug("New: {})".format(tmp[-1].geom.matls[j]))
                    for c in tmp[-1].geom.cells:
                        if c.m == j+1:
                            module_logger.debug("Levy: {})".format(levy))
                            module_logger.debug("mr[Levy]: {})".format(mr[levy]))
                            module_logger.debug("mats[mr[levy].name]: {})".format(mats[mr[levy].name]))
                            module_logger.debug("mats[mr[levy].name].density: {})".format(mats[mr[levy].name].density))
                            c.d=mats[mr[levy].name].density
                    
        elif p > 0.66 and p <= 1.0:
            #Calculate Levy flight based on 14 MeV stopping ratio
            mr.sort(key=lambda x: x.mr_14MeV)
            for j in range(tmp[-1].fixed_mats,len(tmp[-1].geom.matls)-1): #2 structural materials plus 1 void fill on end of list
                # Find current index of material
                ind=next((i for i, item in enumerate(mr) if item.name == tmp[-1].geom.matls[j]), -1)                
                module_logger.debug("Step: {}, {}, {})".format(ind,int(step[i,j-tmp[-1].fixed_mats]),(ind+int(step[i,j-tmp[-1].fixed_mats]))%len(keys)))
                module_logger.debug("Old: {})".format(tmp[-1].geom.matls[j]))
                levy=(ind+int(step[i,j-tmp[-1].fixed_mats]))%len(keys)
                
                # Exclude Fissile Material from Geometry
                while mr[levy].name in exclude:
                    resamp=sm.Levy(max(len(i.geom.matls) for i in x)-tmp[-1].fixed_mats+1,len(x),alpha=S.a,gamma=S.g,n=S.n) #+1 fill 
                    step[i,j-tmp[-1].fixed_mats]=resamp[i,j-tmp[-1].fixed_mats]
                    levy=(ind+int(step[i,j-tmp[-1].fixed_mats]))%len(keys)
                
                # Update material if a new material is selected
                if levy != ind:
                    tmp[-1].geom.matls[j]=mr[levy].name    
                    module_logger.debug("New: {})".format(tmp[-1].geom.matls[j]))
                    for c in tmp[-1].geom.cells:
                        if c.m == j+1:
                            module_logger.debug("Levy: {})".format(levy))
                            module_logger.debug("mr[Levy]: {})".format(mr[levy]))
                            module_logger.debug("mats[mr[levy].name]: {})".format(mats[mr[levy].name]))
                            module_logger.debug("mats[mr[levy].name].density: {})".format(mats[mr[levy].name].density))
                            c.d=mats[mr[levy].name].density
        else: 
            module_logger.error("p is out of bounds.")
            sys.exit()
        
    return tmp

#---------------------------------------------------------------------------------------#
def Cell_Levy_Flights(x,eta,S):
    """
     Cell Levy Flight: Change all cell and foil starting locations and cell deltas based on Levy draw. 
     The parameters modified are $z_{foil}$, $\Delta z_{hc}$, $r_{vc}$, $\Delta r_{vc}$, $z_{vc}$ , 
     and $\Delta z_{vc}$.
   
    Parameters
    ==========
    x : list of parent objects
        The current parents representing system designs
    eta : Object    
        An object representing the constraints for the eta design
    S : Object    
        An object representing the settings for the optimization algorithm
   
    Optional
    ========   
   
    Returns
    =======
    tmp : list of parent objects
        The proposed parents representing new system designs
    """
    assert S.fl>=0 and S.fl <=1, 'The probability that a parent is used for global Levy search must exist on [0,1] & =%d' %S.fl
    
    # Initialize lists
    tmp=[]
    used=[]
    
    for i in range(int(S.fl*S.p)):
        lb=[]
        ub=[]   
        cur_d=[] 
        new_d=[]
        
        # Make a local copy from a random parent
        r=int(np.random.rand()*S.p)
        while r in used:
            r=int(np.random.rand()*S.p)
        used.append(r)
        tmp.append(cp.deepcopy(x[r])) 
        
        # Determine step size using Levy Flight
        step=sm.Levy(1+4*eta.max_vert+eta.max_horiz,len(x),alpha=S.a,gamma=S.g,n=S.n)
        module_logger.debug("The steps for Cell_Levy_Flights are: {}\n".format(step))
        module_logger.debug("Step[0,1]= {}".format(step[0,1]))
        
        # Build design variable set from current parent
        # [foil_z, N_vert*(z, delz, r1, r2), N_horiz*(z)]
        prev_vert=''
        for s in tmp[i].geom.surfaces:
            if s.c=="NAS":
                module_logger.debug("Found NAS. VZ={} and Cell={}".format(s.vz,repr(s)))
                cur_d.append(s.vz)
                
                # Calculate Foil_Z Boundaries
                lb.append(max(eta.tcc_dist+eta.t_w+2*eta.t_nas+sum(eta.t_nas_f)+0.203,eta.tcc_dist+eta.t_w+(eta.r_nas-eta.r_f+eta.t_w)/tan(radians(90-eta.theta))))
                ub.append(eta.snout_dist-eta.t_m-2*eta.t_nas-sum(eta.t_nas_f)-0.203)
                
            elif s.c[0:4]=="vert":
                module_logger.debug("Found {}. VZ={}, HZ={}, and r={} and Cell={}".format(s.c,s.vz,s.hz,s.r,repr(s)))
                if prev_vert==s.c:
                    cur_d.append(s.r)
                    prev_vert=s.c
                else:
                    cur_d.append(s.vz)
                    cur_d.append(s.hz)
                    cur_d.append(s.r)
                    prev_vert=s.c
                    
            elif s.c[0:7]=="horiz #":
                module_logger.debug("Found {}. delZ={} and Cell={}".format(s.c,s.d,repr(s)))
                if s.c=="horiz #1":
                    cur_d.append(s.d-(eta.tcc_dist+eta.t_ds))
                else:
                    cur_d.append(s.d-prev_z)
                prev_z=s.d
                
        # Convert to numpy arrays  
        cur_d=np.asarray(cur_d)
        module_logger.debug("Design Variable set for parent #{} = {}\n".format(tmp[i].ident,cur_d))    
        
        # Update design variable set
        stepsize=1.0/S.sf*step[r,:]
        new_d=cur_d+stepsize
        module_logger.debug("Stepsize ={}".format(stepsize))
        module_logger.debug("Updated Variable set for parent #{} = {}\n".format(tmp[i].ident,new_d))
        
        # Calculate Vertical Cell Boundaries (z, delz, r1, r2)
        for i in range(0,eta.max_vert):
            if new_d[i*4+3] > new_d[i*4+4]:
                t=new_d[i*4+4]
                new_d[i*4+4]=new_d[i*4+3]
                new_d[i*4+3]=t
            lb.append(eta.tcc_dist+eta.t_ds)
            lb.append(0.0000001)
            lb.append(0.0000001)
            if new_d[i*4+3] > 0.0:
                lb.append(new_d[i*4+3])  
            else:
                lb.append(0.0000001) 
            ub.append(eta.snout_dist-eta.t_c-0.00001)
            ub.append(eta.snout_dist-eta.t_c-new_d[i*4+1])
            ub.append(eta.r_o)
            ub.append(eta.r_o)
        
        # Calculate horizontal Z values and bounds
        new_d[-eta.max_horiz]=new_d[-eta.max_horiz]+(eta.tcc_dist+eta.t_ds)
        lb.append(eta.tcc_dist+eta.t_ds)
        ub.append(eta.snout_dist-eta.t_c)
        for i in range(-eta.max_horiz+1,0):
            new_d[i]=new_d[i]+new_d[i-1]
            lb.append(0.0000001)
            ub.append(eta.snout_dist-eta.t_c-new_d[i-1])
        lb=np.array(lb)
        module_logger.debug("Lower Bounds ={}\n".format(lb)) 
        ub=np.array(ub)
        module_logger.debug("Upper Bounds ={}\n".format(ub))
            
        # Applying boundaries check
        new_d=Rejection_Bounds(cur_d,new_d,stepsize,lb,ub,S) 
        module_logger.debug("Post Boundary Variable set for parent #{} = {}\n".format(tmp[i].ident,new_d))
        
        # Update parents with new design set
        # [foil_z, N_vert*(z, delz, r1, r2), N_horiz*(z)]
        prev_vert=''
        prev_z=0.0
        for s in tmp[i].geom.surfaces:
            if s.c[-4:]=="foil" or s.c[-4:]=="TOAD":
                s.vz=s.vz+new_d[0]-cur_d[0]
            if s.c=="NAS":
                s.vz=new_d[0]
                new_d=new_d[1:]
                nas_vz=s.vz
            if s.c=="Holder":
                s.vz=nas_vz-eta.t_h
            if s.c=="Holder Fill":
                s.vz=nas_vz
                
            elif s.c[0:4]=="vert":
                if prev_vert==s.c:
                    s.vz=new_d[0]
                    s.hz=new_d[1]
                    s.r=new_d[3]
                    new_d=new_d[4:]
                    prev_vert=s.c
                else:
                    s.vz=new_d[0]
                    s.hz=new_d[1]
                    s.r=new_d[2]
                    prev_vert=s.c
                    
            elif s.c[0:7]=="horiz #":
                if new_d[0]<prev_z:
                    if prev_z+0.1<eta.snout_dist-eta.t_c:
                        s.d=prev_z+0.1
                    else:
                        s.d=prev_z
                else:
                    s.d=new_d[0]
                new_d=new_d[1:]
                prev_z=s.d
                
        module_logger.debug("For i={}, ident={}, the parent={}\n".format(i,tmp[i].ident,str(tmp[i].geom)))             
    return tmp

#---------------------------------------------------------------------------------------#
def Elite_Crossover(x,mr,eta,mats,S,exclude):
    """
    Change the materials between the top parent and an elite parent based on moderating ratio. 
    The materials will be changed by moderating ratio (for both 1 and 14 MeV).  
    The choice will be based on a random number draw and be 50/50.
   
    Parameters
    ==========
    x : list of parent objects
        The current parent representing system designs
    mr : list of moderating ratio objects
        Contains the moderating ratios for the materials library used to guide the Levy flight
    eta : Object    
        An object representing the constraints for the eta design
    mats : dictionary of material objects
        A materials library containing all relevant nulcear data required to run radiation transport codes.
    S : Object    
        An object representing the settings for the optimization algorithm
    exclude : list    
        A list of materials to be excluded
   
    Optional
    ========   
   
    Returns
    =======
    tmp : list of parent objects
        The proposed parent representing new system design
    """    
    assert S.fe>=0 and S.fe <=1, 'The probability that a parent is elite must exist on [0,1] & =%d' %S.fe

    
    # Initialize variables
    golden_ratio=(1.+sqrt(5))/2.  # Used to bias mutation strategy
    keys=mats.keys()
    module_logger.debug("Keys: {}".format(keys))
        
    # Choose nests for crossover
    top=[]
    top.append(cp.deepcopy(x[0]))
    r=int(np.random.rand()*S.p*S.fe)
    if r==0:
        r=1
    top.append(cp.deepcopy(x[r]))
    
    # Create list of materials for top parent
    t_keys=[]
    t_mr=[]
    for c in top[0].geom.cells:
        module_logger.debug("Top parent #{}={}".format(top[0].ident,repr(c)))
        if c.comment=="vert" or c.comment=="horiz":
            t_keys.append(top[0].geom.matls[c.m-1])
            t_mr.append(next((item for i, item in enumerate(mr) if item.name == t_keys[-1]), -1))
    module_logger.debug("Top Parent[{}] cell material indexes = {}".format(top[0].ident,t_keys))
    module_logger.debug("Moderating ratios for top parent[{}] = {}\n".format(top[0].ident,t_mr))
        
    # Select MR to use
    p=random()
        
    # Create list of materials for random parent
    r_keys=[]
    r_mr=[]
    for c in top[1].geom.cells:
        module_logger.debug("Random parent #{}={}".format(top[1].ident,repr(c)))
        if c.comment=="vert" or c.comment=="horiz":
            r_keys.append(top[1].geom.matls[c.m-1])
            r_mr.append(next((item for i, item in enumerate(mr) if item.name == r_keys[-1]), -1))
    module_logger.debug("Random parent[{}] cell material indexes = {}".format(top[1].ident,r_keys))
    module_logger.debug("Moderating ratios for random parent[{}] = {}\n".format(top[1].ident,r_mr))
        
    # Calculate the mutated material
    new_mat=[]
    for i in range(0,len(t_keys)):
        if p <= 0.5:
            new_mr=t_mr[i].mr_1MeV+(r_mr[i].mr_1MeV-t_mr[i].mr_1MeV)/golden_ratio
            mr.sort(key=lambda x: x.mr_1MeV)
            j=0
            while mr[j].mr_1MeV < new_mr:
                j+=1
            if abs(new_mr-mr[j].mr_1MeV)<=abs(new_mr-mr[j-1].mr_1MeV):
                new_mat.append(mr[j].name)
            elif abs(new_mr-mr[j-1].mr_1MeV)<abs(new_mr-mr[j].mr_1MeV):
                new_mat.append(mr[j-1].name)
                j-=1
            else:
                module_logger.error("Apparently the programmer did not consider all cases.  Sucks to be you")
        elif p <= 1.0:
            new_mr=t_mr[i].mr_14MeV+(r_mr[i].mr_14MeV-t_mr[i].mr_14MeV)/golden_ratio
            mr.sort(key=lambda x: x.mr_14MeV)
            j=0
            while mr[j].mr_14MeV < new_mr:
                j+=1
            if abs(new_mr-mr[j].mr_14MeV)<=abs(new_mr-mr[j-1].mr_14MeV):
                new_mat.append(mr[j].name)
            elif abs(new_mr-mr[j-1].mr_14MeV)<abs(new_mr-mr[j].mr_14MeV):
                new_mat.append(mr[j-1].name)
                j-=1
            else:
                module_logger.error("Apparently the programmer did not consider all cases.  Sucks to be you.")
        else:
            module_logger.error("Moderating Ratio is sampling outside the interval from 0-1")
        
        # Exclude excluded materials from geometry
        while new_mat[-1] in exclude:
            if j+1<len(mr):
                new_mat[-1] = mr[j+1].name
                j+=1
            else:
                new_mat[-1] = mr[0].name
                j=0

    module_logger.debug("The new materials are = {}\n".format(new_mat))
                    
    # Update material if a new material is selected
    j=top[0].fixed_mats
    for c in top[0].geom.cells:
        if c.comment=="vert" or c.comment=="horiz":
            if j < len(top[0].geom.matls):
                top[0].geom.matls[j]=new_mat[0]
            else:
                top[0].geom.add_matls(mats,new_mat[0])
            c.m=j+1
            c.d=mats[new_mat[0]].density
            new_mat=new_mat[1:]    
            j+=1
        elif c.comment=="eta fill":
            if top[0].geom.matls[c.m-1]!=eta.fill_mat:
                module_logger.debug("The materials before are = {}\n".format(top[0].geom.matls))
                top[0].geom.add_matls(mats,eta.fill_mat)
                module_logger.debug("The materials after are = {}\n".format(top[0].geom.matls))
                c.m=j+1
                c.d=mats[eta.fill_mat].density
        
    # Update identity and nps
    top[0].ident=top[1].ident
    top[0].rset.nps=top[1].rset.nps
    tmp=[]
    tmp.append(cp.deepcopy(top[0]))
         
    return tmp

#---------------------------------------------------------------------------------------#
def Partial_Inversion(x,mr,mats,S):
    """
    Invert materials based on moderating ratio gradient.  I.e. Pick random layer l.  
    If layer l+1 has the next highest (or lowest) moderating ratio do nothing.  
    Otherwise invert the layer(s) between layer l and the layer with the next higher (or lower) ratio.
   
    Parameters
    ==========
    x : list of parent objects
        The current parents representing system designs
    mr : list of moderating ratio objects
        Contains the moderating ratios for the materials library used to guide the inversion
    mats : dictionary of material objects
        A materials library containing all relevant nulcear data required to run radiation transport codes.
    S : Object    
        An object representing the settings for the optimization algorithm
   
    Optional
    ========   
   
    Returns
    =======
    tmp : list of parent objects
        The proposed parents representing new system designs
    """    
    tmp=[]
    
    for i in range(0,S.p):
        tmp.append(cp.deepcopy(x[i]))
        
        # Create list of materials and moderating ratios
        p=random()
        keys=[]
        c_mr=[]
        module_logger.debug("The starting matls list is = {}\n".format(tmp[-1].geom.matls))
        for c in tmp[-1].geom.cells:
            if c.comment=="vert" or c.comment=="horiz":
                keys.append(tmp[-1].geom.matls[c.m-1])
                if p<=0.5:
                    c_mr.append(next((item.mr_1MeV for i, item in enumerate(mr) if item.name == keys[-1]), -1))
                elif p<=1.0:
                    c_mr.append(next((item.mr_14MeV for i, item in enumerate(mr) if item.name == keys[-1]), -1))
        module_logger.debug("Parent[{}] cell material indexes = {}".format(tmp[-1].ident,keys))
        module_logger.debug("Moderating ratios for parent[{}] = {}\n".format(tmp[-1].ident,c_mr))
        
        old_keys=cp.copy(keys)
        
        # Select random point, walk through until inversion point found or end 
        loc=int(random()*len(c_mr))
        s=sorted(range(len(c_mr)), key=lambda k: c_mr[k])
        if loc+1 != len(s):
            while s[loc]+1 == s[loc+1]:
                loc+=1
                if loc+1 == len(s):
                    loc=-1
                    break
        else:
            loc=-1
        module_logger.debug("Loc={} and the sorted morderating ratios are = {}\n".format(loc,s))
        
        # Invert materials and correct cell assignments
        try:
            ind=s.index(s[loc]+1)
        except:
            ind=-1
        if ind < loc and loc != -1 and s[loc] != len(s)-1:
            t=keys[loc+1]   
            keys[loc+1]=keys[ind] 
            keys[ind]=t
            module_logger.debug("The index of s[loc]+1={} and the swapped sub-matls list is = {}\n".format(ind,keys))
            tmp[-1].geom.matls[-len(keys)-1:-1]=keys
            module_logger.debug("The reversed matls list is = {}\n".format(tmp[-1].geom.matls)) 
            
            # Update Cell Densities
            j=tmp[-1].fixed_mats
            for c in tmp[-1].geom.cells:
                if c.comment=='horiz' or c.comment=='vert':
                    c.d=mats[tmp[-1].geom.matls[j]].density
                    c.m=j+1
                    j+=1
        elif loc != -1 and s[loc] != len(s)-1:
            keys[loc+1:ind+1]=reversed(keys[loc+1:ind+1])
            module_logger.debug("The index of s[loc]+1={} and the reversed sub-matls list is = {}\n".format(ind,keys))
            tmp[-1].geom.matls[-len(keys)-1:-1]=keys
            module_logger.debug("The reversed matls list is = {}\n".format(tmp[-1].geom.matls))
            
            # Update Cell Densities
            j=tmp[-1].fixed_mats
            for c in tmp[-1].geom.cells:
                if c.comment=='horiz' or c.comment=='vert':
                    c.d=mats[tmp[-1].geom.matls[j]].density
                    c.m=j+1
                    j+=1
        else:
            tmp.pop()
            
    return tmp

#---------------------------------------------------------------------------------------#
def Two_opt(x,S):
    """
    Implement 2_opt by reordering layers for top parents.
   
    Parameters
    ==========
    x : list of parent objects
        The current parents representing system designs
    S : Object    
        An object representing the settings for the optimization algorithm
   
    Optional
    ========   
   
    Returns
    =======
    tmp : list of parent objects
        The proposed parents representing new system designs
    """    
    tmp=[]
    cell_ids=[]
    
    for i in range(0,int(S.fe*S.p)):
        # Make a local copy
        tmp.append(cp.deepcopy(x[i]))
        # Compile list of horizontal cells
        for c in range(0, len(tmp[i].geom.cells)):
            if tmp[i].geom.cells[c].comment=='horiz':
                cell_ids.append(c)
        module_logger.debug("The horizontal cells are at positions = {}\n".format(cell_ids))      
        
        # Select random layer as starting point 
        rand=int(ceil(random()*(len(cell_ids)-3)))
        t_cell=cp.deepcopy(tmp[i].geom.cells[cell_ids[rand]]) 
        tmp[i].geom.cells[cell_ids[rand]]=cp.deepcopy(tmp[i].geom.cells[cell_ids[rand+1]]) 
        tmp[i].geom.cells[cell_ids[rand+1]]=t_cell
        module_logger.debug("Cell[{}] = {}\n".format(rand,tmp[i].geom.cells[cell_ids[rand]])) 
        module_logger.debug("Cell[{}] = {}\n".format(rand+1,tmp[i].geom.cells[cell_ids[rand+1]]))
        
        # Rename moved cells
        t_name=tmp[i].geom.cells[cell_ids[rand]].name
        tmp[i].geom.cells[cell_ids[rand]].name=tmp[i].geom.cells[cell_ids[rand+1]].name 
        tmp[i].geom.cells[cell_ids[rand+1]].name=t_name
        
        # Update the corresponding surface cards
        n_1=int(tmp[i].geom.cells[cell_ids[rand+1]].geom[1:4])
        n_2=int(tmp[i].geom.cells[cell_ids[rand]].geom[1:4])
        n_3=n_2+1
        for s in tmp[i].geom.surfaces:
            if s.name == n_1:
                z_1=s.d
            elif s.name == n_2:
                z_2=s.d
            elif s.name == n_3:
                z_3=s.d
        module_logger.debug("Old: {} = {}, {} = {}, {} = {}\n".format(n_1,z_1,n_2,z_2,n_3,z_3)) 
        del_1=(z_2-z_1)
        del_2=(z_3-z_2)
        z_2=z_1+del_2
        z_3=z_2+del_1
        module_logger.debug("New: {} = {}, {} = {}, {} = {}\n".format(n_1,z_1,n_2,z_2,n_3,z_3)) 
        for s in tmp[i].geom.surfaces:
            if s.name == n_2:
                s.d=z_2
            elif s.name == n_3:
                s.d=z_3
        
        # Update the cell decriptions
        t_geom=cp.deepcopy(tmp[i].geom.cells[cell_ids[rand]].geom) 
        tmp[i].geom.cells[cell_ids[rand]].geom=cp.deepcopy(tmp[i].geom.cells[cell_ids[rand+1]].geom) 
        tmp[i].geom.cells[cell_ids[rand+1]].geom=t_geom
        module_logger.debug("Cell[{}] = {}\n".format(rand,tmp[i].geom.cells[cell_ids[rand]])) 
        module_logger.debug("Cell[{}] = {}\n".format(rand+1,tmp[i].geom.cells[cell_ids[rand+1]]))
                  
    return tmp

#---------------------------------------------------------------------------------------#
def Crossover(x,S):
    """
    For each parent in top S.fe parents, N1, randomly select a parent, N2.  
    Randomly select an overall cell from N1 and copy into N2. Repeat s.pt times.
   
    Parameters
    ==========
    x : list of parent objects
        The current parents representing system designs
    S : Object    
        An object representing the settings for the optimization algorithm
   
    Optional
    ========   
   
    Returns
    =======
    tmp : list of parent objects
        The proposed parents representing new system designs
    """    
    tmp=[]
    cell_ids=[]
    used=[]
    
    for i in range(0,int(S.fe*S.p)):
        # Select random parent 
        rand=int(random()*S.p)
        while x[rand].ident==x[i].ident or rand in used:
            rand=int(random()*S.p)
            
        # Make a local copy
        tmp.append(cp.deepcopy(x[rand]))
        used.append(rand)
        
        # Compile list of possible cells to switch
        for c in range(0, len(x[i].geom.cells)):
            if x[i].geom.cells[c].comment=='horiz' or x[i].geom.cells[c].comment=='vert':
                cell_ids.append(c)
                
        # Select random cell and copy into new geometry
        rand=int(random()*len(cell_ids))
        tmp[i].geom.cells[cell_ids[rand]]=cp.deepcopy(x[i].geom.cells[cell_ids[rand]]) 
        module_logger.debug("The selected cell was cell[{}] = {}\n".format(cell_ids[rand],x[i].geom.cells[cell_ids[rand]])) 
        
        # Copy material into new geometry
        module_logger.debug("The top matls copied is = {}\n".format(x[i].geom.matls[x[i].geom.cells[cell_ids[rand]].m-1])) 
        module_logger.debug("The old matls list is = {}\n".format(tmp[i].geom.matls)) 
        tmp[i].geom.matls[x[i].geom.cells[cell_ids[rand]].m-1]=cp.deepcopy(x[i].geom.matls[x[i].geom.cells[cell_ids[rand]].m-1]) 
        module_logger.debug("The new matls list is = {}\n".format(tmp[i].geom.matls))  
        
        # Update the corresponding surface cards
        if x[i].geom.cells[cell_ids[rand]].comment=='horiz':
            n_1=int(x[i].geom.cells[cell_ids[rand]].geom[1:4])
            n_2=n_1+1
            for s in x[i].geom.surfaces:
                if s.name == n_1:
                    z_1=s.d
                elif s.name == n_2:
                    z_2=s.d
            module_logger.debug("Surface #{}={} and Surface #{}={}".format(n_1,z_1,n_2,z_2))
            for s in tmp[i].geom.surfaces:
                if s.name < n_1 and s.d > z_1:
                    s.d=z_1
                elif s.name == n_1:
                    s.d=z_1
                elif s.name == n_2:
                    s.d=z_2
                elif s.name > n_2 and s.d < z_2:
                    s.d=z_2
                module_logger.debug("Surface #{}={}".format(s.name,s.d))
        elif x[i].geom.cells[cell_ids[rand]].comment=='vert':
            n_1=int(x[i].geom.cells[cell_ids[rand]].geom[1:4])
            n_2=n_1+1  
            for j in range(0,len(x[i].geom.surfaces)):
                if x[i].geom.surfaces[j].name == n_1 or x[i].geom.surfaces[j].name == n_2:
                    tmp[i].geom.surfaces[j]=cp.deepcopy(x[i].geom.surfaces[j])
    return tmp

#---------------------------------------------------------------------------------------#
def Three_opt(x,S):
    """
    Perform for horizontal macrobodies if the number of cells is greater than 6.  
    Reorganizes the cells is all of the possible combinations for each parent.
   
    Parameters
    ==========
    x : list of parent objects
        The current parents representing system designs
    S : Object    
        An object representing the settings for the optimization algorithm
   
    Optional
    ========   
   
    Returns
    =======
    tmp : list of parent objects
        The proposed parents representing new system designs
    """    
    tmp=[]
    
    for i in range(0,int(S.p)):
        tmp.append(cp.deepcopy(x[i]))
        
        # Compile list of horizontal cell objects
        cells=[]
        for c in tmp[i].geom.cells:
            if c.comment=='horiz':
                cells.append(cp.deepcopy(c))     
        
        # Create a copy of corresponding surface objects
        surfs=[]
        for s in tmp[i].geom.surfaces:
            if s.c[0:5]=="horiz":
                surfs.append(cp.copy(s))
        module_logger.debug("The old geom.surfaces are: {}\n".format(tmp[i].geom.surfaces))  
        module_logger.debug("The old surfaces are: {}\n".format(surfs))  
                
        # Calculate the delta values for each cell
        dels=[]
        for j in range(1,len(surfs)):
            dels.append(surfs[j].d-surfs[j-1].d)
        module_logger.debug("The delta values are are: {}\n".format(dels))  
        
        # Select random layer as starting point for 'special A operator'
        rand=int(round(random()*(len(cells)-6)))
        module_logger.debug("The starting point is {} for the cells: {}\n".format(rand,cells)) 
        
        # Modify the original order
        p=random()
        new_cells=[]
        if p <=0.5:
            for a in range(0,rand+1):
                new_cells.append(cells[a])
                surfs[a+1].d=surfs[a].d+dels[a]
            new_cells.append(cells[rand+3])
            surfs[rand+2].d=surfs[rand+1].d+dels[rand+3]
            new_cells.append(cells[rand+4])
            surfs[rand+3].d=surfs[rand+2].d+dels[rand+4]
            new_cells.append(cells[rand+1])
            surfs[rand+4].d=surfs[rand+3].d+dels[rand+1]
            new_cells.append(cells[rand+2])
            surfs[rand+5].d=surfs[rand+4].d+dels[rand+2]
            new_cells.append(cells[rand+5])
            surfs[rand+6].d=surfs[rand+5].d+dels[rand+5]
            for a in range(rand+6,len(cells)):
                new_cells.append(cells[a])
                surfs[a+1].d=surfs[a].d+dels[a]
        elif p <=1.0:
            for a in range(0,rand+1):
                new_cells.append(cells[a])
                surfs[a+1].d=surfs[a].d+dels[a]
            new_cells.append(cells[rand+2])
            surfs[rand+2].d=surfs[rand+1].d+dels[rand+2]
            new_cells.append(cells[rand+1])
            surfs[rand+3].d=surfs[rand+2].d+dels[rand+1]
            new_cells.append(cells[rand+4])
            surfs[rand+4].d=surfs[rand+3].d+dels[rand+4]
            new_cells.append(cells[rand+3]) 
            surfs[rand+5].d=surfs[rand+4].d+dels[rand+3] 
            new_cells.append(cells[rand+5]) 
            surfs[rand+6].d=surfs[rand+5].d+dels[rand+5]
            for a in range(rand+6,len(cells)):
                new_cells.append(cells[a])
                surfs[a+1].d=surfs[a].d+dels[a]
        else:
            module_logger.warning("The modification did not occur for p={} in 3-opt.".format(p))
        module_logger.debug("The new cells are: {}\n".format(new_cells))
        module_logger.debug("The new surfaces are: {}\n".format(surfs)) 
        
        # Copy new cells into geometry
        for j in range(0,len(tmp[i].geom.cells)):
            if tmp[i].geom.cells[j].comment=='horiz':
                tmp[i].geom.cells[j].m=new_cells[0].m
                tmp[i].geom.cells[j].d=new_cells[0].d
                new_cells=new_cells[1:]
        if len(new_cells)!=0:
            module_logger.error("The copy of cells in 3-opt failed. Remaining cells={}".format(new_cells))
        
        # Copy new surfaces into geometry
        for j in range(0,len(tmp[i].geom.surfaces)):
            if tmp[i].geom.surfaces[j].c[0:5]=="horiz":
                tmp[i].geom.surfaces[j]=cp.deepcopy(surfs[0])
                surfs=surfs[1:]
        if len(surfs)!=0:
            module_logger.error("The copy of surfaces in 3-opt failed. Remaining surfaces={}".format(surfs))
        module_logger.debug("The new geom.surfaces are: {}\n".format(tmp[i].geom.surfaces))  
        
    return tmp  

#---------------------------------------------------------------------------------------#
def Discard_Cells(x,mats,S):
    """
    Discard a cell from fd parents. Bias discard towards better solutions.  Only accept if the 
    discard improves the fitness.
   
    Parameters
    ==========
    x : list of parent objects
        The current parents representing system designs
    mats : dictionary of material objects
        A materials library containing all relevant nulcear data required to run radiation transport codes.
    S : Object    
        An object representing the settings for the optimization algorithm
   
    Optional
    ========   
   
    Returns
    =======
    tmp : list of parent objects
        The proposed parents representing new system designs
    """              
    used=[]
    tmp=[]
    
    #Determine weights based on fitness    
    weights=WeightedRandomGenerator(rankdata(range(0,len(x),1)))
    
    for i in range(0,int(S.fd*S.p),1):
        discard=next(weights)    #Linearly weighted random parent to discard; biased towards better solutions
        #discard=int(random()*S.p)
        while discard in used:
            discard=next(weights)
            #discard=int(random()*S.p)
        used.append(discard)
        
        # Make a local copy
        tmp.append(cp.deepcopy(x[discard]))
        
        # Compile list of possible cells to abandon
        cell_ids=[]
        for c in range(0, len(tmp[-1].geom.cells)):
            if tmp[-1].geom.cells[c].comment=='horiz' or tmp[-1].geom.cells[c].comment=='vert':
                cell_ids.append(c)
                
        # Select random cell
        rand=int(random()*len(cell_ids))
        module_logger.debug("The selected parent was parent # {} ranked #{}.  The chosen cell was #{}.\n".format(tmp[-1].ident, discard, cell_ids[rand])) 
        module_logger.debug("The cell details are: {}.\n".format(tmp[-1].geom.cells[cell_ids[rand]]))
        
        # Change materials to 'delete' cell; if last cell, change to vaccuum
        if rand == len(cell_ids)-1:
            tmp[-1].geom.cells[cell_ids[rand]].m=next((i for i, item in enumerate(tmp[-1].geom.matls) if item == tmp[-1].geom.matls[-1]), -1)
            tmp[-1].geom.cells[cell_ids[rand]].d=mats[tmp[-1].geom.matls[-1]].density
            module_logger.debug("The new material is # {}, {}, dens={}".format(tmp[-1].geom.cells[cell_ids[rand]].m, tmp[-1].geom.matls[-1],tmp[-1].geom.cells[cell_ids[rand]].d))
            
        elif tmp[-1].geom.cells[cell_ids[rand]].comment=='horiz':
            module_logger.debug("The old choosen cell material is # {}, {}, dens={}".format(tmp[-1].geom.cells[cell_ids[rand]].m, tmp[-1].geom.matls[tmp[-1].geom.cells[cell_ids[rand]].m-1], tmp[-1].geom.cells[cell_ids[rand]].d))
            tmp[-1].geom.cells[cell_ids[rand]].m=tmp[-1].geom.cells[cell_ids[rand+1]].m
            tmp[-1].geom.cells[cell_ids[rand]].d=tmp[-1].geom.cells[cell_ids[rand+1]].d
            module_logger.debug("The new material is # {}, {}, dens={}".format(tmp[-1].geom.cells[cell_ids[rand]].m, tmp[-1].geom.matls[tmp[-1].geom.cells[cell_ids[rand]].m-1], tmp[-1].geom.cells[cell_ids[rand]].d))
            
        elif tmp[-1].geom.cells[cell_ids[rand]].comment=='vert':
            n_1=int(tmp[-1].geom.cells[cell_ids[rand+1]].geom[1:4])
            for i in range(0,len(tmp[-1].geom.surfaces)):
                if tmp[-1].geom.surfaces[i].name == n_1:
                    module_logger.debug("Old surface #1 = {}.\n".format(tmp[-1].geom.surfaces[i])) 
                    module_logger.debug("Old surface #2 = {}.\n".format(tmp[-1].geom.surfaces[i+1]))
                    tmp[-1].geom.surfaces[i+1].hz=0.0001
                    tmp[-1].geom.surfaces[i].hz=0.0001
                    module_logger.debug("New surface #1 = {}.\n".format(tmp[-1].geom.surfaces[i])) 
                    module_logger.debug("New surface #2 = {}.\n".format(tmp[-1].geom.surfaces[i+1]))
        else:
            module_logger.error("The selected cell #{} was of the incorrect type.  Selected cell = {}\n".format(cell_ids[rand],tmp[-1].geom.cells[cell_ids[rand]]))         
                            
    return tmp

#---------------------------------------------------------------------------------------#
def Mutate(x,eta,S):
    """
    Mutate parent population and build new ones. Bias discard to poor solutions.
    
    Parameters
    ==========
    x : list of parent objects
        The current parents representing system designs
    eta : Object    
        An object representing the constraints for the eta design
    S : Object    
        An object representing the settings for the optimization algorithm
   
    Optional
    ========   
   
    Returns
    =======
    y : list of parent objects
        The proposed parents representing new system designs
    """ 
    tmp=[]
    old=[]
    y=cp.deepcopy(x)
    
    # Build Design vectors
    for i in range(S.p): 
        cur_d=[] 
        
        # Build design variable set from current parent
        # [foil_z, N_vert*(z, delz, r1, r2), N_horiz*(z)]
        prev_vert=''
        for s in x[i].geom.surfaces:
            if s.c=="NAS":
                module_logger.debug("Found NAS. VZ={} and Cell={}".format(s.vz,repr(s)))
                cur_d.append(s.vz)
                
            elif s.c[0:4]=="vert":
                module_logger.debug("Found {}. VZ={}, HZ={}, and r={} and Cell={}".format(s.c,s.vz,s.hz,s.r,repr(s)))
                if prev_vert==s.c:
                    cur_d.append(s.r)
                    prev_vert=s.c
                else:
                    cur_d.append(s.vz)
                    cur_d.append(s.hz)
                    cur_d.append(s.r)
                    prev_vert=s.c
                    
            elif s.c[0:7]=="horiz #":
                module_logger.debug("Found {}. delZ={} and Cell={}".format(s.c,s.d,repr(s)))
                if s.c=="horiz #1":
                    cur_d.append(s.d-(eta.tcc_dist+eta.t_ds))
                else:
                    cur_d.append(s.d-prev_z)
                prev_z=s.d
        old.append(cur_d)
        module_logger.debug("Design Variable set for parent #{} = {}\n".format(x[i].ident,cur_d))
        
    # Convert to numpy arrays  
    old=np.asarray(old)
    tmp=cp.copy(old)
    module_logger.debug("Initial designs ={}\n".format(tmp)) 
        
    #Discover (1-fd); K is a status vector to see if discovered
    K=np.random.rand(len(tmp),len(tmp[0,:]))>S.fd
        
    #Bias the discovery to the worst fitness solutions
    childn1=cp.copy(np.random.permutation(tmp))
    childn2=cp.copy(np.random.permutation(tmp))
    module_logger.debug("Permutation #1 ={}\n".format(childn1)) 
    module_logger.debug("Permutation #2 ={}\n".format(childn2)) 
    
    #New solution by biased/selective random walks
    r=np.random.rand()
    for j in range(0,len(tmp),1):
        lb=[]
        ub=[]
        
        # Mutate designs
        n=np.array((childn1[j,:]-childn2[j,:]))
        step_size=r*n
        tmp[j,:]=tmp[j,:]+step_size*K[j,:]
        
        # Calculate Foil_Z Boundaries
        lb.append(max(eta.tcc_dist+eta.t_w+2*eta.t_nas+sum(eta.t_nas_f)+0.203,eta.tcc_dist+eta.t_w+(eta.r_nas-eta.r_f+eta.t_w)/tan(radians(90-eta.theta))))
        ub.append(eta.snout_dist-eta.t_m-2*eta.t_nas-sum(eta.t_nas_f)-0.203)
        
        # Calculate Vertical Cell Boundaries (z, delz, r1, r2)
        for i in range(0,eta.max_vert):
            if tmp[j,i*4+3] > tmp[j,i*4+4]:
                t=tmp[j,i*4+4]
                tmp[j,i*4+4]=tmp[j,i*4+3]
                tmp[j,i*4+3]=t
            lb.append(eta.tcc_dist+eta.t_ds)
            lb.append(0.0000001)
            lb.append(0.0000001)
            if tmp[j,i*4+3] > 0.0:
                lb.append(tmp[j,i*4+3])  
            else:
                lb.append(0.0000001)   
            ub.append(eta.snout_dist-eta.t_c-0.00001)
            ub.append(eta.snout_dist-eta.t_c-tmp[j,i*4+1])
            ub.append(eta.r_o)
            ub.append(eta.r_o)
            
        # Calculate horizontal bounds
        lb.append(eta.tcc_dist+eta.t_ds)
        ub.append(eta.snout_dist-eta.t_c)
        delta=tmp[j,-eta.max_horiz]
        for i in range(-eta.max_horiz+1,0):
            lb.append(0.0000001)
            ub.append(eta.snout_dist-eta.t_c-delta)
            delta+=tmp[j,i]
        lb=np.array(lb)
        ub=np.array(ub)
                
        module_logger.debug("Lower Bounds ={}\n".format(lb)) 
        module_logger.debug("Upper Bounds ={}\n".format(ub))
        
        module_logger.debug("For parent #{}, pre-bounds ={}\n".format(x[j].ident,tmp[j])) 
        tmp[j]=Simple_Bounds(tmp[j],lb,ub,change_count=1)
        module_logger.debug("For parent #{}, post-bounds ={}\n".format(x[j].ident,tmp[j])) 
        
        # Update parents with new design set
        # [foil_z, N_vert*(z, delz, r1, r2), N_horiz*(z)]
        prev_vert=''
        prev_z=0.0
        new_d=cp.copy(tmp[j,:])
        for s in y[j].geom.surfaces:
            if s.c[-4:]=="foil" or s.c[-4:]=="TOAD":
                s.vz=s.vz+new_d[0]-old[j,0] 
            if s.c=="NAS":
                s.vz=new_d[0]
                new_d=new_d[1:]
                nas_vz=s.vz
            if s.c=="Holder":
                s.vz=nas_vz-eta.t_h
            if s.c=="Holder Fill":
                s.vz=nas_vz
                
            elif s.c[0:4]=="vert":
                if prev_vert==s.c:
                    s.vz=new_d[0]
                    s.hz=new_d[1]
                    s.r=new_d[3]
                    new_d=new_d[4:]
                    prev_vert=s.c
                else:
                    s.vz=new_d[0]
                    s.hz=new_d[1]
                    s.r=new_d[2]
                    prev_vert=s.c
                    
            elif s.c[0:7]=="horiz #":
                if new_d[0]<prev_z:
                    if prev_z+0.1<eta.snout_dist-eta.t_c:
                        s.d=prev_z+0.1
                    else:
                        s.d=prev_z
                else:
                    s.d=new_d[0]
                new_d=new_d[1:]
                prev_z=s.d
                
        module_logger.debug("For i={}, ident={}, the parent={}\n".format(i,y[j].ident,str(y[j].geom)))             

    return y
