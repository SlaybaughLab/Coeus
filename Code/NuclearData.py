#######################################################################################################
#
# Module : NuclearData.py
#
# Contains : Functions and methods used to build materials and handle nuclear data.  Reads cross-sections from
#            NJOY formatted library.  Data structure is made to mimic PARTISN input, but the material data
#            can be extended and adapted to MCNP or other transport codes. 
#
# Author : James Bevins
#
# Last Modified: 17Aug16
#
#######################################################################################################

import logging
module_logger = logging.getLogger('Coeus.NuclearData')

import pyne
import warnings
warnings.simplefilter("ignore")

from pyne.dbgen.materials_library import make_matslib, make_elements
from pyne.xs.data_source import SimpleDataSource
from pyne.data import atomic_mass
from math import log
from scipy.constants import N_A

import os
import sys

#---------------------------------------------------------------------------------------# 
def Build_Matlib(mat_path='/home/pyne-user/Dropbox/UCB/Research/ETAs/META-CODE/MCNP/pyne/eta_materials_compendium.csv', remove_gases=True, remove_liquids=True, remove_expensive=True):
    """
    Builds and initializes a library of elements and materials provided by user using PyNE material library 
    functions.  

    Parameters
    ==========

    Optional
    ========
    mat_path : str
        Absolute path to the location of the user supplied materials compendium
    remove_gases : boolean
        Allows the user to selectively remove gases from the elements library
    remove_liquids : boolean
        Allows the user to selectively remove liquids from the elements library
    remove_expensive : boolean
        Allows the user to selectively remove expensive materials from the elements library.  "Expensive" encompasses from a 
        materials hazard and cost perspective

    Returns
    =======
    mat_lib : dictionary of material objects
        A materials library containing all relevant nulcear data required to run radiation transport codes.  Isotopic densities 
        are in atoms/b-cm
        """
    
    # Test path for materials compendium file. Build materials library if file exists; only build element library if not
    if os.path.isfile(mat_path): 
        module_logger.info("Loading materials compendium file located at: {}\n".format(mat_path))
        mat_lib = make_matslib(mat_path)
    else:
        module_logger.info("No user supplier materials file located.  Will build elemental materials library instead.\n")
        mat_lib = make_elements()

    # Set Elemental Material Densities
    mat_lib=Set_Density(mat_lib)
    
    # Test to ensure all elements have material densities
    for i in mat_lib:
        if mat_lib[i].density==-1:
            sys.exit("ERROR: {} does not have a density.".format(i)) 
            
        
    # Trim the materials list down by removing engineered challenged materials
    mat_lib=Strip_Undesireables(mat_lib, remove_gases, remove_liquids, remove_expensive)
        
    return mat_lib

#---------------------------------------------------------------------------------------# 
def Set_Density(mat_lib):
    """
    Initialized the material density for the elemental library. 

    Parameters
    ==========
    mat_lib : dictionary of material objects
        A materials library containing all relevant nulcear data required to run radiation transport codes.  Isotopic densities 
        are in atoms/b-cm

    Optional
    ========

    Returns
    =======
    mat_lib : dictionary of material objects
        An updated materials library containing all relevant nulcear data required to run radiation transport codes.  
    """
    
    dens=dict(H=0.0000899, He=0.0001785, Li=0.535, Be=1.848, B=2.370, C=2.260, N=0.001251, O=0.001429, F=0.001696, Ne=0.000900,
              Na=0.968,Mg=1.738, Al=2.7, Si=2.330, P=1.823, S=1.960, Cl=0.003214, Ar=0.001784, K=0.856, Ca=1.550, Sc=2.985, 
              Ti=4.507, V=6.110, Cr=7.190, Mn=7.470, Fe=7.874, Co=8.9, Ni=8.908, Cu=8.960, Zn=7.140, Ga=5.904, Ge=5.323, 
              As=5.727, Se=4.819, Br=3.120, Kr=0.00375, Rb=1.532, Sr=2.630, Y=4.472, Zr=6.511, Nb=8.570, Mo=10.280,
              Ru=12.370, Rh=12.450, Pd=12.023, Ag=10.490, Cd=8.650, In=7.310, Sn=7.310, Sb=6.697, Te=6.240, I=4.940, Xe=0.0059,
              Cs=1.879, Ba=3.510, La=6.146, Ce=6.689, Pr=6.640, Nd=7.010, Sm=7.353, Eu=5.244, Gd=7.901, Tb=8.219, 
              Dy=8.551, Ho=8.795, Er=9.066, Tm=9.320, Yb=6.570, Lu=9.841, Hf=13.310, Ta=16.650, W=19.250, Re=21.020, Os=22.59,
              Ir=22.56, Pt=21.450, Au=19.3, Hg=13.534, Tl=11.850, Pb=11.340, Bi=9.780, Th=11.724, U=19.050, Pa=15.4)
    
    for n, d in dens.iteritems():
        if n in mat_lib:
            mat_lib[n].density=d 
        else:
            module_logger.warning("{} not found in the materials library.".format(n))
        
    return mat_lib

#---------------------------------------------------------------------------------------# 
def Strip_Undesireables(mat_lib, remove_gases, remove_liquids, remove_expensive):
    """
    Removes materials from library that don't work from an engineering, safety, or cost perspective. 

    Parameters
    ==========
    mat_lib : dictionary of material objects
        A materials library containing all relevant nulcear data required to run radiation transport codes.  Isotopic densities 
        are in atoms/b-cm
    remove_gases : boolean
        Allows the user to selectively remove gases from the materials library
    remove_liquids : boolean
        Allows the user to selectively remove liquids from the elements library
    remove_expensive : boolean
        Allows the user to selectively remove expensive materials from the materials library.  "Expensive" encompasses from a 
        materials hazard and cost perspective

    Optional
    ========

    Returns
    =======
    mat_lib : dictionary of material objects
        An updated materials library containing all relevant nulcear data required to run radiation transport codes.  
        Isotopic densities are in atoms/b-cm
    """
    
    if "Pa" in mat_lib:
        del mat_lib["Pa"] 
    else:
        module_logger.warning("'Pa' does not exist in the material library.")
    
    if remove_gases==True:
        lst=['H', 'He', 'N', 'O', 'F', 'Ne', 'Cl', 'Ar', 'Kr', 'Xe']
        for i in lst:
            if i in mat_lib:
                del mat_lib[i] 
            else:
                module_logger.warning("{} not found in the materials library.".format(i))
    
    if remove_liquids==True:
        lst=['Br','Hg','Cs']
        for i in lst:
            if i in mat_lib:
                del mat_lib[i] 
            else:
                module_logger.warning("{} not found in the materials library.".format(i))
    
    if remove_expensive==True:
        lst=['Sc','Ge','As','Se','Rb','Pd','Ag','Ho','Tm','Yb','Lu','Re','Os','Ir','Rh','Pt','Tl','Th', 'U']
        for i in lst:
            if i in mat_lib:
                del mat_lib[i] 
            else:
                module_logger.warning("{} not found in the materials library.".format(i))
        
    return mat_lib

#---------------------------------------------------------------------------------------# 
def Calc_Moderating_Ratio(mats):
    """
    Calculated and returns the moderating ratio for each material in a materials library. 
    Currently limited to 1 and 14 MeV and the EAS data source.  More sophisticated approaches are 
    possible but not implemented. The moderating ratio is calculated as:
    
         MR={1 - (A-1)^2/2A * ln[(A+1)(A-1)]} * Sima_s / Sigma_a
         
    The EAS data source does not have any cross-sections for Zn, Dy, or Er.

    Parameters
    ==========
    mats : dictionary of material objects
        A materials library containing all relevant nulcear data required to run radiation transport codes.  Isotopic densities 
        are in atoms/b-cm
    
    Optional
    ========

    Returns
    =======
    mr : list of Moderating_Ratio objects
        A list containing the 1 and 14 MeV moderating ratios for the input material library
    """
    # Initialize output and key list
    mr=[]
    key_lst=mats.keys()
    
    for i in key_lst:
        rho=mats[i].density
            
        #Find A for elements
        try:
            A=atomic_mass(mats[i].metadata['name'])
            
        #Find A for compounds
        except RuntimeError as r:
            A=0
            for k in mats[i].comp.keys():
                A+=atomic_mass(k)*mats[i].comp[k]
            
        # Calculate Lethargy
        xi=1 - (A-1)**2/(2*A) * log((A+1)/(A-1))
        
        # Get cross-section Data (Reaction #2=elastic scattering, #4=inel scattering, #16=n,2n, #27=absorption
        sds = SimpleDataSource()
        sig_el1,sig_inl1,sig_a1=0,0,0
        sig_el14,sig_inl14,sig_a14=0,0,0
        for k in mats[i].comp.keys():
            try:
                sig_el14 += sds.reaction(k, 2)[0]*mats[i].comp[k]   #Index 0 is 14 MeV, 1 is 1 MeV, 2 is thermal
                sig_inl14 += sds.reaction(k, 4)[0]*mats[i].comp[k] 
                sig_a14 += sds.reaction(k, 27)[0]*mats[i].comp[k]
                sig_el1 += sds.reaction(k, 2)[1]*mats[i].comp[k]   
                sig_inl1 += sds.reaction(k, 4)[1]*mats[i].comp[k] 
                sig_a1 += sds.reaction(k, 27)[1]*mats[i].comp[k]
            except TypeError as t:
                module_logger.warning("{}({}) cross-section not found in EAS data.".format(i,k))

        # Calculate macroscopic cross-sections
        sig_el14 = sig_el14*N_A*mats[i].density/A
        sig_inl14 = sig_inl14*N_A*mats[i].density/A
        sig_a14 = sig_a14*N_A*mats[i].density/A
        sig_el1 = sig_el1*N_A*mats[i].density/A
        sig_inl1 = sig_inl1*N_A*mats[i].density/A
        sig_a1 = sig_a1*N_A*mats[i].density/A
        
        # Calculate moderating ratio 
        try:
            mr.append(Moderating_Ratio(i, xi*(sig_el1+sig_inl1)/sig_a1, xi*(sig_el14+sig_inl14)/sig_a14))
        except ZeroDivisionError as z:
            module_logger.warning("Divide by zero error.  No absorption cross section for {} in EAS data.".format(i))
            mr.append(Moderating_Ratio(i,0.0,0.0))

    return mr

#-------------------------------------------------------------------------------------------------------------#  
class Moderating_Ratio:
    """
    Creates a moderating ratio object.  
   
    Attributes
    ==========
    name : str
        The material name 
    mr_1MeV : int
        The moderating ratio at 1 MeV
    mr_14MeV : int
        The moderating ratio at 14 MeV
        
    Returns
    =======
    None
    """
    
    def __init__(self, name, mr_1MeV=0, mr_14MeV=0):  
        
        assert isinstance(mr_1MeV, float)==True, 'mr_1MeV must be of type float. {} given.'.format(mr_1MeV)
        assert isinstance(mr_14MeV, float)==True, 'mr_14MeV must be of type float. {} given.'.format(mr_14MeV)
        
        self.name=name
        self.mr_1MeV=mr_1MeV
        self.mr_14MeV=mr_14MeV
               
    def __repr__(self):
        return "{0}, 1MeV={1}, 14MeV={2}"\
                .format(self.name, self.mr_1MeV, self.mr_14MeV)
    
    def __str__(self):
        header = ["Moderating Ratios Instance:"]
        header += ["Name:".format(self.name)]  
        header += ["1 MEV MR:".format(self.mr_1MeV)] 
        header += ["14 MEV MR:".format(self.mr_14MeV)]  
        header ="\n".join(header)+"\n"
        s = header 
        return s    
        

