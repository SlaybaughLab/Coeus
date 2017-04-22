
from Gnowee_Utilities import Gnowee_Settings

_FUNC_DICT = {"RelativeLeastSquares": RelativeLeastSquares, "LeastSquares": LeastSquares, "Uopt": Uopt}

## Calculate the U-optimality 
# @param c the candidate design
# @param d the objective design
# @return The u-optimality design based fitness
def Uopt(c,d):
    assert len(c)==len(d), "The length of the candidate and objective design must be equal in Uopt."  
   
    return np.sum(abs(d-c))

    
## Calculate the U-optimality 
# @param c the candidate design
# @param d the objective design
# @return The least-squares design based fitness
def LeastSquares(c,d):
    
    assert len(c)==len(d), "The length of the candidate and objective design must be equal in LeastSquares."  
   
    return np.sum((d-c)**2)

## Calculates the relative least squares.  Assumes a normalized input candidate and objective spectrum to simplify  
#    calculation (i.e. sum of bins should equal 1). Not valid for unnormalized spectra.  
# @param c the candidate design
# @param o the objective design
# @return The least-squares design based fitness
def RelativeLeastSquares(c,o):  
    assert len(c)==len(o), "The length of the candidate and objective design must be equal in RelativeLeastSquares."  
    rls=(o-c)**2/o
    
    # For bins with no tally results, project the fitness 
    loc=len(c)-6
    while c[loc]!=0.0:
        loc-=1
        if loc == -1:
            break
    rls[0:loc+1]=np.array([np.average(rls[loc+1:loc+4])]*(loc+1))
    return np.sum(rls) 


