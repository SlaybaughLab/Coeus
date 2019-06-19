"""!
@file SamplingMethods.py
@package Coeus

@defgroup SamplingMethods SamplingMethods

@brief Different methods to perform phase space sampling and random walks and
visualization tools.

@author James Bevins

@date 16June19
"""

import logging
import math
import argparse
import random

import numpy as np

import scipy
from scipy import integrate
from Utilities import Switch
from DOE import lhs

module_logger = logging.getLogger('Coeus.SamplingMethods')

#-----------------------------------------------------------------------------#
def Initial_Samples(lb, ub, method, n=25):
    """!
    Generate a set of samples in a given phase space.

    @param lb: \e array \n
        The lower bounds of the design variable(s). \n
    @param ub: \e array \n
        The upper bounds of the design variable(s). \n
    @param method: \e string \n
        String representing the chosen sampling method. \n
    @param n: \e integer \n
        The number of samples to be generated.  Ignored for nolh algorithms.
        (Default=25). \n

    @return \e array: The list of coordinates for the sampled phase space. \n
    """

    assert len(ub) == len(lb), 'Boundaries best have different # of design ' \
                               'variables in Initial_Samples function.'
    if method != "random":
        assert len(ub) >= 2 and len(ub) <= 29, 'The Phase space dimensions ' \
                                               'are outside of the bounds ' \
                                               'for Initial_Samples.'
    assert method in ['random', 'nolh', 'nolh-rp', 'nolh-cdr', 'lhc'], \
                    'An invalid method was specified for the initial sampling.'

    for case in Switch(method):
        # Random draw from uniform distribution
        if case('random'):
            s = np.zeros((n, len(lb)))
            for i in range(0, n, 1):
                s[i, :] = (lb+(ub-lb)*np.random.rand(len(lb)))
            break

        # Standard nearly-orthoganal latin hypercube (LHC) call
        if case('nolh'):
            dim = len(ub)
            m, q, r = params(dim)
            conf = range(q)
            if r != 0:
                remove = range(dim - r, dim)
                nolh = NOLH(conf, remove)
            else:
                nolh = NOLH(conf)
            s = np.array([list(lb+(ub-lb)*nolh[i, :]) for i in range(
                len(nolh[:, 0]))])
            break

        # Nearly-orthoganal LHC with random permutation for removed colummns
        if case('nolh-rp'):
            dim = len(ub)
            m, q, r = params(dim)
            conf = random.sample(range(q), q)
            if r != 0:
                remove = random.sample(range(q-1), r)
                nolh = NOLH(conf, remove)
            else:
                nolh = NOLH(conf)
            s = np.array([list(lb+(ub-lb)*nolh[i, :]) for i in range(len(
                nolh[:, 0]))])
            break

        # Nearly-orthoganal LHC with Cioppa and De Rainville permutations
        if case('nolh-cdr'):
            dim = len(ub)
            m, q, r = params(dim)
            (conf, remove) = Get_CDR_Permutations(len(ub))
            if remove != []:
                nolh = NOLH(conf, remove)
            else:
                nolh = NOLH(conf)
            s = np.array([list(lb+(ub-lb)*nolh[i, :]) for i in range(len(
                nolh[:, 0]))])
            break

        # Latin hypercube sampling
        if case('lhc'):
             #Alt valid criterion are 'corr','center','maximum','centermaximum'
            tmp = lhs(len(lb), samples=n, criterion="center")
            s = np.array([list(lb+(ub-lb)*tmp[i, :]) for i in range(len(
                tmp[:, 0]))])
            break

        # Catch all
        if case():
            module_logger.warning('Somehow you evaded my assert statement - '
                                  'good job!  However, you still need to use a'
                                  ' valid method string.')

        module_logger.debug('Initial Samples: {}'.format(s))
    return s

#-----------------------------------------------------------------------------#
def Levy_Function(bins, alpha=1.5, gamma=1):
    """!
    Generate the levy function.

    @param bins: \e array \n
        The bin values used to generate the disribution. \n
    @param alpha: \e float \n
        Levy exponent - defines the index of the distribution and controls
        scale properties of the stochastic process (Default: 1.5). \n
    @param gamma: \e scalar \n
        Gamma - Scale unit of process for Levy flights (Default: 1). \n

    @return \e array: Array representing the levy flights for each member. \n
    """

    def _integrand(x, a, b, g):
        return np.exp(-g*x**(a))*np.cos(x*b)

    assert len(bins) >= 1, 'The length of the bins must be positive and at ' \
                            'least one.'
    assert 0.3 < alpha < 1.99, 'Valid range for alpha is [0.3:1.99].'
    assert gamma >= 0, 'Gamma must be positive'

    levy = np.zeros_like(bins)
    for i in range(0, len(bins)):
        levy[i] = 1.0/math.pi*integrate.quad(_integrand, 0, np.inf,
                                          args=(alpha, bins[i], gamma))[0]
    return levy

#-----------------------------------------------------------------------------#
def Levy(nc, nr=0, alpha=1.5, gamma=1, n=1):
    """!
     Sample the Levy distribution to generate Levy flights steps using the
     Mantegna Algorithm: "Fast, accurate algorithm for numerical simulation
     of Levy stable stochastic processes"

    @param nc: \e integer \n
        The number of columns of Levy values for the return array. \n
    @param nr: \e integer \n
        The number of rows of Levy values for the return array
        (Default: 0). \n
    @param alpha: \e float \n
        Levy exponent - defines the index of the distribution and controls
        scale properties of the stochastic process (Default: 1.5). \n
    @param gamma: \e scalar \n
        Gamma - Scale unit of process for Levy flights (Default: 1). \n

    @return \e array: Array representing the levy flights for each member. \n
    """

    assert 0.3 < alpha < 1.99, 'Valid range for alpha is [0.3:1.99].'
    assert gamma >= 0, 'Gamma must be positive'
    assert n >= 1, 'n Must be positive'

    invalpha = 1./alpha
    sigx = ((scipy.special.gamma(1.+alpha)*np.sin(np.pi*alpha/2.)) / \
            (scipy.special.gamma((1.+alpha)/2)* \
             alpha*2.**((alpha-1.)/2.)))**invalpha
    if nr > 0:
        v = sigx*np.random.randn(n, nr, nc)/ \
            (abs(np.random.randn(n, nr, nc))**invalpha)
    else:
        v = sigx*np.random.randn(n, nc)/(abs(np.random.randn(n, nc))**invalpha)

    kappa = (alpha*scipy.special.gamma((alpha+1.)/(2.*alpha)))/ \
            scipy.special.gamma(invalpha)* \
            ((alpha*scipy.special.gamma((alpha+1.)/2.))/ \
            (scipy.special.gamma(1.+alpha)*np.sin(np.pi*alpha/2.)))**invalpha
    p = [-17.7767, 113.3855, -281.5879, 337.5439, -193.5494, 44.8754]
    c = np.polyval(p, alpha)
    w = ((kappa-1.)*np.exp(-abs(v)/c)+1.)*v

    if n > 1:
        z = (1/n**invalpha)*sum(w)
    else:
        z = w

    z = gamma**invalpha*z
    if nr > 0:
        z = z.reshape(nr, nc)
    else:
        z = z.reshape(nc)

    module_logger.debug('In Levy flight algorithm: \n 1/alpha: {}\n X '
                        'Standard Deviation: {}\n K(alpha): {}\n C(alpha):  '
                        '{}'.format(invalpha, sigx, kappa, c))

    return z

#-----------------------------------------------------------------------------#
def TLF(alpha=1.5, gamma=1., numSamp=1, cutPoint=20.):
    """!
    Generates and samples from a truncated levy flight distribution. Produces a
    Levy Equivalent distribution on the interval (0,1).

    @param alpha: \e float \n
        Levy exponent - defines the index of the distribution and controls
        scale properties of the stochastic process (Default: 1.5). \n
    @param gamma: \e scalar \n
        Gamma - Scale unit of process for Levy flights (Default: 1). \n
    @param numSamp: \e integer \n
        Number of Levy flights to sample (Default: 1). \n
    @param cutPoint: \e integer \n
        Point at which to cut sampled Levy values and resample
        (Default: 20). \n

    @return \e array: Array representing the levy flights on the interval
        (0,1). \n
    """

    # Draw numSamp samples from the Levy distribution
    levy = abs(Levy(1, numSamp, alpha, gamma)/cutPoint).reshape(numSamp)

    module_logger.debug('Prior to resampling, the maximum sampled value is: {}'
                        '. {} of the samples are above the cut point.'.format(
                        np.max(levy), (levy > 1.0).sum()/numSamp))

    # Resample values above the range (0,1)
    for i in range(len(levy)):
        while levy[i] > 1:
            levy[i] = abs(Levy(1, 1, alpha, gamma)/cutPoint).reshape(1)

    return levy

#-----------------------------------------------------------------------------#
def NOLH(conf, remove=None):
    """!
    This library allows to generate Nearly Orthogonal Latin Hypercubes (NOLH)
    according to Cioppa (2007) and De Rainville et al. (2012) and reference
    therein.

    https://pypi.python.org/pypi/pynolh

    Constructs a Nearly Orthogonal Latin Hypercube (NOLH) of order *m* from a
    configuration vector *conf*. The configuration vector may contain either
    the numbers in $[0 q-1]$ or $[1 q]$ where $q = 2^{m-1}$. The columns to be
    *removed* are also in $[0 d-1]$ or $[1 d]$ where $d = m + \binom{m-1}{2}$
    is the NOLH dimensionality.

    The whole library is incorporated here with minimal modification for
    commonality and consolidation of methods.

    @param conf: \e array \n
        Configuration vector. \n
    @param remove: \e array \n
        Array containing the indexes of the colummns to be removed from conf
        vetor (Default: NONE). \n

    @return \e array: Array containing nearly orthogonal latin hypercube
        sampling. \n
    """

    I = np.identity(2, dtype=int)
    R = np.array(((0, 1), (1, 0)), dtype=int)

    if 0 in conf:
        conf = np.array(conf) + 1

        if remove is not None:
            remove = np.array(remove) + 1

    q = len(conf)
    m = math.log(q, 2) + 1
    s = m + (math.factorial(m - 1) / (2 * math.factorial(m - 3)))

    # Factorial checks if m is an integer
    m = int(m)

    A = np.zeros((q, q, m - 1), dtype=int)
    for i in range(1, m):
        Ai = 1
        for j in range(1, m):
            if j < m - i:
                Ai = np.kron(Ai, I)
            else:
                Ai = np.kron(Ai, R)

        A[:, :, i-1] = Ai

    M = np.zeros((q, s), dtype=int)
    M[:, 0] = conf

    col = 1
    for i in range(0, m - 1):
        for j in range(i + 1, m):
            if i == 0:
                M[:, col] = np.dot(A[:, :, j-1], conf)
            else:
                M[:, col] = np.dot(A[:, :, i-1], np.dot(A[:, :, j-1], conf))
            col += 1

    S = np.ones((q, s), dtype=int)
    v = 1
    for i in range(1, m):
        for j in range(0, q):
            if j % 2**(i-1) == 0:
                v *= -1
            S[j, i] = v

    col = m
    for i in range(1, m - 1):
        for j in range(i + 1, m):
            S[:, col] = S[:, i] * S[:, j]
            col += 1

    T = M * S

    keep = np.ones(s, dtype=bool)
    if remove is not None:
        keep[np.array(remove) - 1] = [False] * len(remove)

    return (np.concatenate((T, np.zeros((1, s)), -T), axis=0)[:, keep] + q) / \
           (2.0 * q)

#-----------------------------------------------------------------------------#
def params(dim):
    """!
    Returns the NOLH order $m$, the required configuration length $q$
    and the number of columns to remove to obtain the desired dimensionality.

    @param dim: \e integer \n
        The dimension of the phase space. \n

    @return \e array: Array containing nearly orthogonal latin hypercube
        sampling. \n
    """

    m = 3

    # Original version has three here, but this failed each time the # of
    # samples required switched (ie at dim=3,7,11,etc)
    s = 1
    q = 2**(m-1)

    while s < dim:
        m += 1
        s = m + math.factorial(m - 1) / (2 * math.factorial(m - 3))
        q = 2**(m-1)

    return m, q, s - dim

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=("Compute a Nearly "
        "Orthogonal Latin hypercube from a configuration vector."))
    parser.add_argument("conf", metavar="C", type=int, nargs="+",
        help="The configuration vector given as a list N1 N2 ... Nm")
    parser.add_argument("-r", "--remove", metavar="R", type=int, nargs="+",
        help="Columns to remove as a list N1 N2 ... Nm")

    args = parser.parse_args()
    module_logger.debug((NOLH(conf=args.conf, remove=args.remove)))

#-----------------------------------------------------------------------------#
def Get_CDR_Permutations(dim):
    """!
    Generate a set of samples in a given phase space.

    @param dim: \e integer \n
        The dimension of the phase space. \n

    @return \e array: Configuration vector. \n
    @return \e array: Array containing the indexes of the colummns to be
        removed from conf vector. \n
    """

    assert 2 <= dim <= 29, 'The Phase space dimensions are outside of the ' \
                                'bounds for CDR Permutations.'

    # Permutation and columns to remove given by Cioppa
    C_CONF = {
        2 : ([1, 2, 8, 4, 5, 6, 7, 3], [1, 3, 4, 6, 7]),
        3 : ([1, 2, 8, 4, 5, 6, 7, 3], [1, 2, 3, 6]),
        4 : ([1, 2, 8, 4, 5, 6, 7, 3], [1, 3, 6]),
        5 : ([1, 2, 8, 4, 5, 6, 7, 3], [1, 6]),
        6 : ([1, 2, 8, 4, 5, 6, 7, 3], [1]),
        7 : ([1, 2, 8, 4, 5, 6, 7, 3], [])
    }

    # Permutation and columns to remove given by De Rainville et al.
    EA_CONF = {
        8  : ([4, 14, 1, 2, 16, 13, 5, 8, 12, 9, 6, 7, 11, 3, 15, 10],
              [1, 3, 10]),
        9  : ([4, 14, 1, 2, 16, 13, 5, 8, 12, 9, 6, 7, 11, 3, 15, 10],
              [6, 10]),
        10 : ([4, 14, 1, 2, 16, 13, 5, 8, 12, 9, 6, 7, 11, 3, 15, 10], [10]),
        11 : ([4, 14, 1, 2, 16, 13, 5, 8, 12, 9, 6, 7, 11, 3, 15, 10], []),

        12 : ([5, 13, 19, 23, 28, 10, 12, 32, 17, 2, 30, 15, 6, 31, 21, 8, 24,
               29, 9, 14, 11, 22, 18, 25, 3, 1, 20, 7, 27, 16, 26, 4],
              [2, 4, 5, 11]),
        13 : ([5, 13, 19, 23, 28, 10, 12, 32, 17, 2, 30, 15, 6, 31, 21, 8, 24,
               29, 9, 14, 11, 22, 18, 25, 3, 1, 20, 7, 27, 16, 26, 4],
               [3, 6, 14]),
        14 : ([5, 13, 19, 23, 28, 10, 12, 32, 17, 2, 30, 15, 6, 31, 21, 8, 24,
               29, 9, 14, 11, 22, 18, 25, 3, 1, 20, 7, 27, 16, 26, 4], [4, 5]),
        15 : ([5, 13, 19, 23, 28, 10, 12, 32, 17, 2, 30, 15, 6, 31, 21, 8, 24,
               29, 9, 14, 11, 22, 18, 25, 3, 1, 20, 7, 27, 16, 26, 4], [6]),
        16 : ([5, 13, 19, 23, 28, 10, 12, 32, 17, 2, 30, 15, 6, 31, 21, 8, 24,
               29, 9, 14, 11, 22, 18, 25, 3, 1, 20, 7, 27, 16, 26, 4], []),

        17 : ([7, 8, 51, 3, 40, 44, 29, 19, 61, 43, 26, 48, 20, 52, 4, 49, 2,
               57, 31, 30, 24, 23, 56, 50, 18, 59, 63, 37, 38, 21, 54, 9, 46,
               27, 36, 1, 10, 42, 13, 55, 15, 25, 22, 45, 41, 39, 53, 34, 6, 5,
               2, 58, 16, 28, 64, 14, 47, 33, 12, 35, 62, 17, 11, 60],
               [8, 11, 12, 14, 17]),
        18 : ([7, 8, 51, 3, 40, 44, 29, 19, 61, 43, 26, 48, 20, 52, 4, 49, 2,
               57, 31, 30, 24, 23, 56, 50, 18, 59, 63, 37, 38, 21, 54, 9, 46,
               27, 36, 1, 10, 42, 13, 55, 15, 25, 22, 45, 41, 39, 53, 34, 6, 5,
               2, 58, 16, 28, 64, 14, 47, 33, 12, 35, 62, 17, 11, 60],
               [8, 11, 12, 17]),
        19 : ([7, 8, 51, 3, 40, 44, 29, 19, 61, 43, 26, 48, 20, 52, 4, 49, 2,
               57, 31, 30, 24, 23, 56, 50, 18, 59, 63, 37, 38, 21, 54, 9, 46,
               27, 36, 1, 10, 42, 13, 55, 15, 25, 22, 45, 41, 39, 53, 34, 6, 5,
               2, 58, 16, 28, 64, 14, 47, 33, 12, 35, 62, 17, 11, 60],
               [10, 15, 22]),
        20 : ([7, 8, 51, 3, 40, 44, 29, 19, 61, 43, 26, 48, 20, 52, 4, 49, 2,
               57, 31, 30, 24, 23, 56, 50, 18, 59, 63, 37, 38, 21, 54, 9, 46,
               27, 36, 1, 10, 42, 13, 55, 15, 25, 22, 45, 41, 39, 53, 34, 6, 5,
               2, 58, 16, 28, 64, 14, 47, 33, 12, 35, 62, 17, 11, 60],
               [8, 12]),
        21 : ([7, 8, 51, 3, 40, 44, 29, 19, 61, 43, 26, 48, 20, 52, 4, 49, 2,
               57, 31, 30, 24, 23, 56, 50, 18, 59, 63, 37, 38, 21, 54, 9, 46,
               27, 36, 1, 10, 42, 13, 55, 15, 25, 22, 45, 41, 39, 53, 34, 6, 5,
               2, 58, 16, 28, 64, 14, 47, 33, 12, 35, 62, 17, 11, 60], [15]),
        22 : ([7, 8, 51, 3, 40, 44, 29, 19, 61, 43, 26, 48, 20, 52, 4, 49, 2,
               57, 31, 30, 24, 23, 56, 50, 18, 59, 63, 37, 38, 21, 54, 9, 46,
               27, 36, 1, 10, 42, 13, 55, 15, 25, 22, 45, 41, 39, 53, 34, 6, 5,
               2, 58, 16, 28, 64, 14, 47, 33, 12, 35, 62, 17, 11, 60], []),

        23 : ([9, 108, 39, 107, 62, 86, 110, 119, 46, 43, 103, 71, 123, 91, 10,
               13, 126, 63, 83, 47, 100, 54, 23, 16, 124, 45, 27, 4, 93, 74,
               76, 90, 30, 81, 77, 53, 116, 49, 104, 6, 70, 82, 26, 118, 55,
               79, 32, 109, 57, 31, 22, 101, 44, 87, 121, 7, 37, 56, 89, 115,
               25, 92, 85, 20, 58, 52, 3, 11, 106, 17, 117, 38, 78, 28, 59, 96,
               18, 97, 50, 114, 112, 60, 84, 1, 12, 61, 98, 128, 14, 42, 64,
               105, 68, 75, 111, 34, 141, 65, 99, 2, 19, 33, 35, 94, 51, 122,
               127, 36, 125, 80, 73, 8, 24, 21, 88, 48, 69, 66, 40, 15, 29,
               113, 72, 5, 95, 120, 6, 102], [18, 20, 21, 24, 27, 29]),
        24 : ([9, 108, 39, 107, 62, 86, 110, 119, 46, 43, 103, 71, 123, 91, 10,
               13, 126, 63, 83, 47, 100, 54, 23, 16, 124, 45, 27, 4, 93, 74,
               76, 90, 30, 81, 77, 53, 116, 49, 104, 6, 70, 82, 26, 118, 55,
               79, 32, 109, 57, 31, 22, 101, 44, 87, 121, 7, 37, 56, 89, 115,
               25, 92, 85, 20, 58, 52, 3, 11, 106, 17, 117, 38, 78, 28, 59, 96,
               18, 97, 50, 114, 112, 60, 84, 1, 12, 61, 98, 128, 14, 42, 64,
               105, 68, 75, 111, 34, 141, 65, 99, 2, 19, 33, 35, 94, 51, 122,
               127, 36, 125, 80, 73, 8, 24, 21, 88, 48, 69, 66, 40, 15, 29,
               113, 72, 5, 95, 120, 6, 102], [4, 15, 18, 24, 27]),
        25 : ([9, 108, 39, 107, 62, 86, 110, 119, 46, 43, 103, 71, 123, 91, 10,
               13, 126, 63, 83, 47, 100, 54, 23, 16, 124, 45, 27, 4, 93, 74,
               76, 90, 30, 81, 77, 53, 116, 49, 104, 6, 70, 82, 26, 118, 55,
               79, 32, 109, 57, 31, 22, 101, 44, 87, 121, 7, 37, 56, 89, 115,
               25, 92, 85, 20, 58, 52, 3, 11, 106, 17, 117, 38, 78, 28, 59, 96,
               18, 97, 50, 114, 112, 60, 84, 1, 12, 61, 98, 128, 14, 42, 64,
               105, 68, 75, 111, 34, 141, 65, 99, 2, 19, 33, 35, 94, 51, 122,
               127, 36, 125, 80, 73, 8, 24, 21, 88, 48, 69, 66, 40, 15, 29,
               113, 72, 5, 95, 120, 6, 102], [21, 26, 27, 29]),
        26 : ([9, 108, 39, 107, 62, 86, 110, 119, 46, 43, 103, 71, 123, 91, 10,
               13, 126, 63, 83, 47, 100, 54, 23, 16, 124, 45, 27, 4, 93, 74,
               76, 90, 30, 81, 77, 53, 116, 49, 104, 6, 70, 82, 26, 118, 55,
               79, 32, 109, 57, 31, 22, 101, 44, 87, 121, 7, 37, 56, 89, 115,
               25, 92, 85, 20, 58, 52, 3, 11, 106, 17, 117, 38, 78, 28, 59, 96,
               18, 97, 50, 114, 112, 60, 84, 1, 12, 61, 98, 128, 14, 42, 64,
               105, 68, 75, 111, 34, 141, 65, 99, 2, 19, 33, 35, 94, 51, 122,
               127, 36, 125, 80, 73, 8, 24, 21, 88, 48, 69, 66, 40, 15, 29,
               113, 72, 5, 95, 120, 6, 102], [26, 27, 29]),
        27 : ([9, 108, 39, 107, 62, 86, 110, 119, 46, 43, 103, 71, 123, 91, 10,
               13, 126, 63, 83, 47, 100, 54, 23, 16, 124, 45, 27, 4, 93, 74,
               76, 90, 30, 81, 77, 53, 116, 49, 104, 6, 70, 82, 26, 118, 55,
               79, 32, 109, 57, 31, 22, 101, 44, 87, 121, 7, 37, 56, 89, 115,
               85, 20, 58, 52, 3, 11, 106, 17, 117, 38, 78, 28, 59, 96, 18, 97,
               50, 114, 112, 60, 84, 1, 12, 61, 98, 128, 14, 42, 64, 105, 68,
               75, 111, 34, 141, 65, 99, 2, 19, 33, 35, 94, 51, 122, 127, 36,
               125, 80, 73, 8, 24, 21, 88, 48, 69, 66, 40, 15, 29, 113, 72, 5,
               95, 120, 6, 102], [27, 29]),
        28 : ([9, 108, 39, 107, 62, 86, 110, 119, 46, 43, 103, 71, 123, 91, 10,
               13, 126, 63, 83, 47, 100, 54, 23, 16, 124, 45, 27, 4, 93, 74,
               76, 90, 30, 81, 77, 53, 116, 49, 104, 6, 70, 82, 26, 118, 55,
               79, 32, 109, 57, 31, 22, 101, 44, 87, 121, 7, 37, 56, 89, 115,
               25, 92, 85, 20, 58, 52, 3, 11, 106, 17, 117, 38, 78, 28, 59, 96,
               18, 97, 50, 114, 112, 60, 84, 1, 12, 61, 98, 128, 14, 42, 64,
               105, 68, 75, 111, 34, 141, 65, 99, 2, 19, 33, 35, 94, 51, 122,
               127, 36, 125, 80, 73, 8, 24, 21, 88, 48, 69, 66, 40, 15, 29,
               113, 72, 5, 95, 120, 6, 102], [20]),
        29 : ([9, 108, 39, 107, 62, 86, 110, 119, 46, 43, 103, 71, 123, 91, 10,
               13, 126, 63, 83, 47, 100, 54, 23, 16, 124, 45, 27, 4, 93, 74,
               76, 90, 30, 81, 77, 53, 116, 49, 104, 6, 70, 82, 26, 118, 55,
               79, 32, 109, 57, 31, 22, 101, 44, 87, 121, 7, 37, 56, 89, 115,
               25, 92, 85, 20, 58, 52, 3, 11, 106, 17, 117, 38, 78, 28, 59, 96,
               18, 97, 50, 114, 112, 60, 84, 1, 12, 61, 98, 128, 14, 42, 64,
               105, 68, 75, 111, 34, 141, 65, 99, 2, 19, 33, 35, 94, 51, 122,
               127, 36, 125, 80, 73, 8, 24, 21, 88, 48, 69, 66, 40, 15, 29,
               113, 72, 5, 95, 120, 6, 102], [])
    }

    # Create dictionary
    CONF = dict()
    CONF.update(C_CONF)
    CONF.update(EA_CONF)

    return CONF[dim][0], CONF[dim][1]
