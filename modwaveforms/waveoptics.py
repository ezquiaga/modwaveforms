import numpy as np
import mpmath as mp
from scipy.special import gamma
from bilby.gw.source import *
from .utils import point_mass as pm
from .utils.constants import *

""" Point lens waveform """

"Diffraction integral"
#Amplification function
def laguerre(n,a,x):
    """ Laguerre polynomial of degree n with parameter a at point x 
    
    Parameters
    ----------
    n : int
        Degree of the Laguerre polynomial
    a : float
        Parameter of the Laguerre polynomial
    x : float
        Point at which to evaluate the Laguerre polynomial

    Returns
    -------
    float
        Value of the Laguerre polynomial at x
    """

    nplaguerre = np.frompyfunc(mp.laguerre,3,1)
    return nplaguerre(n,a,x)
vlaguerre = np.vectorize(laguerre)

def F_pointlens(f,ML,y):
    """ Amplification function of a point lens

    Follows definition of Eq.16 https://arxiv.org/abs/2005.10702

    Parameters
    ----------
    fs : float
        Frequency of the gravitational wave in Hz
    ML : float
        Redshifted mass of the lens in solar masses
    y : float
        Dimensionless impact parameter of the lens

    Returns
    -------
    complex
        Amplification function of the point lens
    """

    w = 2.*np.pi * (4*TSUN*ML) * f

    lague = np.array(vlaguerre(-0.5j*w, 0, 0.5j*w*y**2.),dtype='complex')
    ampl = np.power(-0.5j,1.+0.5j*w)*np.power(w,1.+0.5j*w)*gamma(-0.5j*w)*lague
    ampl[w==0. + 0.j] = 1. #gamma of 0 is inf but F(0,b)=1.

    t_plus_pm = pm.t_delay_geom_plus(y)* pm.t_ref(ML) #+ image time delay in seconds
    ampl *= np.exp(-1j*t_plus_pm*2*np.pi*f) #eliminate global time delay of the + image

    # complex conjugate because of different sign convention in the Fourier transform of 
    # LALSuite (engineering) compared to the one used in wave optics (physics)
    return np.conjugate(ampl)
