import scipy.constants as constants
import numpy as np
import astropy
from astroquery.jplhorizons import Horizons

h = constants.h
kb = constants.k
c = constants.c

def brightness(T , d):
    '''
    Insert useful comment
    '''

    B = 2*h*(v**2)*k*T/(c**2)*4*np.pi*(d**2)

    return B

def planck_func(nu, T=2.7255, deriv=False):
    '''
    The Planck blackbody function

    Returns
    -------
    Bnu : float
      The Planck blackbody function for spectral radiance

    Keyword arguments
    ------------------
    T : float 
      (Default : 2.7255)
    deriv : bool
      (Default : False)

    '''

    x = h*nu/kb/T
    if deriv:
        Bnu = (2*kb*nu**2)/c**2 * ((x**2 * np.exp(x))/(np.exp(x)-1)**2)
    else:
        Bnu = (2*h*nu**3)/c**2 / (np.exp(x)-1)

    return Bnu

def rj_func(nu, T=1.0):
    '''
    The Rayleigh-Jeans limit of the Planck blackbody function

    Returns
    -------
    Bnu : float
      The RJ function for spectral radiance

    Keyword arguments
    ------------------
    T : float 
      (Default : 1.0)

    '''

    return 2*(nu**2)*kb*T/c**2

def Inu(nu,alpha=-1):

    return nu**alpha

def arcsecsq2srad(arcsecsq):
    '''
    Converts arcsec^2 to steradian
    '''

    srad = arcsecsq/(41252.96/4/np.pi)/3600/3600
    return srad

def srad2arcsecsq(srad):
    '''
    Converts steradian squared to arcsec^2
    '''

    arcsecsq = srad*(41252.96/4/np.pi)*3600*3600
    return arcsecsq

def arcminsq2srad(arcminsq):
    '''
    Converts arcmin^2 to steradian
    '''

    srad = arcminsq/(41252.96/4/np.pi)/3600
    return srad

def srad2arcminsq(srad):
    '''
    Converts steradian squared to arcmin^2
    '''

    arcminsq = srad*(41252.96/4/np.pi)*3600
    return arcminsq


def main():
    pass

if __name__ == '__main__':

    main()
