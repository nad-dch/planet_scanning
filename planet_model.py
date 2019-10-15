import scipy
import numpy as np
import astropy
from astroquery.jplhorizons import Horizons

h = scipy.constants.h

def brightness(T , d):
    '''
    Insert useful comment
    '''

    h = 6.2606014e-034
    B = 2*h*(v**2)*k*T/(c**2)*4*np.pi*(d**2)

    return B

def main():

    T = 2.7
    d = 10.

if __name__ == '__main__':

    main()
