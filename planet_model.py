import numpy as np
import astropy
from astroquery.jplhorizons import Horizons
import os

def brightness(T , d):

	h = 6.2606014e-034
	#v= 
	B = 2*h*(v**2)*k*T/(c**2)*4*np.pi*(d**2)

	return B



	
