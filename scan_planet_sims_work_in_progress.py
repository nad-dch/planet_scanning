import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import os
import numpy as np
import tod
import noise as nm
import beam_model as bm
import planet_model as pm
import example_jpl as jpl
opj = os.path.join

def gaus(x,a,x0,sigma):
	'''
	The Gaussian probability function

	Returns
	-------
	Array of Gaussian curve values

	'''
    return a*exp(-(x-x0)**2/(2*sigma**2))

def get_interpolated_theta(ra,dec):
	'''
	Interpolates ra and dec to correspond to 
	the given sample length

	Returns
	-------
	theta_int : float
		The angles for the given sample length

	'''

	x = np.arange(len(ra))
    xi = np.linspace(0, len(x), nsamp)
    rai = np.interp(xi, x, ra)
    deci = np.interp(xi, x, dec)
    theta_int = np.sqrt(ra**2 + dec**2)

    return(theta_int)

def estimate_bias(expected,measured):
	'''
	Returns the bias of the estimation
	'''

	bias = ((measured-expected)**2)/expected
	bias = np.sum(bias)
	return(bias)

def scan_planet(nsamp=1000, planet_id=5, nu=100e9, T=110, 
				beam_type='gaussian', noise_presence=True, 
				noise_type='onef', compute_bias=True):

	'''
	A function that outputs the signal tod from a planet scanning

	Arguments:
		nsamp : sample length
		planet_id : define the planet you want to scan
		nu : frequency
		T : blackbody temperature of the planet (110K for Jupyter)
		beam_type : gaussian, physical optics, elliptical gaussian
		noise_presence : bool
		noise_type : string 'gaussian','onef'
		compute_bias : bool. default value = True

	Returns:
	-------
		signal or fitted_values: Array of amplitude signal or 
									noisy signals fitted values
		bias : the bias of the estimation

	'''

	ra,dec = get_planet_timelines()
	theta = get_interpolated_theta(ra,dec)
    amplitude = pm.planck_func(nu,T)

    if beam_type='gaussian':
    	signal = amplitude * bm.gaussian_beam(theta)

    if noise_presence == False:
		return(signal)

    else:
    	#assure that you also change noise_type when
    	#turning noise_presence into True
    	if noise_type == None:
			print('you forgot to set noise_type')

    	if noise_type == 'white':
    		noise = nm.white(nsamp,sigma=0.01)

		else:
			noise = nm.noise_rel(nsamp, fsamp=5)

    	n_signal += noise

    	#fit the curve in the presence of noise
    	mean = sum(ra*n_signal)/n                   
		sigma = sum(n_signal*(ra-mean)**2)/n 
    	popt,pcov = curve_fit(gaus,ra,n_signal,p0=[1,mean,sigma])
    	fitted_values = gaus(ra,*popt)

    	#plt.plot(ra,n_signal,label='actual signal')
    	#plt.plot(ra,fitted_values,label='fitted signal')

    	if compute_bias == True:
    		bias = estimate_bias(signal,fitted_values)
    		#print(bias)

		return(fitted_values,bias)


def parametrizing_bias(sigma0,beam_width,ra,dec,time):
	'''
	Work in progress

	Compute bias for all different input features

	Construct a model that can predict bias behavior

	*use deep learning forecasting techniques
	'''

	bias_sigma = []
	sigma0 = np.array(np.linspace(0.01,1.0,100))
	for sigma in sigma0:
		fin, psdin = noise_psd(sigma0=1.0, fknee=0.020, alpha=1.5)
		fitted,bias = scan_planet()
		bias_sigma.append(bias)
	#plt.plot(sigma0,bias_sigma)


def exclude_fitting_bias():
	'''
	Work in progress

	add more noise each time considering the fitted curve
	as the true one (just signal) and then re-estimate the 
	bias. If bias is increasing then part of it is probably
	also caused by the fitting and not only noise systematics

	'''





