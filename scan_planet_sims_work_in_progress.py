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
	return a*np.exp(-(x-x0)**2/(2*sigma**2))

def get_interpolated_theta(ra,dec,nsamp):
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
	theta_int = np.sqrt(rai**2 + deci**2)

	return(rai,deci,theta_int)

def estimate_bias(expected,measured):
	'''
	Returns the bias of the estimation
	'''

	bias = ((measured-expected)**2)/expected
	bias = np.sum(bias)
	return(bias)

def scan_planet(real_data=False, nsamp=1000, planet_id=5, nu=100e9, T=110, 
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

	if real_data == False:
		ra, dec = tod.fake_raster(nsamp)
		theta = np.sqrt(ra**2 + dec**2)
		amplitude = 1.0

	else:
	
		ra,dec = jpl.get_planet_timelines(5,20100101,20110101,1)
		ra,dec,theta = get_interpolated_theta(ra,dec,nsamp)
		amplitude = pm.planck_func(nu,T)

	if beam_type=='gaussian':
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

			n_signal = signal+noise

			#fit the curve in the presence of noise
			n_signal = np.real(n_signal)
			mean = sum(ra*n_signal)/len(ra)                   
			sigma = sum(n_signal*(ra-mean)**2)/len(ra) 
			popt,pcov = curve_fit(gaus,ra,n_signal,p0=[1,mean,sigma])
			fitted_values = gaus(ra,*popt)

			#plt.plot(ra,n_signal,label='actual signal')
			#plt.plot(ra,fitted_values,label='fitted signal')

			if compute_bias == True:
				bias = estimate_bias(signal,fitted_values)
				#print(bias)

		print(signal[100],n_signal[100],fitted_values[100],bias)
		return(fitted_values,bias)


def parametrizing_bias(len_sigma0, len_beam_width, len_ra,
						len_dec, nsamp):
	'''
	Work in progress

	Compute bias for all different input features

	Construct a model that can predict bias behavior

	*use deep learning forecasting techniques
	'''

	bias = np.zeros((nsamp,nsamp,nsamp,nsamp))
	sigma0 = np.array(np.linspace(0, len_sigma0,nsamp))
	beam_width = np.array(np.linspace(0,len_beam_width,nsamp))
	ra = np.array(np.linspace(0,len_ra,nsamp))
	dec = np.array(np.linspace(0,len_dec,nsamp))

	for sigma,bw,r,d in sigma0:
		fin, psdin = noise_psd(sigma0, fknee=0.020, alpha=1.5)
		fitted,bias = scan_planet(nsamp,ra=r,dec=d)
		bias_matrix[sigma,bw,r,d] = bias

	np.save('/Users/nadia/Downloads/bias_matrix')
	#plt.plot(sigma0,bias_sigma)


def exclude_fitting_bias():
	'''
	Work in progress

	add more noise each time considering the fitted curve
	as the true one (just signal) and then re-estimate the 
	bias. If bias is increasing then part of it is probably
	also caused by the fitting and not only noise systematics

	'''

scan_planet()





