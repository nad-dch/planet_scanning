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


def gaus_2D((x, y), mean_x, mean_y, sigma_x, sigma_y, xo, yo, amplitude=1, 
        theta=0, offset=0): 
    '''
    The Gaussian probability function
    Returns
    -------
    Array of Gaussian curve values
    '''
    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
    + c*((y-yo)**2)))

    return (g.ravel())

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
        beam_type='gaussian', noise_type=None, compute_bias=True):

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


    if noise_type == None:
        print('Noise type is None')
        noise = np.zeros_like(signal)

    if noise_type == 'white':
        noise = nm.white(nsamp,sigma=0.01)

    elif noise_type =='onef':
        noise = nm.noise_rel(nsamp, fsamp=5)

    else:
        raise ValueError('Noise is either None or in [white, onef]')

    n_signal = signal + noise

    #fit the curve in the presence of noise
    n_signal = np.real(n_signal)
    mean_ra = sum(ra*n_signal)/len(ra)
    sigma_ra = sum(n_signal*(ra-mean_ra)**2)/len(ra)
    mean_dec = sum(dec*n_signal)/len(dec)
    sigma_dec = sum(n_signal*(dec-mean_dec)**2)/len(dec)
    initial_guess = (1,mean_ra,mean_dec,sigma_ra,sigma_dec,0,0)
    print(len(ra),len(dec),len(n_signal))
    popt,pcov = curve_fit(gaus_2D,(ra,dec),n_signal,p0=initial_guess)
    fitted_values = gaus_2D((ra,dec),*popt)

    print(np.shape(ra))
    plt.plot(ra, n_signal, ls='', marker='.', label='actual signal')
    plt.plot(ra, fitted_values, ls='', marker='.', label='fitted signal')
    plt.legend()
    plt.show()


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

    for sigma,bw,r,d in zip(sigma0,beam_width,ra,dec):
        fin, psdin = noise_psd(sigma0, fknee=0.020, alpha=1.5)
        fitted,bias = scan_planet(nsamp,ra=r,dec=d)
        bias_matrix[sigma,bw,r,d] = bias

    np.save('/Users/nadia/Downloads/bias_matrix')
    #plt.plot(sigma0,bias_sigma)


def exclude_fitting_bias(max_iterations):
    '''
    Work in progress

    add more noise each time considering the fitted curve
    as the true one (just signal) and then re-estimate the 
    bias. If bias is increasing then part of it is probably
    also caused by the fitting and not only noise systematics

    '''
    nsamp,fsamp = 1000,50
    noise = nm.noise_rel(nsamp,fsamp)
    fitted_values,bias = scan_planet()
    rel_bias = []
    for i in range (max_iterations):
        bias0 = bias
        true_signal = fitted_values 
        noisy_signal = true_signal+noise 
        fitted_values = fit_noise()
        bias = estimate_bias(true_signal,fitted_values)
        rel_bias.append(bias-bias0)
        
    plt.plot(rel_bias)


def main():

    fitted_values, bias = scan_planet(noise_type='onef')

if __name__ == '__main__':
    main()
#print(fitted_values)
#print(bias)