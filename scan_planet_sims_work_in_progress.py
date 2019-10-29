import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import os
import numpy as np
import tod
import noise as nm
import beam_model as bm
import planet_model as pm
import plotting_tools as pt
import example_jpl as jpl
opj = os.path.join

def get_interpolated_theta(ra,dec,numel=[101, 101, 0]):
    '''
    Interpolates ra and dec to correspond to
    the given sample length (numel)
    Returns
    -------
    theta_int : float
        The angles for the given sample length
    '''

    x = np.arange(len(ra))
    xi = np.linspace(0, len(x), numel[0]*numel[1])
    rai = np.interp(xi, x, ra)
    deci = np.interp(xi, x, dec)
    theta_int = np.sqrt(rai**2 + deci**2)

    return(rai,deci,theta_int)

def estimate_bias(expected, measured):
    '''
    Returns the bias of the estimation
    '''

    bias = ((measured-expected)**2)/expected
    bias = np.sum(bias)
    return(bias)

def scan_planet(real_data=False, nsamp=1000, planet_id=5, nu=100e9, T=110,
        cr=[-1, -1, 1, 1], numel=[101, 101, 0], beam_type='elliptical', p=None,
        fwhm=40.0, noise_type=None, compute_bias=True):

    '''
    A function that outputs the signal tod from a planet scanning
    Arguments:
        real_data : bool. If True use planet timelines otherwise 
            create fake ra and dec through tod.py functions.
        nsamp : sample length
        planet_id : define the planet you want to scan
        nu : frequency
        T : blackbody temperature of the planet (110K for Jupyter)
        cr : individual limits for the parameters in degrees
        numel : number of values inside the linear spaces
        beam_type : gaussian, physical optics, elliptical gaussian
        p : initial guess
        fwhm : full-width-half-maximum that characterizes the beam
        noise_presence : bool.If True choose noise_type.
        noise_type : string 'gaussian','onef'
        compute_bias : bool. default value = True
    Returns:
    -------
        signal or fitted_values: Array of amplitude signal or
                                    noisy signals fitted values
        bias : the bias of the estimation
    '''

    if not real_data:
        nsamp = numel[0]*numel[1]
        ra, dec = tod.fake_raster(cr=cr, numel=numel)
        theta = np.sqrt(ra**2 + dec**2)
        amplitude = 1.0

    else:
        ra,dec = jpl.get_planet_timelines(5,20100101,20110101,1)
        ra,dec,theta = get_interpolated_theta(ra, dec, numel=numel)
        amplitude = pm.planck_func(nu,T)

    if beam_type=='gaussian':
        signal = amplitude * bm.gaussian_beam(theta, fwhm=fwhm)

    elif beam_type == 'elliptical':

        if p == None:
            p = [0,0,0,fwhm,0.01,0]

        signal = amplitude * bm.eg(p,ra,dec)

    if noise_type == None:
        print('Noise type is None')
        noise = np.zeros_like(signal)

    if noise_type == 'white':
        noise = nm.white_noise(nsamp, sigma=0.01)

    elif noise_type =='onef':
        noise = nm.noise_rel(nsamp, fsamp=5)

    else:
        raise ValueError('Noise is either None or in [white, onef]')

    n_signal = signal + noise

    #fit the curve in the presence of noise
    mean_ra = np.sum(ra*n_signal)/len(ra)
    sigma_ra = np.sum(n_signal*(ra-mean_ra)**2)/len(ra)
    mean_dec = np.sum(dec*n_signal)/len(dec)
    sigma_dec = np.sum(n_signal*(dec-mean_dec)**2)/len(dec)
    initial_guess = (1, mean_ra, mean_dec, sigma_ra, sigma_dec, 0, 0)

    n_signal = np.reshape(n_signal, (numel[0], numel[1]))
    noise = np.reshape(noise, (numel[0], numel[1]))

    cx, cy, sx, sy, angle, e, cr_eg, numel_eg, model_out = \
        bm.gfit(cr, numel, n_signal, gfwhm=40./60,
            gell=0.01, fit_radius=1.0, return_model=True,
            verbose=True)

    # Plotting raw (noisy) beam map
    pt.plot_beam(cr, numel, n_signal, fname='noisy_signal', log=False)

    pt.plot_beam(cr, numel, model_out, fname='fit', log=False)

    pt.plot_beam(cr, numel, n_signal-model_out, vmin=-0.02, vmax=0.02,
        fname='diff', log=False)

    pt.plot_beam(cr, numel, noise, vmin=-0.02, vmax=0.02,
        fname='noise', log=False)

    # Plotting timelines as a function of az
    plt.plot(ra, n_signal.flatten(), ls='', marker='.', label='actual signal')
    plt.plot(ra, model_out.flatten(), ls='', marker='.',
        label='fitted signal', alpha=0.5)
    plt.legend()
    plt.savefig(opj('img/', 'comparison.png'), dpi=300)

    # Needs work
    if compute_bias == True:
        bias = estimate_bias(signal, fitted_values)
        print(signal[100],n_signal[100],fitted_values[100],bias)
        return(fitted_values,bias)

    return


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

    scan_planet(noise_type='white', fwhm=35.0/60.,
        compute_bias=False)

if __name__ == '__main__':
    main()
#print(fitted_values)
#print(bias)