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
import pickle
opj = os.path.join

def get_interpolated_theta(ra, dec, numel):
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

    return rai,deci,theta_int

def estimate_bias(expected, measured):
    '''
    Returns the bias of the estimation
    '''

    bias = ((measured-expected)**2)/expected
    bias = np.sum(bias)
    return(bias)

def create_data(dates=None, real_data=False,cr=[-1, -1, 1, 1],
    numel=[101, 101, 0]):

    '''
    Create a real or fake dataset
    '''

    if not real_data:
        nsamp = numel[0]*numel[1]
        ra, dec = tod.fake_raster(cr=cr, numel=numel)
        theta = np.sqrt(ra**2 + dec**2)
        amplitude = 1.0

    else:

        if dates is None:
            ra,dec = jpl.get_planet_timelines(5,20100101,20110101,1)
        else:    
            ra,dec = jpl.get_planet_timelines(5,dates[0],dates[1],1)

        ra,dec,theta = get_interpolated_theta(ra, dec, numel=numel)
        amplitude = pm.planck_func(nu,T)

    return ra,dec,theta,amplitude

def fit_signal(ra, dec, signal, noise, n_signal, cr=[-1, -1, 1, 1],
                numel=[101, 101, 0], make_plots=False):
    '''
    '''
    #fit the curve in the presence of noise
    mean_ra = np.sum(ra*n_signal)/len(ra)
    sigma_ra = np.sum(n_signal*(ra-mean_ra)**2)/len(ra)
    mean_dec = np.sum(dec*n_signal)/len(dec)
    sigma_dec = np.sum(n_signal*(dec-mean_dec)**2)/len(dec)
    initial_guess = (1, mean_ra, mean_dec, sigma_ra, sigma_dec, 0, 0)

    n_signal = np.reshape(n_signal, (numel[0], numel[1]))
    noise = np.reshape(noise, (numel[0], numel[1]))

    cx, cy, sx, sy, angle, e = \
    bm.gfit(cr, numel, n_signal, gfwhm=40./60,
    gell=0.01, fit_radius=1.0, return_model=False,
    verbose=False)

    if make_plots:
        # Plotting raw (noisy) beam map
        pt.plot_beam(cr, numel, n_signal, fname='noisy_signal', log=False)

        pt.plot_beam(cr, numel, model_out, fname='fit', log=False)

        pt.plot_beam(cr, numel, n_signal-model_out, vmin=-0.02, vmax=0.02,
        fname='diff', log=False)

        pt.plot_beam(cr, numel, noise, vmin=-0.02, vmax=0.02,
        fname='noise', log=False)

        # Plotting timelines as a function of az
        plt.plot(ra, n_signal.flatten(), ls='', marker='.', label='signal+noise')
        plt.plot(ra, model_out.flatten(), ls='', marker='.',
        label='fitted signal', alpha=0.5)
        plt.plot(ra, signal.flatten(), ls='', marker='.', label='signal')
        plt.legend()
        plt.savefig(opj('img/', 'comparison.png'), dpi=300)

    return cx, cy, sx, sy, np.sqrt(sx*sy)*np.sqrt(8*np.log(2)), angle, e


def scan_planet(Data=None, real_data=False, numel=[101, 101, 0], planet_id=5, nu=100e9, T=110,
                beam_type='gaussian', p=None, fwhm=40./60., angle=30., ec=1.5, 
                noise_type='random', fit=True, compute_bias=False):

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
    if Data == None:
        ra,dec,theta,amplitude = create_data()

    else:
        ra,dec,theta,amplitude = Data[0], Data[1], Data[2], Data[3]

    if beam_type=='gaussian':
        signal = amplitude * bm.gaussian_beam(theta, fwhm=fwhm)

    elif beam_type == 'elliptical':

        if p == None:
            p = [0, 0, angle, fwhm, ec, amplitude]

        signal = bm.eg(p, ra, dec)

    if noise_type == None:
        print('Noise type is None')
        noise = np.zeros_like(signal)

    if noise_type == 'white':
        noise = nm.white_noise(numel[0]*numel[1], sigma=0.01)

    if noise_type =='onef':
        noise = nm.noise_rel(numel[0]*numel[1], fsamp=5)

    elif noise_type == 'random':
        noise = np.random.randn(len(signal))


    else:
        raise ValueError('Noise is either None or in [white, onef, random]')

    n_signal = signal + noise

    if fit:

        cx, cy, sx, sy, fwhm, angle, e = fit_signal(ra, dec, signal, noise, n_signal)
        return cx, cy, sx, sy, fwhm, angle, e

        # Needs work
        if compute_bias:

            bias = estimate_bias(signal, fitted_values)
            return fitted_values,bias

    #the fitted model compared to the noisy model to check the fitting
    #the input signal compared to the output after the fitting
    


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

    return rel_bias

def make_dirs():

    if not os.path.exists('img/'):
        os.makedirs('img/')

def run_sims(sim_number, pace=1, parameter='noise', plot_comparison=True):

    comparison = np.zeros((sim_number,7))

    make_dirs()

    if parameter == 'noise':

        for i in range(sim_number):
            print('  sim {}/{}'.format(i, sim_number))
            comparison[i,:] = scan_planet(noise_type='random')

    if parameter == 'fwhm':

        fwhm = 38./60

        for i in range(sim_number):
            comparison[i,:] = scan_planet(fwhm=fwhm)
            comparison[i,5] -= fwhm
            fwhm = fwhm + pace


    if parameter == 'nsamp':

        nsamp = 1000

        for i in range(sim_number):
            comparison[i,:] = scan_planet(numel=[np.sqrt(nsamp),
                                                np.sqrt(nsamp),0])
            nsamp = nsamp + pace

    elif parameter == 'ra_dec':

        start_date = 20100101
        end_date = 20110101

        for i in range(sim_number):

            # needs further work

            ra,dec = jpl.get_planet_timelines(5,start_date,end_date,1)
            ra,dec,theta = get_interpolated_theta(ra, dec, numel=numel)
            amplitude = pm.planck_func(nu,T)
            Data = [ra,dec,theta,amplitude]
            comparison[i,:] = scan_planet(Data=Data)
            start_date = start_date + pace
            end_date = end_date + pace


    with open('run_sims.pkl','wb') as f:
        pickle.dump(comparison, f)

    if plot_comparison:

        with open('run_sims.pkl','rb') as f:
            comparison = pickle.load(f)
        [n,bins] = np.histogram(comparison[:,5], bins=31)
        plt.xlabel('fwhm difference')
        plt.ylabel('counts')
        plt.plot(bins[:-1], n, label='fwhm of the fits')

        plt.savefig(opj('img/', 'comparison_fwhm.png'))

    return comparison



def main():

    run_sims(100)

if __name__ == '__main__':
    main()


