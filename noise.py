import numpy as np
import scipy
from scipy.fftpack import fft, ifft
import matplotlib.pyplot as plt

def white_noise(nsamp, sigma=1.0, mu=0.0):
    '''
    '''

    return sigma*np.random.randn(nsamp)+mu

#method 1: analytical computations using scipy.fft

def onef_noise(nsamp, sigma=1.0, fknee=20e-03):
    '''
    Generate the desired 1/f noise tod

    Arguments:
	    nsamp : length of the sample
	    fknee : knee frequency, e.g., value at which the frequency equals
	    		2 times the white noise floor
		fs : the frequency sample 1/nsamp
    '''
    C_k=[]
    fs = 2*fknee*np.random.sample(nsamp)
    fk = scipy.fftpack.fft(np.array(fs))
    phase = 2*np.pi*np.random.sample(len(fk))
    C_k = 1/np.abs(fk)
    C_k = C_k*np.exp(1j*phase)
    noise = scipy.fftpack.ifft(C_k)
    return noise
    #pass



# method 2: using psd_noise_model_no_h and converting psd -->tod



def psd_noise_model_no_h(f = np.array(np.linspace(0,30e-03,100)), white = 10e-05,fknee = 20e-03, alpha =1, fmin = 0.001):
    '''
    The same as psd_noise_model, but with know hknee.

    Arguments:
        f : Array of frequencies
        white : White noise floor
        fknee : Knee frequency of the "1/f noise"
        alpha : Power-law slope of the "1/f noise"

    Optional Args:
        fmin : Minimum frequency.
            This can either describe a spectrum that levels out at low
            frequency, or be set to a small value to prevent diverging
            as f approaches 0.    

    Returns:
        psd (float): The psd model as a function of f.
    '''
    
    white = white**2
    #alpha = 2*alpha
    return white * (fknee**alpha) / (f**alpha + fmin**alpha)

def get_tod_from_psd(psd, check = True):
    '''
    The purpose of this function is to return the timeline data given a psd

    Arguments:
    	psd   : Array of psd values
    	check : bool,optional. Performs the inverse transformation to check
    			the output. 
	Returns:
		timeline data : array of floats
	'''
	amplitude = np.sqrt(2*psd)
	phase = 2*np.pi*np.random.sample()
	z = amplitude*np.exp(j*phase)
	onef_tod = scipy.fftpack.ifft(z)

	if check = True:

		if scipy.fftpack.fft(onef_tod) = psd:
			print('the function is working')
		else:
			print('change your function')

	return(onef_tod)

# method 3: calculate a 1/f Fourier filter kernel amd interpolate it to the frequencies of the TOD

def notch_kernel(f, fc, w, ns=5):
    """
    A very stupid (gaussian) notch filter kernel

    Args:
        f: array of frequencies at which to sample [Hz]
        fc: the notch center [Hz]
        w: the width of the notch [Hz]

    Returns:
        Array of kernel values at each frequency
    """

    f_low = fc - w*ns
    f_hi = fc + w*ns
    if (f_low<0):
        f_low=0.0

    krnl = np.ones_like(f)
    cond = (f > f_low) & (f < f_hi)
    krnl[cond] = 1.0 - np.exp(-0.5*((f[cond]-fc)/w)**2)
    
    if not np.all(np.isfinite(krnl)):
        raise ValueError("Invalid value in notch_kernel")

    return krnl

def interp_kernel(freq , krnl, nsamp=None, fsamp=None, new_freq=None):
    """
    Interpolate a Fourier filter kernel.
    It takes the kernel provided, which can be of arbitrary size
    and interpolates it to the frequencies of the TOD samples.

    Args:
        freq: (array [arbitrary len]) frequencies of input kernel
        krnl: (array [len of freq]) kernel, can be complex if you like
        nsamp: length of frequency array for new samples
        fsamp: sampling frequency of the tod.
               Used with nsamp to generate new frequencies with fftfreq
        new_freq: (array) array of new frequencies. Overrides nsamp, fsamp.

    Returns:
        Interpolated kernel, joy and happiness.

    Note:
        Currently assumes that freq > 0 and that krnl is conjugate-symmetric.
    """
    if new_freq is not None:
        f = new_freq
    elif nsamp is not None and fsamp is not None:
        dt = 1./fsamp
        #Generate the frequencies of the input TOD
        #for even input, packing is 0 to nyquist, then -(nyquist-1) to 1
        #for odd input, packing is 0 to nyquist, -nyquist to 1
        f = np.fft.fftfreq(nsamp, dt)
    else:
        raise ValueError("Need either nsamp and fsamp, or new_freq")
    
    #Check that the tabulated frequencies are monotonic
    if not np.all(np.diff(freq)):
        raise ValueError("input kernel frequencies are not monotonic.")

    #Interpolate the input kernel to the sample frequencies
    #Note: might want to make this a polynomial rather than linear
    #Interpolate the real and imaginary bits separately
    k_re = np.interp(np.abs(f), freq, krnl.real)
    k_im = np.interp(np.abs(f), freq, krnl.imag)
    k_im[f < 0] *= (-1)         #Negative frequencies are conjugate 
    return k_re + (1j)*k_im

# test

def test_noise_model():
	f = np.array(np.linspace(0,30e-02,10))
	noise_psd = psd_noise_model_no_h()
	plt.plot(f,noise_psd)
	plt.show()
	noisy_timeline = get_tod_from_psd(noise_psd, check = True)
        plt.plot(f,noisy_timeline)
        plt.show()

def main():

    test_noise_model()

if __name__ == '__main__':
    main()
