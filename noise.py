import numpy as np
import scipy
from scipy.fftpack import fft, ifft
import matplotlib.pyplot as plt

def white_noise(nsamp, sigma=1.0, mu=0.0):
    '''
    '''

    return sigma*np.random.randn(nsamp)+mu

#method 1: analytical computations using scipy.fft

def onef_noise(nsamp, sigma=1.0, fknee=20e-04):
    '''
    Generate the desired 1/f noise tod

    Arguments:
	    nsamp : length of the sample
	    fknee : knee frequency, e.g., value at which the frequency equals
	    		2 times the white noise floor
		fs : the frequency sample 1/nsamp
    '''
    fs = 20e-04*np.random_sample(nsamp)
    fk = scipy.fftpack.fft(fs, n=None, axis=-1, overwrite_x=False)
    phase = 2*np.pi*np.random_sample(nsamp)
    C_k[0] = 0
    C_k = 1/np.abs(fk)
    C_k = C_k*np.exp(i*phase)
    noise = scipy.fftpack.ifft(C_k, n=None, axis=-1, overwrite_x=False)

    return noise
    #pass



# method 2: using psd_noise_model_no_h and converting psd -->tod



def psd_noise_model_no_h(f, white = 10e-05,fknee = 20e-03, alpha = 1, fmin = 0.1):
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
    alpha = 2*alpha

    return white * ((f**alpha + fknee**alpha) / (f**alpha + fmin**alpha))


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
	phase = np.random_sample()
	z = amplitude*np.exp(j*phase)
	onef_tod = scipy.fftpack.ifft(z)

	if check = True:

		if scipy.fftpack.fft(onef_tod) = psd:
			print('the function is working')
		else:
			print('change your function')

	return(onef_tod)



def test_noise_model():

	f = np.array(np.linspace(0,0.040,10))
	oise1 = psd_noise_model_no_h()
	plt.plot(f,oise1)
	plt.show()





def main():

    test_noise_model()

if __name__ == '__main__':
    main()
