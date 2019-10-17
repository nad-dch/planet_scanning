import numpy as np
from scipy.fftpack import fft, ifft

def white_noise(nsamp, sigma=1.0, mu=0.0):
    '''
    '''

    return sigma*np.random.randn(nsamp)+mu

def onef_noise(nsamp, sigma=1.0, fknee=20e-04):
    '''
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

def kdnoise_gen(fk, krnl, res, fsamp, seed=None, prng=None, use_numpy=False):
    """
    Gaussian noise generator, for arbitrary input spectrum

    Args:
        fk: (array [arbitrary len]) frequencies of input kernel
        krnl: (array [len of fk]) Noise power spectrum
        res: number of desired output samples
        fsamp: sampling frequency of the tod

    Keywords:
        seed (int): Accepts a seed keyword, which will fix the seed used in the 
            numpy stuff.  This might impact sims, or whatever else would rely
            on seedy things.  As long as a seed value is never presumed, all
            should be ok.
        prng:  If specified, overrides seed. Uses prng is the RandomState object
            from which to get random numbers.
        use_numpy (bool): Uses numpy.fft for Fourier operations if set to True.
            Default is False. 

    Returns:
        A noise TOD, with accompanying depression and sadness.
    """

    dt=float(1)/fsamp

    #If the input resolution is odd, pad it by one for later ease
    if ((res % 2) != 0):
        res += 1
        odd_flag = 1
    else:
        odd_flag = 0

    duration=res*dt

    #Populate the frequency vector
    if use_numpy:
        nps_f=np.fft.fftfreq(res,dt)
    else:
        nps_f=fftpack.fftfreq(res,dt)

    #Interpolate the input spectrum to these frequencies, 
    #and for now enforce a real valued noise power spectrum
    #(nps_f, k) - k complex : useful for cross correlations 
    #(nps_f,nps) - nps real : used here

    #Check that the tabulated frequencies are monotonic
    if (~np.all(np.diff(fk) > 0)):
        print("KDFILTER: input kernel frequencies are not monotonic. Sorting it.")
        inds = np.argsort(fk)
        krnl = krnl[inds]
        fk = fk[inds]

    #Interpolate the input kernel to the sample frequencies
    #Note: might want to make this a polynomial rather than linear
    #Interpolate the real and imaginary bits separately
    krnl_re=np.real(krnl)
    krnl_im=np.imag(krnl)
    k_re=np.interp(np.abs(nps_f),fk,krnl_re)
    k_im=np.interp(np.abs(nps_f),fk,krnl_im)
    #Negative frequencies are conjugate
    idx=np.where(nps_f < 0)
    k_im[idx]*=(-1)
    
    k = k_re + (1j)*k_im
    
    #For now, enforce a real, positive definite kernel
    nps=np.absolute(k)

    #Null the DC level
    nps[0]=0

    #Deal with the seed
    if (seed is not None):
        if prng is None:
            np.random.seed(seed=seed)

    #Generate the random variates
    if prng is not None:
        w = prng.randn(res)
    else:
        w = np.random.randn(res)

    ntilde=np.full(res,0+0j,dtype=complex)
    #Amplitude at DC
    ntilde[0]=w[0]*np.sqrt(nps[0]/duration)
    #Amplitude at Nyquist
    ntilde[res/2]=w[res-1]*np.sqrt(nps[res/2]/duration)
    #Amplitude at positive frequencies have uncorrelated real
    #and imaginary parts

    ctmp_amp = np.sqrt(nps[1:res/2]/duration)/2.0
    ctmp_re = w[1:res/2]*ctmp_amp
    ctmp_im = w[res/2+1:] *ctmp_amp
    ctmp = ctmp_re + (1.j *ctmp_im)

    ntilde[1:res/2]=ctmp
    #Amplitude at negative frequencies are the conjugate of the positive
    ntilde[res/2+1:]=np.conj(ctmp[::-1])

    if use_numpy:
        s=np.fft.fft(ntilde)
    else:
        s=fftpack.fft(ntilde)
    
    if (odd_flag == 1):
        s=s[:res-1]
        
    return np.real(s)


def generate_correlated_noise(noise_cormat, ws=None, Vs=None, seed=None,
    prng=None):
    '''
    Returns N timelines with correlated noise described by noise_cormat.

    noise_cormat (float array): An array of size n_det x n_det x n_samp
        that describes the correlation between detectors.

    Optional Args:
    ws (float array): The eigenvalues of noise_cormat. This is allows the
        algorithm to be more efficient because it does not need to calculate
        them.
    Vs (float array): The eigvenvetors of noise_cormat. This is allows the
        algorithm to be more efficient because it does not need to calculate
        them.
    seed (int): seed the random generator.

    Returns:
    noise (float array): An array of size n_det x n_samp of correlated noise
        timestreams. 
    '''
    ndet, ndet, nfft = np.shape(noise_cormat)

    # Generate random variables 
    if seed is not None:
        if prng is None:
            np.random.seed(seed)

    if prng is not None:
        R = (prng.randn(ndet, nfft)  +
            1.j*prng.randn(ndet, nfft))/np.sqrt(2)
    else:
        R = (np.random.randn(ndet, nfft)  +
            1.j*np.random.randn(ndet, nfft))/np.sqrt(2)

    # Get eigenvals and eigenvecs if not supplied.  
    if (ws is None) or (Vs is None):
        ws, Vs = get_eig(noise_cormat)

    # Set negative eigenvalues to 0 
    if (ws<0).any():
        warnings.warn('Some negative eigenvalues. Setting them to zero.')
        ws[ws<0] = 0

    N = np.zeros((ndet, nfft), dtype=np.complex)
    for i in np.arange(nfft):
        N[:,i] = np.dot(Vs[:,:,i], np.dot(np.diag(np.sqrt(ws[:,i])),
            R[:,i]))

    ntilde = np.concatenate((np.zeros((ndet, 1)), N[:,1:], 
        np.zeros((ndet, 1)), np.conj(np.fliplr(N[:,1:]))), axis=1)
    return np.real(fftpack.ifft(ntilde, axis=1))



def psd_noise_model_no_h(f, white, fknee, alpha, fmin=.02):
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
    return white * (f**alpha + fknee**alpha) / (f**alpha + fmin**alpha)


def test_noise_model():

   	[f,white,fknee,alpha] = [np.random_sample(20e-02),0.00001,20e-04,1]
	noise1 = psd_noise_model_no_h(f, white, fknee, alpha, fmin=.02)
	noise2 = onef_noise(nsamp, sigma=1.0, fknee=20e-04)





def main():

    test_noise_model()

if __name__ == '__main__':
    main()
