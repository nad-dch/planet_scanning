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
def psd_noise_model_no_h(f = np.array(np.linspace(0,30e-03,100)),
    white = 10e-05,fknee = 20e-03, alpha =1, fmin = 0.001):
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

def get_tod_from_psd(psd, check=True):
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

    if check == True:

        if scipy.fftpack.fft(onef_tod) == psd:
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


# Method 4: Jon's solution
def extrap(x, xp, yp):
    """
    np.interp function with linear extrapolation
    """
    y = np.interp(x, xp, yp)
    y = np.where(x<xp[0], yp[0]+(x-xp[0])*(yp[0]-yp[1])/(xp[0]-xp[1]), y)
    y = np.where(x>xp[-1], yp[-1]+(x-xp[-1])*(yp[-1]-yp[-2])/(xp[-1]-xp[-2]), y)
    return y

def noise_psd(sigma0=0.01, fknee=0.020, alpha=1.5,
    f=np.logspace(-3, 2, 100)):
    '''
    Generates a noise PSD from a three-paramater noise model
    '''

    psd = (sigma0 ** 2) * (1 + (fknee/f)**alpha)

    return f, psd


def noise_rel(nsamp, fsamp, fin=None, psdin=None):
    '''
    A generic noise generation function

    Args:

        nsamp: is the length of the output timestream mod(nsamp,2) = 0 (a must)
        fsamp: is the sampling frequency [Hz]
        fin: is the frequency values corresponding to the PSD
        psdin: is the input psd

    Returns:
        noise: The noise timeline

    '''

    # Ensure that input parameters are interpreted as floats
    nsamp = int(nsamp)
    fsamp = int(fsamp)

    # Using default values
    if fin is None or psdin is None:
        fin, psdin = noise_psd()

    hnsamp = nsamp/2 # Half the length of the output timestream
    hnsamp=int(hnsamp)
    dt = 1/fsamp
    duration = nsamp*dt

    # The frequency used as the abscissa for interpolation
    f = np.zeros(nsamp)
    f[0]= 0
    f[1:(hnsamp+1)] = np.linspace((fsamp/nsamp), (fsamp/2), (hnsamp))
    ftmp = f[1:(hnsamp)]
    ftmp = ftmp[::-1] # Flipping the vector around
    f[(hnsamp+1):] = ftmp

    # Creating a new psd by interpolating the input psd
    psdi = extrap(f,fin,psdin) # extrap comes from the pdata module
    psdi[0] = 0 # force zero mean

    # Gaussian random vector with zero mean and unit variance
    np.random.seed()
    w = np.random.randn(nsamp)
    ntilde = np.zeros(int(nsamp), dtype='complex64')
    ntilde[0] = w[0]*np.sqrt(psdi[0]/duration) # Zero by construction
    ntilde[hnsamp] = w[-1]*np.sqrt(psdi[hnsamp]/duration)

    # Creating the complex coefficients of the Fourier transform
    ntilde[1:hnsamp] = (w[1:hnsamp]+w[(hnsamp+1):]*1j) \
        *np.sqrt(psdi[1:hnsamp]/duration)/2.0

    # The complex conjugate
    ntmp = np.conjugate(ntilde[1:hnsamp])
    ntmp = ntmp[::-1] # Flipping the vector around
    ntilde[(hnsamp+1):] = ntmp
    # ntilde is now conjugate symmetric
    noise = np.fft.ifft(ntilde)*nsamp

    return noise

def test_noise_model():

    noise_psd = psd_noise_model_no_h()
    f = np.array(np.linspace(0, 30e-2, len(noise_psd)))
    plt.plot(f, noise_psd)
    plt.show()
    noisy_timeline = get_tod_from_psd(noise_psd, check = True)
    plt.plot(f,noisy_timeline)
    plt.show()

def test_option4():

    nsamp = 100000.
    fsamp = 50.
    f, psd = noise_psd()
    noise_out = noise_rel(nsamp, fsamp)

    # Noise timestream
    plt.plot(noise_out)
    #plt.savefig('img/noise_out.png', dpi=300)
    plt.show()

    # Noise PSD
    plt.loglog(f, psd)
    #plt.savefig('img/noise_psd.png', dpi=300)
    plt.show()
def main():

    # test_noise_model()
    test_option4()

if __name__ == '__main__':
    main()
