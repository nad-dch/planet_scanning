import numpy as np
import scipy
from scipy.fftpack import fft, ifft
import matplotlib.pyplot as plt

def white_noise(nsamp, sigma=1.0, mu=0.0):
    '''
    '''

    return sigma*np.random.randn(nsamp)+mu


# Jon's solution for 1/f noise

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
    cut_last=False
    nsamp = float(nsamp)
    fsamp = float(fsamp)

    # Using default values
    if fin is None or psdin is None:
        fin, psdin = noise_psd()

    if nsamp%2 == 1:
        nsamp += 1
        cut_last=True

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
    if cut_last:
        noise = noise[:-1]

    return np.real(noise)


def test_noise_model():

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
    test_noise_model()

if __name__ == '__main__':
    main()
