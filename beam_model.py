import numpy as np

def gaussian_beam(theta, fwhm=40.0):
    '''

    Arguments:
    theta:

    kwargs:
    fwhm : beam FWHM in arcmin
    '''

    sigma = fwhm / np.sqrt(8*np.log(2))

    return np.exp(-0.5*theta**2/sigma**2)

def H_beam():

    '''
    construct Hiroaki's beam(band-averaged 
    truncated Gaussian)
    '''
def po_beam():

    '''
    '''
def main():

    return

if __name__ == '__main__':

    main()
