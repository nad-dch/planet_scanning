
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import numpy as np
import tod
import noise as nm
import beam_model as bm
import planet_model as pm

opj = os.path.join

# read planet timelines
# generate beam time using beam_model.py and planet_model.py
# add noise
# plot timeline

def test_fake_scan(nsamp=1000, nu=100e9, T=100):
    '''
    The purpose of this function is to test simple
    '''

    ra, dec = tod.fake_raster(nsamp)

    theta = np.sqrt(ra**2 + dec**2)

    # amplitude = pm.planck_func(nu, T)
    amplitude = 1.0
    signal = amplitude * bm.gaussian_beam(theta)
    noise = nm.white_noise(nsamp, sigma=0.01)
    signal += noise

    print(amplitude)

    plt.plot(ra, signal, lw=0, marker='.', ms=8)
    plt.ylabel('Signal amplitude')
    plt.xlabel('RA [arcmin]')
    plt.savefig(opj('img/', 'tod.png'), dpi=300, bbox_inches='tight')
    plt.close()

    plt.semilogy(ra, signal, lw=0, marker='.', ms=8)
    plt.ylabel('Signal amplitude')
    plt.xlabel('RA [arcmin]')
    plt.savefig(opj('img/', 'tod_log.png'), dpi=300, bbox_inches='tight')
    plt.close()

def main():

    test_fake_scan(nsamp=1000)

if __name__ == '__main__':
    main()