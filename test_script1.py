
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import numpy as np
import tod
import noise as nm
import beam_model as bm

opj = os.path.join

# read planet timelines
# generate beam time using beam_model.py and planet_model.py
# add noise
# plot timeline

def test_fake_scan():
    '''
    The purpose of this function is to test simple
    '''

    nsamp = 1000
    ra, dec = tod.fake_raster(nsamp)

    theta = np.sqrt(ra**2 + dec**2)
    signal = bm.gaussian_beam(theta)
    noise = nm.white_noise(nsamp, sigma=0.1)
    signal += noise

    #print(np.sum)

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

    test_fake_scan()

if __name__ == '__main__':
    main()