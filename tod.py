import numpy as np

def fake_raster(nsamp, box=60):
    '''
    Creates a fake raster scan pointing timeline to use for generating timelines
    '''

    N = np.sqrt(nsamp) + 1
    RA, DEC = np.meshgrid(np.linspace(-box, box, N), np.linspace(-box, box, N))
    ra, dec = RA.flatten(), DEC.flatten()

    x = np.arange(len(ra))
    xi = np.linspace(0, len(x), nsamp)
    rai = np.interp(xi, x, ra)
    deci = np.interp(xi, x, dec)

    return rai, deci