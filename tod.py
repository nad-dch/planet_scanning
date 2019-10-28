import numpy as np

def fake_raster(cr=[-60, -60, 60, 60], numel=[1001, 1001, 0]):
    '''
    Creates a fake raster scan pointing timeline to use for generating timelines
    '''

    N = np.sqrt(numel[0]*numel[1]) + 1
    RA, DEC = np.meshgrid(np.linspace(cr[0], cr[2], numel[0]),
        np.linspace(cr[1], cr[3], numel[1]))
    ra, dec = RA.flatten(), DEC.flatten()

    x = np.arange(len(ra))
    xi = np.linspace(0, len(x), numel[0]*numel[1])
    rai = np.interp(xi, x, ra)
    deci = np.interp(xi, x, dec)

    return rai, deci


