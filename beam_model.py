import time
import numpy as np
from scipy import optimize
import beam_tools as bt

def gaussian_beam(theta, fwhm=40./60.):
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

def eg(p, x, y, fit4dc=False):

    # X offset
    beamx = p[0]
    # Y offset
    beamy = p[1]
    # rotation angle
    gam = p[2]
    # mean FWHM
    fwhm_mean = p[3]
    #ellipticity
    ellip = p[4]
    # amplitude
    A = p[5]
    # baseline offset from zero
    if fit4dc:
        offset = p[6]
    else:
        offset = 0.0

    p[2] = np.fmod(p[2], 180.0)
    sgam = np.sin(gam/180.0*np.pi)
    cgam = np.cos(gam/180.0*np.pi)

    fwhm_x = fwhm_mean/np.sqrt(ellip)
    sig_x = fwhm_x /np.sqrt(8*np.log(2))
    sig_y = sig_x*ellip

    u = x-beamx
    v = y-beamy

    xp = cgam*u+sgam*v
    yp = -sgam*u+cgam*v

    return A*np.exp(-(xp/sig_x)**2 /2 - (yp/sig_y)**2 / 2) + offset

def eg_res(p, x, y, data, sigma, fit4dc=False):

    out = ((data-eg(p, x, y, fit4dc=fit4dc))/sigma).flatten()

    return out

def egfit(Data, fit_radius=0.5, fit4dc=False, pguess=None, benchmark=False,
    force_center=False):

    x, y, data = Data['x'], Data['y'], Data['data']
    if fit_radius is not None:

        if (pguess is None) or force_center:
            r = np.sqrt(x**2 + y**2)
        else:
            r = np.sqrt((x-pguess[0])**2 + (y-pguess[1])**2)

        ifit = (r < fit_radius)
        if len(ifit) == 0:
            return 0.0, None, 0.0, 0.0

        sigma = np.std(data[np.where(r > fit_radius)])
        if fit4dc and pguess is None:
            pguess = [0.0, 0.0, 1.0, 0.5, 1.0, data[ifit].max(), 0.0]
        elif pguess is None:
            pguess = [0.0, 0.0, 1.0, 0.5, 1.0, data[ifit].max()]

        po, cov, infodict, errmsg, success = optimize.leastsq(eg_res, pguess,
            args=(x[ifit], y[ifit], data[ifit], sigma, fit4dc),
            full_output=1)

        model = eg(po, x, y)
        res = np.sum(np.abs(data[ifit]-model[ifit]))/len(model[ifit]) \
            /np.mean(data[ifit])

    else:
        sigma = np.std(data)
        if fit4dc and pguess is None:
            pguess = [0.0, -5.0, 1.0, 0.5, 1.0, data.max(), 0.0]
        elif pguess is None:
            pguess = [0.0, -5.0, 1.0, 0.5, 1.0, data.max()]

        po, cov, infodict, errmsg, success = optimize.leastsq(eg_res, pguess,
            args=(x, y, data, sigma, fit4dc), full_output=0)

        model = eg(po,x,y)
        res = np.sum(np.abs(data-model))/float(len(model))/np.mean(data)

    return po, cov, model, res

def gfit(cr, numel, arr1, verbose=True, gfwhm=0.4, gell=0.01, fit_radius=1.0,
    return_model=False, benchmark=False):

    t0 = time.time()

    xx, yy, dA = bt.get_mesh(cr, numel)
    x, y = bt.get_lin(cr, numel, mult=1)

    t2 = time.time()

    t3 = time.time()

    Data = {'data':arr1, 'x':xx, 'y':yy}

    srln = np.sqrt(8*np.log(2))
    pguess1 = [0., 0., 0., gfwhm, gell, np.max(arr1)]
    po, cov, model_out1, res1 = egfit(Data, fit_radius=fit_radius,
        fit4dc=False, pguess=pguess1)

    t4 = time.time()

    fwhm_x = po[3]/np.sqrt(po[4])
    sx = fwhm_x /np.sqrt(8*np.log(2))
    sy = sx*po[4]

    sx = np.abs(sx)
    sy = np.abs(sy)

    cx = po[0]
    cy = po[1]

    ellipticity = (np.max([sx, sy])-np.min([sx, sy])) / \
        (sx+sy)

    t5 = time.time()

    if verbose:

        print('Elliptical Gaussian beam fit results:')
        print('  Amplitude      : {:5.2f}'.format(po[5]))
        print('  X centroid     : {:5.3f} deg'.format(po[0]))
        print('  Y centroid     : {:5.3f} deg'.format(po[1]))
        print('  FWHM              : {:5.3f} arcmin'.format(np.abs(po[3])*60))
        print('  Rotation angle : {:5.2f} deg'.format(po[2]))
        print('  Ellipticity    : {:5.4f}'.format(ellipticity))

    if benchmark:
        print('Timing results:')
        for t in [t0, t1, t2, t3, t4, t5]:
            print('  t = {}'.format(t-t0))

    if return_model:
        return cx, cy, sx, sy, np.abs(po[3]), ellipticity, cr, numel, model_out1

    return cx, cy, sx, sy, np.abs(po[3]), ellipticity

def main():

    return

if __name__ == '__main__':

    main()
