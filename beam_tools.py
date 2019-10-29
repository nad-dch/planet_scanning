import sys, os
import matplotlib
matplotlib.use('Agg')
import numpy as np
import scipy.interpolate
from scipy import optimize
import time

def center_idxs(arr):
    '''
    Returns the indexes corresponding to the location of the peak pixel value

    The first dimension returned is y-axis (vertical)
    The second dimension returned is x-axis (horizontal)
    '''

    return np.unravel_index(np.argmax(arr), np.shape(arr))

def center_map(cr, numel, arr, ix=None):
    '''
    Centers a map on the peak value
    '''

    X, Y, _ = get_mesh(cr, numel)

    if ix is None:

        idxs, dists = dist2edge(cr, numel, arr)
        sh = np.min(dists)
        xl, xh = (idxs[1]-sh), (idxs[1]+sh-1)
        yl, yh = (idxs[0]-sh), (idxs[0]+sh-1)

    else:

        xl, xh = ix[2], ix[3]
        yl, yh = ix[0], ix[1]

    cr2use = [X[0,xl], Y[yl,0], X[0,xh], Y[yh,0]]
    numel2use = [yh-yl, xh-xl, 0]
    arr2use = arr[yl:yh, xl:xh]

    return [yl, yh, xl, xh], cr2use, numel2use, arr2use

def get_lin(cr, numel, mult=10):

    '''
    get the linear space for given cr and numel
    '''

    x2 = np.linspace(cr[0], cr[2], mult*numel[0])
    y2 = np.linspace(cr[1], cr[3], mult*numel[1])

    return x2, y2

def get_da(cr, numel, mult=1):

    # return ((cr[2]-cr[0])/(mult*numel[0]))**2
    return ((float(cr[2])-float(cr[0]))/(mult*float(numel[0])))**2

def get_mesh(cr, numel, mult=1.0):

    '''
    Return shaped 2-dimensional arrays
    '''

    xx, yy = np.meshgrid(np.linspace(cr[0], cr[2], mult*numel[0]),
        np.linspace(cr[1], cr[3], mult*numel[1]))

    dA = get_da(cr, numel, mult=1)

    return xx, yy, dA

def dist2edge(cr, numel, arr):

    idxs = center_idxs(arr)
    #dists = np.array([idxs[0], numel[1]-idxs[0], idxs[1], numel[0]-idxs[1]])
    dists = np.array([idxs[0], numel[0]-idxs[0], idxs[1], numel[1]-idxs[1]])

    return idxs, dists

def centerofmass(arr):
    '''
    Estimates the center of mass and the width of the Gaussian distribution
    '''

    nx, ny = np.shape(arr)
    mx, my = np.sum(arr, axis=0), np.sum(arr, axis=1)
    mx = mx*(mx > 0).astype(float)
    my = my*(my > 0).astype(float)

    x, y = np.arange(nx), np.arange(ny)
    cx, cy = np.sum(mx*x)/np.sum(mx), np.sum(my*y)/np.sum(my)
    sx = np.sqrt(np.sum(mx*np.abs(x-cx)**2)/np.sum(mx))
    sy = np.sqrt(np.sum(my*np.abs(y-cy)**2)/np.sum(my))

    return cx, cy, sx, sy
