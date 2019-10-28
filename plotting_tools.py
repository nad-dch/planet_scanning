import os
import matplotlib
import itertools
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import patches
from matplotlib.collections import PatchCollection
import matplotlib.cm as cm
import scipy
import numpy as np
import beam_tools as bt


cmap = plt.get_cmap('viridis')
cmap_log = plt.get_cmap('magma')
cmap_grey = plt.get_cmap('Greys')

def configure_labels(ax=None, title=None, xlab=None, ylab=None,
    plot_displacement=False):

    if ax is None:

        if title is None:
            plt.title('Field amplitude')
        else:
            plt.title(title)

        if plot_displacement:
            plt.xlabel('Displacement [m]')
            plt.ylabel('Displacement [m]')

        if xlab is None or ylab is None:
            plt.xlabel('Az [deg]')
            plt.ylabel('El [deg]')
        else:
            plt.xlabel(xlab)
            plt.ylabel(ylab)

    else:

        if title is None:
            ax.set_title('Field amplitude')
        else:
            ax.set_title(title)

        if plot_displacement:
            ax.set_xlabel('Displacement [m]')
            ax.set_ylabel('Displacement [m]')

        if xlab is None or ylab is None:
            ax.set_xlabel('Az [deg]')
            ax.set_ylabel('El [deg]')
        else:
            ax.set_xlabel(xlab)
            ax.set_ylabel(ylab)

def plot_beam(cr, numel, arr, title=None, fname='beam', figsize=None,
    dpi=150, xlab=None, ylab=None, vmin=None, vmax=None,
    vminl=None, vmaxl=None, interp=True, mult=5, logrange=None, aspect='equal',
    scatter_plot=False, xs=None, ys=None, scatter_arr=None,
    truncate=False, truncr=0.2, clab=None, save=True, plot_pix=False,
    cmapl=plt.get_cmap('magma'), cmap=plt.get_cmap('viridis'),
    log=True, imgd='img/', nticks=5, xlim=None, ylim=None, relative_levels=True,
    add_contour=False, contlabel=True, levels=None, lw=0.5, add_grids=False,
    plot_circ=False, circle_color='white', circle_ls='-', plot_ellip=False, crad=1,
    ellipse=[0, 0, 1, 1, 0], ellipse_pars=None,
    plot_squares=False, square_width=87.5,
    xsq=[-2, -1, 0, 1, -2, -1, 0, 1],
    ysq = [-1, -1, -1, -1, 0, 0, 0, 0],
    anno_str=None):
    '''
    Plots a beam
    '''

    xsq = square_width * np.array(xsq)
    ysq = square_width * np.array(ysq)

    extent = [cr[0], cr[2], cr[1], cr[3]]
    xx, yy, dA = bt.get_mesh(cr, numel)

    xx0, yy0 = xx.copy(), yy.copy()

    if not os.path.exists(imgd):
        os.makedirs(imgd)

    fig = plt.figure(figsize=figsize)
    ax = plt.subplot(111)
    if not interp:
        mult = 1
        arr2use = arr.copy()
    else:
        xx, yy, dA = bt.get_mesh(cr,numel)
        xx2, yy2 = bt.get_lin(cr,numel, mult=mult)

        fint = scipy.interpolate.RectBivariateSpline(xx[0,:], yy[:,0], arr)
        arr2use = fint(xx2, yy2)
        xx, yy, dA = bt.get_mesh(cr, numel, mult=mult)

    if log:
        maxv = 10*np.log10(arr2use).max()
        if logrange is not None:
            vmaxl = maxv
            vminl = maxv - logrange

        im = plt.imshow(10*np.log10(arr2use), cmap=cmapl, aspect=aspect,
            norm=plt.Normalize(vmin=vminl, vmax=vmaxl, clip=True),
            extent=extent)

        cbar = plt.colorbar(im)
        cbar.ax.set_xlabel(clab) if clab is not None else cbar.ax.set_xlabel('dB')
        # cbar = plt.colorbar(cmap=cmapl)
        # cbar = plt.colorbar(sm, boundaries=[vminl, vmaxl],
        #     ticks=[-70, -60, -50, -40, -30, -20, -10])

    else:

        if vmin is None:
            vmin = np.nanmin(arr2use)

        if vmax is None:
            vmax = np.nanmax(arr2use)

        im = plt.imshow(arr2use, cmap=cmap, aspect=aspect, vmin=vmin, vmax=vmax,
            extent=extent, alpha=1.0)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax)

        plt.sca(ax)

        if clab is not None and len(clab) > 7:
            cbar.ax.set_ylabel(clab)

        else:
            cbar.ax.set_xlabel(clab) if clab is not None else cbar.ax.set_xlabel('Power')

    if add_contour:

        x, y = bt.get_lin(cr, numel, mult=mult)
        if log:

            if levels is None:
                levels = [maxv-20, maxv-13, maxv-10, maxv-3]
            else:
                levels = levels
                if relative_levels:
                    levels += maxv

            c = plt.gca().contour(x, y, 10*np.log10(np.flipud(arr2use)),
                cmap=cmap_grey, linewidths=lw, levels=levels)

            if contlabel:
                plt.clabel(c, fmt="%5.0f")
        else:
            c = plt.gca().contour(x, y, np.flipud(arr2use), cmap=cmap_grey,
                 linewidths=0.5, origin='lower', extent=extent)

            if contlabel:
                plt.clabel(c, fmt="%5.2f")

    if plot_pix:
        plot_pixels()

    if plot_circ:
        plot_circle(plt.gca(), crad, color=circle_color,
            ls=circle_ls)

    if plot_ellip:
        if len(np.shape(ellipse)) == 1:

            if ellipse_pars is not None:
                col2use = ellipse_pars['color']
                ls2use = ellipse_pars['ls']
            else:
                col2use, ls2use = 'white', '-'

            plot_ellipse(plt.gca(), ellipse,
                color=col2use, ls=ls2use)

        else:
            for ei, ellip in enumerate(ellipse):
                if ellipse_pars is not None:
                    col2use = ellipse_pars['color'][ei]
                    ls2use = ellipse_pars['ls'][ei]
                else:
                    col2use, ls2use = 'white', '-'

                plot_ellipse(plt.gca(), ellip,
                    color=col2use, ls=ls2use)

    if plot_squares:

        squares = []
        for xsqi, ysqi in zip(xsq, ysq):
            # square_verts.append(fp.squareVert([xsqi, ysqi], square_width))
            square = patches.Rectangle([xsqi, ysqi], square_width, square_width,
                edgecolor='white', linewidth=2, facecolor='None', alpha=0.2)

            squares.append(square)
            ax.add_artist(square)

        N = len(xsq)
        square_mask = merge_masks(map(mask_patch,
            squares, itertools.repeat(xx, N), itertools.repeat(yy, N)))

    configure_labels(title=title, xlab=xlab, ylab=ylab)

    if anno_str is not None:
        plt.annotate(anno_str, xy=(0.03, 0.95), xycoords='axes fraction',
            color='white')

    fname2use = 'co' if fname is None else fname
    if log:
        fname2use = fname+'_log'

    if scatter_plot:

        if log:
            plt.scatter(xs, ys, c=scatter_arr, s=50, marker='h',
                cmap=cmapl, edgecolor='black',
                vmin=np.nanmin(10*np.log10(arr2use)),
                vmax=np.nanmax(10*np.log10(arr2use)))
        else:
            plt.scatter(xs, ys, c=scatter_arr, s=50, marker='h',
                cmap=cmap, edgecolor='black', vmin=vmin, vmax=vmax)

    if add_grids:
        ax = plt.gca()

        ylim = ax.get_ylim()
        xlim = ax.get_xlim()
        ax.grid(color='w', linestyle='-', linewidth=1, alpha=0.2)
        ax.set_xticks(np.arange(np.round(xlim[0], decimals=1),
            np.round(xlim[1], decimals=1), 0.1), minor=True)

        ax.set_xticks(np.arange(np.round(ylim[0], decimals=1),
            np.round(ylim[1], decimals=1), 0.1), minor=True)

        #ax.set_yticks(np.arange(-3, 3.1, 0.1), minor=True)
        ax.grid(which='minor', color='w', linestyle='-', linewidth=0.1, alpha=0.3)

    if save:
        plt.savefig(imgd+fname2use+'.png', dpi=dpi, bbox_inches='tight')

    plt.close()