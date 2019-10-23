import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.optimize as opt
from spt3g import core, mapmaker, util

"""
This script plots Mars maps from disk and fits a 2-d gaussian beam
profile to them.
"""
show_plots = True
save_plots = True

plt.style.use('bmh')

# These are for identifying the name of the maps on disk
poly_order = 3
mars_offs = False
flag_sat = False

if mars_offs:
    mtag = '_marsoffs'
else:
    mtag = ''
if flag_sat:
    ftag = '_flagsat'
else:
    ftag = ''

map_dir = os.path.join('/spt', 'user', 'agambrel', 'beams', 
                       '201903', 
                       'mars_maps_source_centered_poly{}'.format(poly_order))
file_name = 'mars_latest_yhack_{tag}.g3'
fig_dir = os.path.join('/spt', 'user', 'agambrel', 'beams', '201903', 'plots',
                       'source_centered')
if not os.path.exists(fig_dir):
    os.mkdir(fig_dir)

tags = ['54868035', '55211360']

for tag in tags:
    fl = core.G3File(
        os.path.join(map_dir, file_name.format(tag=tag)))
    fig, ax = plt.subplots(1, 3, figsize=(11,4))
    fig.suptitle(tag)
    i = 0
    for fr in fl:
        print(fr)
        if 'T' in fr and fr['Id'] not in ['PointSourceMask', 'bsmap', 
                                          'PointSourceMap']:
            Tmap = np.asarray(
                mapmaker.mapmakerutils.remove_weight_t(fr['T'], fr['Wunpol']))
            res = fr['T'].res / core.G3Units.rad
            sidelength = Tmap.shape[-1]
            ax[i].matshow(Tmap, cmap='inferno', 
                          interpolation='none')
            ax[i].set_title('{}'.format(fr['Id'].split('-')[-1]))
            
            params = util.maths.fit_gaussian2d(Tmap)
            (amp, x, y, width_x, width_y, rota, height) = params
            fit = util.maths.gaussian2d(*params)
            ax[i].contour(fit(*np.indices(Tmap.shape)), 
                          levels=[amp/10, amp/2],
                          cmap=plt.cm.winter)
            
            fac = 2 * np.sqrt(2 * np.log(2))
            width_x = np.abs(width_x)
            width_y = np.abs(width_y)
            fwhm_pix = np.sqrt(((fac*width_x)**2. + (fac*width_y)**2)/2.)
            fwhm_arcmin = np.degrees(fwhm_pix*res)*60.
            fwhm_x = 60 * np.degrees(fac*width_x*res)
            fwhm_y = 60 * np.degrees(fac*width_y*res)
            ell = np.sqrt(np.abs(width_x**2-width_y**2) /
                          np.max([width_x**2, width_y**2])) # ellipticity
            a = np.max([width_x, width_y])
            b = np.min([width_x, width_y])
            ar = (a - b) / a # aspect ratio
            ax[i].text(0.05, 0.05, 
                     "({:0.1f},{:0.1f})\nfwhm: {:0.3f}'\nell: {:0.3f}".format(
                         np.degrees((x-sidelength/2.)*res)*60,
                         np.degrees((y-sidelength/2.)*res)*60,fwhm_arcmin, ell),
                     verticalalignment='bottom', horizontalalignment='left',
                     transform=ax[i].transAxes, color='white')
            ax[i].text(0.95, 0.05, 
                     "fwhm_x: {:0.3f}'\nfwhm_y: {:0.3f}'".format(fwhm_x, 
                                                                fwhm_y),
                     verticalalignment='bottom', horizontalalignment='right',
                     transform=ax[i].transAxes, color='white')

            ax[i].set_xlim(int(x-15),int(x+15))
            ax[i].set_ylim(int(y-15),int(y+15)) 
            xlocs = ax[i].get_xticks()
            ylocs = ax[i].get_yticks()
            ax[i].set_xticklabels(['{:.0f}'.format(x) for x in 
                        60*np.degrees(((xlocs-xlocs[0])*res))])
            ax[i].set_yticklabels(['{:.0f}'.format(x) for x in 
                        60*np.degrees(((ylocs-ylocs[0])*res))])
            ax[i].plot(x, y, '+', color='k')
            
            ax[i].set_xlabel('arcmin')
            ax[i].set_ylabel('arcmin')
            i+=1
    if save_plots:
        plt.savefig(
            os.path.join(fig_dir, 
                         '{}_beam_maps_poly{}{}{}.png'.format(tag, poly_order,
                                                              mtag, ftag)))
    if show_plots:
        plt.show()
    plt.close(fig)
