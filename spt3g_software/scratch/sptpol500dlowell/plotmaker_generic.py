'''
This just spits out T,Q,U maps and power spectra. 

run it as 'python looping_plotmaker_generic.py [filename] [short phrase for titling output files, i.e. name of observation or some short string]'

WARNING: assumed 10 arcminute pixels for the ell axis at the end, and that the second frame in a .g3 map file will be the 150s data.

'''
import argparse as ap
from spt3g import core, std_processing, mapmaker, dfmux, calibration, xtalk, todfilter, coordinateutils
from spt3g.timestreamflagging.glitchfinding import get_num_glitches
from spt3g.pointing import offline_pointing as op
from spt3g.std_processing import weighting
import scipy.stats
import pdb
import os, numpy
import pickle
from spt3g.coordinateutils import FlatSkyMap
import copy
from pylab import *
import sys

map_unfiltered1 = sys.argv[1]

m_uf1 = list(core.G3File(map_unfiltered1))

m_uf_use1 = m_uf1[1] #Note this assumes that the second map in the frame is going to be the 150s data, which is consistent with the current committed mapmaking script

from spt3g.coordinateutils import FlatSkyMap
import copy

def change_apod_mask():
    bl = pickle.load(open('/home/javva/spt3g/spt3g_software/scratch/javva/iceland/ptsrc_mask_10.0arcmin.pkl'))
    bl = numpy.ascontiguousarray(bl[::-1,::-1])
    old_map = FlatSkyMap(bl, res = 0.0001454441043328608, proj =coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea, alpha_center = 0, delta_center = -59.0333*core.G3Units.degrees, pol_type = core.MapPolType.T,coord_ref = core.MapCoordReference.Equatorial)
    new_map = copy.copy(m_uf_use1['T'])
    mapmaker.rebin_map(old_map,new_map)
    return new_map

apod_map = change_apod_mask()


#Apodize Maps 
def load_files(d):
    map_unfiltered = d
    m_uf  = list(core.G3File(map_unfiltered))

    m_uf_use = m_uf[1]

    return m_uf_use


def remove_w(fr):
    return mapmaker.mapmakerutils.remove_weight(fr['T'],fr['Q'],fr['U'],fr['Wpol'])

def return_ps(fr):
    apodized_unfiltered_T = (fr)*apod_map

    new_thing_uf = numpy.nan_to_num(np.asarray(apodized_unfiltered_T))

    two_d_fft_uf =numpy.fft.fft2(new_thing_uf)

    two_d_ps_uf = abs(two_d_fft_uf**2)

    reso_rad = 10*core.G3Units.arcmin/core.G3Units.radians

    el_x, ell_y = np.meshgrid(2*np.pi*np.fft.fftfreq(two_d_ps_uf.shape[1], reso_rad),2*np.pi*np.fft.fftfreq(two_d_ps_uf.shape[0\
], reso_rad))

    ell_x = el_x[0][0:len(el_x[0])/2]
    one_d_ps_uf = two_d_ps_uf[two_d_ps_uf.shape[0]/2,:two_d_ps_uf.shape[1]/2]

    avg_one_d_ps_uf = np.average(two_d_ps_uf[:,:two_d_ps_uf.shape[1]/2], axis = 0)
    return avg_one_d_ps_uf, ell_x


d = sys.argv[2]
m_uf_use = load_files(sys.argv[1])
rw_uf = remove_w(m_uf_use)
imshow(rw_uf[0]*apod_map, vmin = -1e-5, vmax = 1e-5)
plt.title('T Map, %s'%d)
plt.savefig('data_maps_T_%s.eps'%d)
clf()

imshow(rw_uf[1]*apod_map, vmin = -2e-5, vmax = 2e-5)
plt.title('Q Map, %s'%d)
plt.savefig('data_maps_Q_%s.eps'%d)                    
clf()

imshow(rw_uf[2]*apod_map, vmin = -1e-5, vmax = 1e-5)
plt.title('U Map, %s'%d)
plt.savefig('data_maps_unfiltered_U_%s.eps'%d)                    
clf()
print 'Taking PSD'
uf_t_ps, ell_x = return_ps(rw_uf[0])
plt.plot(ell_x,uf_t_ps,'r')
plt.semilogx()
plt.semilogy()
plt.ylabel('PSD')
plt.xlabel(' Either ell or ell/2pi')
plt.title('TT Power Spectrum %s'%d)
plt.savefig('data_maps_TT_%s.eps'%d)
plt.clf()

uf_t_ps, ell_x = return_ps(rw_uf[1])
plt.plot(ell_x,uf_t_ps,'r')
plt.semilogx()
plt.semilogy()
plt.title('QQ Power Spectrum %s'%d)
plt.ylabel('PSD')
plt.xlabel(' Either ell or ell/2pi')
plt.savefig('data_maps_QQ_%s.eps'%d)
plt.clf()

uf_t_ps, ell_x = return_ps(rw_uf[2])
plt.plot(ell_x,uf_t_ps, 'r')
plt.semilogx()
plt.semilogy()
plt.ylabel('PSD')
plt.xlabel(' Either ell or ell/2pi')
plt.title('UU Power Spectrum %s'%d)
plt.savefig('data_maps_UU_%s.eps'%d)
plt.clf()
