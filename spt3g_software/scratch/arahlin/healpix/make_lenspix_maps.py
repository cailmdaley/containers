import numpy as np
import healpy as hp
import healpix_tools as hpt
import os
from spt3g import core

def mask_1500d(nside, coord=None):
    return hpt.latlon_mask(
        nside,
        latrange=(-74, -38),
        lonrange=(-68, 68),
        coord=coord,
    )

fileroot = '/spt/user/arahlin/lenspix_alms'
outfileroot = fileroot.replace('_alms', '_maps')
filefmt = 'lensed_cmb_lmax7000_nside8192_interp0.3_method1_pol_1_sim_{:d}_lensed_{:s}.fits'

nside = 8192
pixwin = False
lmax = 7000
mmax = None
verbose = True
units = core.G3Units.uK  # ensure maps are stored with correct unit information

mask = mask_1500d(nside)

for idx in range(500):
    print('Sim', idx)

    mapfile = os.path.join(outfileroot, filefmt.format(idx, 'map'))
    mapfile = mapfile.replace('nside8192', 'nside{:d}'.format(nside))
    if os.path.exists(mapfile):
        continue

    # read in alm's
    print('Reading in alms')
    almfiles = [os.path.join(fileroot, filefmt.format(idx, 'alm' + x)) for x in 'TGC']
    alms = tuple(hp.read_alm(f).astype(np.complex128) for f in almfiles)

    # create map
    print('Creating map')
    m = np.asarray(hp.alm2map(
        alms,
        nside,
        pixwin=pixwin,
        verbose=verbose,
    )) * units
    del alms

    # store map in 3G-compatible format
    print('Writing map')
    hpt.write_map(
        mapfile,
        m,
        coord='C',
        mask=mask,
        column_names=['T', 'Q', 'U'],
        extra_header={'WEIGHTED': False,
                      'UNITS': 'Tcmb'},
    )
    del m
