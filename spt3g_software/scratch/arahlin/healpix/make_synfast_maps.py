import numpy as np
import healpy as hp
import healpix_tools as hpt
import os
# from spt3g import core

def mask_1500d(nside, coord=None):
    return hpt.latlon_mask(
        nside,
        latrange=(-74, -38),
        lonrange=(-68, 68),
        coord=coord,
    )

camb_file = os.path.expandvars('$SPT3G_SOFTWARE_PATH/simulations/camb/'
                               'base_plikHM_TTTEEE_lowl_lowE_lensing_lensedCls.dat')
cls = np.loadtxt(camb_file, unpack=True)
ell, cls = cls[0], cls[1:]
ellfac = ell * (ell + 1.) / 2. / np.pi
cls /= ellfac

nside = 8192
# nside = 512
lmax = min([3 * nside - 1, len(ell) - 1])
mmax = lmax
verbose = True
pixwin = False

units = 0.001 # core.G3Units.uK  # ensure maps are stored with correct unit information
cls *= units * units

fileroot = '/spt/user/arahlin/synfast'
filefmt = 'cmb_lmax{:d}_nside{:d}_sim{{:04d}}_lensed_{{:s}}.fits'.format(lmax, nside)

fwhm = 0. # no beam

mask = mask_1500d(nside)

for idx in range(500): # range(1):
    print('Sim', idx)

    almfile = os.path.join(fileroot, filefmt.format(idx, 'alm'))
    mapfile = os.path.join(fileroot, filefmt.format(idx, 'map'))
    if os.path.exists(mapfile):
        continue

    # read in alm's
    if not os.path.exists(almfile):
        print('Creating alms')
        # use a fixed seed per sim index so that realization doesn't
        # change even if map parameters do
        alms = hpt.synalm(cls, lmax=lmax, mmax=mmax, new=True, seed=idx)
        hpt.write_alm(almfile, alms)
    else:
        print('Reading in alms')
        alms = tuple(hp.read_alm(almfile, hdu=hdu) for hdu in [1,2,3])

    # TODO:
    # lenspix correction
    # point sources
    # nominal beam

    # create map
    print('Creating map')
    m = np.asarray(hp.alm2map(
        alms,
        nside,
        lmax=lmax,
        mmax=mmax,
        pixwin=pixwin,
        verbose=verbose,
        pol=True,
        fwhm=fwhm, # for a real beam, use an almxfl call before this step instead
    ))
    del alms
    print(m.max(axis=-1), m.min(axis=-1))

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
