import os
import numpy as np
from spt3g.util import healpix_tools as hpt
from spt3g import core

high_pass = 300
high_pass_width = 150

def mask_1500d(nside, coord=None):
    return hpt.latlon_mask(
        nside,
        latrange=(-74, -38),
        lonrange=(-68, 68),
        coord=coord
    )

def hpf(cutoff, width, nside):
    # high pass in l-space
    lmax = 3 * nside - 1
    ell = np.arange(lmax+1)
    trans = np.logical_and(ell > cutoff-width/2, ell < cutoff+width/2)
    kernel = np.zeros_like(ell, dtype=float)
    kernel[ell <= cutoff-width/2] = 0
    kernel[ell >= cutoff+width/2] = 1
    kernel[trans] = 0.5 - 0.5 * np.cos(ell[trans]*np.pi/width - 
                                       np.pi*(cutoff-width/2)/width)
    return ell, kernel
    
mask = mask_1500d(2048)
map_path = '/spt/user/agambrel/planck_maps'
mname = 'HFI_SkyMap_{}_2048_R3.01_halfmission-{}.fits'

for freq in [100, 143, 217]:
    for hm in [1, 2]:
        if high_pass is not None:
            mname_out = 'HFI_SkyMap_{}_2048_R3.01_halfmission-{}_cut_C_G3Units_hpl{}.fits'.format(freq, hm, high_pass)
        else:
            mname_out = 'HFI_SkyMap_{}_2048_R3.01_halfmission-{}_cut_C_G3Units.fits'.format(freq, hm)

        if os.path.exists(os.path.join(map_path, mname_out)):
            continue

        print('{} GHz HM{}'.format(freq, hm))
        m = hpt.read_map(os.path.join(map_path, mname.format(freq, hm)),
                         field=None)
        if high_pass is not None:
            print('Applying high pass filter')
            ell, bl_hp = hpf(high_pass, high_pass_width, 2048)
            m[:3] = hpt.smoothing(m[:3], beam=np.tile(bl_hp, 3))
        print('Rotating from G to C')
        mr = hpt.rotate_map(m[:3], ['G', 'C'])
        mr[:,mask] *= core.G3Units.K
        print('Writing output map')
        hpt.write_map(
            os.path.join(map_path, mname_out.format(freq, hm)),
            mr,
            coord='C',
            mask=mask,
            column_names=['T', 'Q', 'U'],
            extra_header={'WEIGHTED': False,
                          'UNITS': 'Tcmb'},
        )
