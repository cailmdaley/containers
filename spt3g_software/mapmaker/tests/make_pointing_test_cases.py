from sptpol_software.observation.sky import pix2Ang, ang2Pix
from spt3g import core, coordinateutils
import random
import numpy as np

un = core.G3Units
projs = coordinateutils.MapProjection

projections = [projs.Proj0,
               projs.Proj1,
               projs.Proj2,
               projs.Proj4,
               projs.Proj5
           ]

n_random = 10000
ra_min = -355
ra_max = 355
dec_min = -80
dec_max = 80
reso_arcmin_min = 0.05
reso_arcmin_max = 2.5
n_x_min = 100
n_x_max = 1000
n_y_min = 500
n_y_max = 1000

out_ds = []

for proj in projections:
    for i in range(n_random):
        ra0  = ra_min + (ra_max-ra_min) * random.random()
        dec0 = dec_min + (dec_max-dec_min) * random.random() 
        n_x = random.randint(n_x_min, n_x_max)
        n_y = random.randint(n_y_min, n_y_max)
        reso_arcmin = (reso_arcmin_max-reso_arcmin_min) * random.random() + reso_arcmin_min
        
        while(True):
            ra  = ra0  + n_x * reso_arcmin/60.0 * ( .5 - random.random())
            dec = dec0 + n_y * reso_arcmin/60.0 * ( .5 - random.random())
            if (dec > 85 or dec < -85) or (ra > 355 or ra < -355):
                continue

            ss_pix = ang2Pix(np.array([ra, dec]),
                             np.array([ra0, dec0]),
                             reso_arcmin,
                             np.array([n_y, n_x]),
                             proj=proj,
                             round=True, bin_center_zero=True,
                             return_validity=True,use_c_code=False)
            pix_x = ss_pix[0][1][0]
            pix_y = ss_pix[0][0][0]

            in_ss_pix = ss_pix[1][0][0] and (pix_x < n_x) and (pix_y < n_y) and (pix_x >= 0 ) and (pix_y >= 0)

            ss_ang = pix2Ang( np.array([[pix_y], [pix_x]]), np.array([ra0, dec0]), reso_arcmin, np.array([n_y, n_x]),  proj )
            ra_out = ss_ang[0][0]
            dec_out = ss_ang[1][0]
            break            
        store_d = {'pix_x':pix_x, 'pix_y':pix_y, 'ra0':ra0, 'dec0':dec0, 'proj':proj,
                   'n_x':n_x, 'n_y':n_y, 'ra_in':ra, 'dec_in': dec, 'reso_arcmin':reso_arcmin,
                   'ra_out':ra_out, 'dec_out':dec_out, 'in_ss_pix':in_ss_pix  }
        out_ds.append(store_d)

import cPickle as pickle
pickle.dump(out_ds, open('test_pointing_information.pkl','w'), protocol=2)
