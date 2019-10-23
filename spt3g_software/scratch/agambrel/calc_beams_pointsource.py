import matplotlib.pyplot as plt
import numpy as np
import os
from spt3g import core
from spt3g.mapspectra import map_analysis as ma
from spt3g.mapspectra import apodmask as am
from itertools import combinations

plt.style.use('bmh')
# full 2018 data coadd
#field_map = '/spt/user/ddutcher/coadds/lpf9000_tonly_0.5res.g3'
# online pointing partial data coadd
field_map = os.path.join('/spt', 'user', 'ddutcher', 'coadds', 
                         'lpf18k_150t_0.5res_online_20190308_51.g3')
ps_list_file = os.path.join('/home', 'ddutcher', 'code', 'spt3g_software',
                            'scratch', 'ddutcher', 
                            '1500d_3band_10sigma_ptsrc.txt')
fig_dir = os.path.join('/spt', 'user', 'agambrel', 'beams', '201903', 'plots')
# for debugging
view_maps = False
# for binning
lmax = 8000
lmin = 120
dl = 125

# The list of nsource brightest sources is crossed
nsource = 10
ps_list = np.loadtxt(ps_list_file)[:nsource]

# half-size of cutout in pixels
pm_size = 100 

#ps_maps = {'90': [], '150': [], '220': []}
#ps_weights = {'90': [], '150': [], '220': []}
#ps_specs = {'90': [], '150': [], '220': []}
# online pointing only has 150 GHz data
ps_maps = {'150': []}
ps_weights = {'150': []}
ps_specs = {'150': []}
for p, ps in enumerate(ps_list):
    fl = core.G3File(field_map)
    ra = ps[1] * core.G3Units.deg
    dec = ps[2] * core.G3Units.deg
    for fr in fl:
        if 'T' in fr and fr['Id'] != 'PointSourceMask':
            res=fr['T'].res
            for freq in ps_maps.keys():
                if freq in fr['Id']:
                    freq0 = freq
            Tmap = np.asarray(fr['T'])
            pix = fr['T'].angle_to_pixel(ra,dec)
            x, y = np.unravel_index(pix, Tmap.shape)
            xsl = slice(x-pm_size, x+pm_size)
            ysl = slice(y-pm_size, y+pm_size)
            for k in fr.keys():
                if k!= 'Id':
                    if k=='Wunpol':
                        wup = np.asarray(fr[k].TT)
                        wup = wup[xsl, ysl]
                        ps_weights[freq0].append(wup)
                    else:
                        Tmap0 = np.asarray(fr[k])
                        Tmap0 = Tmap0[xsl, ysl]
                        if view_maps:
                            plt.imshow(Tmap0/wup, interpolation=None)
                            plt.show()
            ps_maps[freq0].append([Tmap0, wup])

for freq, maps in ps_maps.items():
    print(freq)
    n=0
    apod_mask = None
    for comb in combinations(maps, 2):
        print(n)
        n+=1
        # Cross all the maps
        if apod_mask is None:
            apod_mask = am.make_border_apodization(
                comb[0][1], apod_type='cos', res=res,
                radius_arcmin=0.05*pm_size*res/core.G3Units.arcmin)
        cls = ma.calculateCls(comb[0][0], cross_map=comb[1][0], 
                              weight=comb[0][1], cross_weight=comb[1][1], 
                              res=res, t_only=True, apod_mask=apod_mask, 
                              ell_min=lmin, ell_max=lmax, delta_ell=dl)
        ell_binned = cls['ell']
        # Bl is square root of crosses-- currently unnormalized
        ps_specs[freq].append(cls['TT']**(1/2.))
#np.savez('/spt/user/agambrel/beams/201903/ps_beams_online.npz', **ps_specs)

for freq, specs in ps_specs.items():
    mean_spec = np.mean(specs, axis=0)
    std_spec = np.std(specs, axis=0)
    plt.errorbar(ell_binned, mean_spec, yerr=std_spec/np.sqrt(len(specs)))
    plt.title(freq)
    plt.show()

