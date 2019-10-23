"""
Using a set of Mars maps, get the shape of the beam function from
cross spectra of multiple observations. If desired, correct for extra 
smoothing due to pointing jitter computed from cross spectra of point 
sources in the field. If desired, use cross spectra between SPT and Planck
to get a check on the low ell beam. This is only used for plotting purposes.
Normalize to unity between ell of 600 and 1000. Interpolate between bins 
with a cubic spline. Write product to disk as a text file. Write to disk plots 
that show the beams computed from the three different sources. Defaults are set
to what are used for the final SPT beam product.
"""
import matplotlib.pyplot as plt
import numpy as np
import os
import healpy as hp
import getpass
import scipy.optimize as opt
from scipy import stats
from scipy import interpolate
from itertools import combinations
import argparse as ap
from spt3g import core, coordinateutils
from spt3g.mapspectra import map_analysis as ma
from spt3g.mapspectra import pixel_window
from spt3g.mapspectra import basicmaputils as bmu
from spt3g.mapmaker import mapmakerutils as mmu
from spt3g.util import maths
from spt3g.beams import beam_analysis as ba

P = ap.ArgumentParser(description='Compute beam window functions',
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('--no-do-planck', dest='do_planck', action='store_false',
               default='True', 
               help="Don't use Planck-SPT cross spectra to check low ell beam")
P.add_argument('--no-do-ps', dest='do_ps', action='store_false',
               default='True', 
               help="Don't use point source cross spectra to correct for "
               "pointing jitter or to check high ell beam.")
P.add_argument('--no-compute-transfer', dest='compute_transfer', 
               action='store_false', default='True', 
               help="Don't compute the transfer function from a point-source"
               "sim and correct Mars/point source beams for it.")
P.add_argument('--planet-map-dir', action='store', 
               default=os.path.join('/spt', 'user', 'agambrel', 'beams', 
                                    'beam00_maps'),
               help='Location of planet maps and mock planet maps')
P.add_argument('--field-map-dir', action='store', 
               default=os.path.join('/spt', 'user', 'agambrel', 'beams', 
                                    'beam00_maps'),
               help='Location of mock observed Planck maps and field maps')
P.add_argument('--ps-map-dir', action='store', 
               default=os.path.join('/spt', 'user', 'agambrel', 'beams', 
                                    'beam00_maps'),
               help='Location of point source maps')
P.add_argument('--planet-obsids', action='store', nargs='+',
               default=['54868035', '55211360', '55383550', '55643104', 
                        '56489370'],
               help='Obsids of maps to use for cross spectra')
P.add_argument('--planet-map-name', action='store', default='mars_poly3_{}.g3',
               help='Name of the map files with obsid replaced with brackets')
P.add_argument('--planet-sim-name', action='store', 
               default='mars_poly3_sim_{}.g3',
               help='Name of the sim map files with obsid replaced with '
               'brackets')
P.add_argument('--output-path', action='store', 
               default=os.path.join(os.getenv('SPT3G_SOFTWARE_PATH'), 'beams', 
                                    'products'),
               help='Location for writing out beam00.txt file')
P.add_argument('--fig-dir', action='store',
               default=os.path.join('/spt', 'user', getpass.getuser(),
                                    'figures'),
               help='Directory in which to write figures')
P.add_argument('--force-recompute', action='store_true', default=False,
               help='Recompute cross spectra even if they have been computed '
               'before and can be read from disk.')
P.add_argument('--ptsrc-mask', action='store', 
               default=os.path.join(os.getenv('SPT3G_SOFTWARE_PATH'), 
                                    'beams', 'mapmaking_scripts', 
                                    '1500d_3band_10sigma_ptsrc.txt'),
               help='List of point sources to mask')
args = P.parse_args()

### Set up constants ###
point_sources = range(10)

if not os.path.exists(args.fig_dir):
    os.mkdir(args.fig_dir)

lmin = 50
lmax = 20000
dl = 50

ells = np.arange(lmin, lmax+1, dl)
bins = bmu.get_reg_spaced_ell_bins(ells)
ell_bins = [x[0] for x in bins]
ell_bins.append(bins[-1][-1])
########################

headers = [
    'beam00\n May 2019 \n For 2018 SPT-3G Data \n'
    'Beam product computed from the last 5 Mars observations of 2018. '
    'Mars maps are made in source-centered coordinates. Spectra are '
    'corrected by the 15 arcsecond size of Mars, the 0.2 arcmin pixel window '
    'function, and the transfer function for the third order polynomial '
    'filter computed from mock observations. This beam is convolved with'
    ' a 42 arcsecond Gaussian to account for pointing jitter in the '
    'field observations as determined from taking the cross spectra '
    'of the ten brightest point sources in the 2018 field observations '
    'of ID 47* to 57*. The jitter beam calculation is found by dividing '
    'the mean point source beam by the mean Mars beam, and fitting a '
    'Gaussian. The mean Gaussian across the three frequencies is used. '
    'The per-ell beam is a cubic interpolation between points from '
    'ell={} to ell={} in bin size={}'.format(lmin, lmax, dl) + 'It is '
    'zero above ell={}'.format(lmax) + ' and the constant at the lowest bin '
    'value below ell={}.'.format(lmin) + 'The beam is normalized to 1 '
    'between ell=600 and ell=1000.\n']

plt.style.use('bmh')


#### compute all the pixel window function needed ####
pw_mars = pixel_window.pixwin_flatsky(range(lmax+1), 0.2)
pw_mars_binned = stats.binned_statistic(np.arange(len(pw_mars)), 
                                        pw_mars, bins=ell_bins)[0]

# final map pixwin, source map pixwin, source map smoothing
pw_mars_sim = (pixel_window.pixwin_flatsky(range(lmax+1), 0.2) * 
               pixel_window.pixwin_flatsky(range(lmax+1), 0.1) * 
               hp.gauss_beam(np.deg2rad(1./60), lmax=lmax))
pw_sim_binned = stats.binned_statistic(np.arange(len(pw_mars_sim)), 
                                       pw_mars_sim, bins=ell_bins)[0]

# just input pixel window function -- output cancels in ratio with SPT

pw_ps = pixel_window.pixwin_flatsky(range(lmax+1), 0.1)
pw_ps_binned = stats.binned_statistic(np.arange(len(pw_ps)), 
                                      pw_ps, bins=ell_bins)[0]

########################################################

#### Set up data structures ####
freqs = [90, 150, 220]

colors = {90: ['#084594', '#6baed6', '#c6dbef'],
          150: ['#7a0177', '#dd3497', '#f768a1'],
          220: ['#005a32', '#41ab5d', '#a1d99b']}

starting_maps = np.zeros([len(freqs), len(args.planet_obsids)], dtype='object')

mars_maps = dict(zip(freqs, starting_maps.copy()))
sim_maps = dict(zip(freqs, starting_maps.copy()))
weight_maps = dict(zip(freqs, starting_maps.copy()))
sim_weight_maps = dict(zip(freqs, starting_maps.copy()))
mars_xspecs = {}
xspecs_sim = {}
f_l = {}

if args.do_ps:
    starting_maps_ps = np.zeros([len(freqs), len(point_sources)], 
                                dtype='object')
    ps_maps = dict(zip(freqs, starting_maps_ps.copy()))
    ps_weight_maps = dict(zip(freqs, starting_maps_ps.copy()))
    ps_xspecs = {}

##################################################################

### SPT x Planck ###
if args.do_planck:
    print('Getting beam from SPT x Planck')
    files_exist = (
        os.path.exists(os.path.join(args.field_map_dir, 'field_specs.npy')) &
        os.path.exists(os.path.join(args.field_map_dir, 'field_specs_lo.npy')) &
        os.path.exists(os.path.join(args.field_map_dir, 'field_specs_hi.npy')))

    if not args.force_recompute and files_exist:
        field_xspecs = np.load(os.path.join(args.field_map_dir, 
                                            'field_specs.npy'))
        field_xspecs = field_xspecs.item()
        field_xspecs_lo = np.load(os.path.join(args.field_map_dir, 
                                               'field_specs_lo.npy'))
        field_xspecs_lo = field_xspecs_lo.item()
        field_xspecs_hi = np.load(os.path.join(args.field_map_dir, 
                                               'field_specs_hi.npy'))
        field_xspecs_hi = field_xspecs_hi.item()
    else:
        field_file = os.path.join(args.field_map_dir, 'spt_field_lr.g3')
        field_maps_l, field_mask, field_mask_lo, field_mask_hi = \
                ba.get_field_maps(field_file, tag='Left',
                                  return_mask=True,
                                  return_mask_split_dec=True)
        field_maps_r = ba.get_field_maps(field_file, tag='Right')
        # point source mask
        ptsrc = ma.apodmask.make_apodized_ptsrc_mask(
            field_maps_l[150], args.ptsrc_mask, radius_arcmin=10.,
            zero_border_arcmin=0.)
        field_mask *= ptsrc
        field_mask_hi *= ptsrc
        field_mask_lo *= ptsrc

        hm1_file = os.path.join(args.field_map_dir, 
                                'hm1_lr.g3')
        hm1_maps_l = ba.get_field_maps(hm1_file, tag='Left')
        hm1_maps_r = ba.get_field_maps(hm1_file, tag='Right')
        hm2_file = os.path.join(args.field_map_dir, 
                                'hm2_lr.g3')
        hm2_maps_l = ba.get_field_maps(hm2_file, tag='Left')
        hm2_maps_r = ba.get_field_maps(hm2_file, tag='Right')

        field_xspecs = ba.get_field_xspecs(
            field_maps_l, field_maps_r, hm1_maps_l, hm1_maps_r, hm2_maps_l, 
            hm2_maps_r, field_mask, lmin=lmin, lmax=lmax, dl=dl, 
            freqs=freqs, 
            out_file=os.path.join(args.field_map_dir, 'field_specs.npy'))
        field_xspecs_lo = ba.get_field_xspecs(
            field_maps_l, field_maps_r, hm1_maps_l, hm1_maps_r, hm2_maps_l, 
            hm2_maps_r, field_mask_lo, lmin=lmin, lmax=lmax, dl=dl, 
            freqs=freqs, 
            out_file=os.path.join(args.field_map_dir, 'field_specs_lo.npy'))
        field_xspecs_hi = ba.get_field_xspecs(
            field_maps_l, field_maps_r, hm1_maps_l, hm1_maps_r, hm2_maps_l, 
            hm2_maps_r, field_mask_hi, lmin=lmin, lmax=lmax, dl=dl, 
            freqs=freqs, 
            out_file=os.path.join(args.field_map_dir, 'field_specs_hi.npy'))
    
### Point source transfer function ###
if args.compute_transfer:
    print('Getting F_l for point source/planet obs')

    if not args.force_recompute and \
            os.path.exists(os.path.join(args.planet_map_dir, 
                                        'mars_sim_specs.npy')):
        xspecs_sim = np.load(os.path.join(args.planet_map_dir, 
                                          'mars_sim_specs.npy'))
        xspecs_sim = xspecs_sim.item()
    else:
        min_width_x = 1e9
        min_width_y = 1e9
        for t, tag in enumerate(args.planet_obsids):
            fl = core.G3File(os.path.join(args.planet_map_dir, 
                                          args.planet_sim_name.format(tag)))
            for fr in fl:
                if 'T' in fr and fr['Id'] not in ['PointSourceMask', 
                                                  'PointSourceMap',
                                                  'bsmap']:
                    d = ba.recenter_maps(fr, sim_maps, sim_weight_maps, t, 
                                         min_width_x, min_width_y)
                    sim_maps, sim_weight_maps, min_width_x, min_width_y = d
                    res_mars = fr['T'].res

        for t, tag in enumerate(args.planet_obsids):
            for j, f in enumerate(freqs):
                # reshape each map to have size of smallest map
                sim_maps[f][t] = maths.recenter_map(
                    sim_maps[f][t], xcenter=sim_maps[f][t].shape[0]/2,
                    ycenter=sim_maps[f][t].shape[1]/2, widthx=min_width_x, 
                    widthy=min_width_y)

                sim_weight_maps[f][t] = maths.recenter_map(
                    sim_weight_maps[f][t], 
                    xcenter=sim_weight_maps[f][t].shape[0]/2,
                    ycenter=sim_weight_maps[f][t].shape[1]/2, 
                    widthx=min_width_x, widthy=min_width_y)

        mean_sim_weight_map = np.mean(np.asarray(
            [wmap for freqs, wmap in sim_weight_maps.items()]).ravel(),
                                      axis=0)
        mask_sim_mars = ma.apodmask.make_border_apodization(
            mean_sim_weight_map, res=res_mars, radius_arcmin=20.,
            smooth_weights_arcmin=60, apod_threshold=0.05)

        for j, f in enumerate(freqs):
            # for each frequency
            for i, m in enumerate(combinations(sim_maps[f], 2)):
                # for each possible pair to cross (excluding autos)
                ell_sim, cls = bmu.simple_cls(m[0], m[1], 
                                              apod_mask=mask_sim_mars, 
                                              ell_min=lmin, ell_max=lmax, 
                                              delta_ell=dl, res=res_mars)

                # correct for pixel window function and beam of source, leaving
                # just filter transfer function
                spec = cls / pw_sim_binned**2
                spec /= np.max(spec)
                if f not in xspecs_sim:
                    xspecs_sim[f] = spec
                else:
                    xspecs_sim[f] = np.vstack([xspecs_sim[f], spec])
        np.save(os.path.join(args.planet_map_dir, 'mars_sim_specs.npy'), 
                xspecs_sim)

for i, f in enumerate(freqs):
    if args.compute_transfer:
        f_l[f] = np.nanmean(xspecs_sim[f], axis=0)
    else:
        f_l[f] = np.ones(len(pw_mars_binned))

### Field point sources ###
if args.do_ps:
    print('Getting beams from point sources')
    min_width_x = 1e9
    min_width_y = 1e9
    
    if not args.force_recompute and \
            os.path.exists(os.path.join(args.ps_map_dir, 'ps_specs.npy')):
        ps_xspecs = np.load(os.path.join(args.ps_map_dir, 'ps_specs.npy'))
        ps_xspecs = ps_xspecs.item()
    else:
        for ps in point_sources:
            print('Point source {}'.format(ps))
            fl = core.G3File(
                os.path.join(args.ps_map_dir, 'map_ps{}.g3'.format(ps)))
            for fr in fl:
                if 'T' in fr and fr['Id'] not in ['PointSourceMask', 
                                                  'PointSourceMap',
                                                  'bsmap']:
                    d = ba.recenter_maps(fr, ps_maps, ps_weight_maps, ps, 
                                         min_width_x, min_width_y)
                    ps_maps, ps_weight_maps, min_width_x, min_width_y = d
                    res_ps = fr['T'].res

        for ps in point_sources:
            for j, f in enumerate(freqs):
                # reshape each map to have size of smallest map
                ps_maps[f][ps] = maths.recenter_map(
                    ps_maps[f][ps], xcenter=ps_maps[f][ps].shape[0]/2,
                    ycenter=ps_maps[f][ps].shape[1]/2, widthx=min_width_x, 
                    widthy=min_width_y)

                ps_weight_maps[f][ps] = maths.recenter_map(
                    ps_weight_maps[f][ps], 
                    xcenter=ps_weight_maps[f][ps].shape[0]/2,
                    ycenter=ps_weight_maps[f][ps].shape[1]/2,
                    widthx=min_width_x, widthy=min_width_y)

        mean_weight_map_ps = np.mean(np.asarray(
            [wmap for freqs, wmap in ps_weight_maps.items()]).ravel(),
                                     axis=0)
        mask_ps = ma.apodmask.make_border_apodization(mean_weight_map_ps, 
                                                    res=res_ps,
                                                    radius_arcmin=20.,
                                                    smooth_weights_arcmin=60,
                                                    apod_threshold=0.05)
        for j, f in enumerate(freqs):
            # for each frequency
            for i, m in enumerate(combinations(ps_maps[f], 2)):
                # for each possible pair to cross, except autos
                ell_ps, cls = bmu.simple_cls(m[0], m[1], apod_mask=mask_ps, 
                                             ell_min=lmin, ell_max=lmax, 
                                             delta_ell=dl, res=res_ps)

                # correct for pixel window function and transfer function
                spec = cls / pw_ps_binned**2. / f_l[f]
                if f not in ps_xspecs:
                    ps_xspecs[f] = spec
                else:
                    ps_xspecs[f] = np.vstack([ps_xspecs[f], spec])
        np.save(os.path.join(args.ps_map_dir, 'ps_specs.npy'), ps_xspecs)

### Mars observations ###
print('Getting beams from Mars')
if not args.force_recompute and \
        os.path.exists(os.path.join(args.planet_map_dir, 'mars_specs.npy')):
    mars_xspecs = np.load(os.path.join(args.planet_map_dir, 
                                       'mars_specs.npy'))
    mars_xspecs = mars_xspecs.item()
else:
    min_width_x = 1e9
    min_width_y = 1e9
    for t, tag in enumerate(args.planet_obsids):
        print('Mars map {}'.format(t))
        fl = core.G3File(
            os.path.join(args.planet_map_dir, args.planet_map_name.format(tag)))
        for fr in fl:
            if 'T' in fr and fr['Id'] not in ['PointSourceMask', 
                                              'PointSourceMap',
                                              'bsmap']:
                d = ba.recenter_maps(fr, mars_maps, weight_maps, t, 
                                     min_width_x, min_width_y)
                mars_maps, weight_maps, min_width_x, min_width_y = d

                res_mars = fr['T'].res

    for t, tag in enumerate(args.planet_obsids):
        for j, f in enumerate(freqs):
            # reshape each map to have size of smallest map
            mars_maps[f][t] = maths.recenter_map(
                mars_maps[f][t], xcenter=mars_maps[f][t].shape[0]/2,
                ycenter=mars_maps[f][t].shape[1]/2, widthx=min_width_x, 
                widthy=min_width_y)
            weight_maps[f][t] = maths.recenter_map(
                weight_maps[f][t], xcenter=weight_maps[f][t].shape[0]/2,
                ycenter=weight_maps[f][t].shape[1]/2, widthx=min_width_x, 
                widthy=min_width_y)

    mean_weight_map = np.mean(np.asarray(
        [wmap for freqs, wmap in weight_maps.items()]).ravel(),
                              axis=0)
    mask_mars = ma.apodmask.make_border_apodization(mean_weight_map, res=res_mars,
                                                  radius_arcmin=20.,
                                                  smooth_weights_arcmin=60,
                                                  apod_threshold=0.05)
    for j, f in enumerate(freqs):
        # for each frequency
        for i, m in enumerate(combinations(mars_maps[f], 2)):
            # for each possible pair to cross, except autos
            ell_mars, cls = bmu.simple_cls(m[0], m[1], 
                                           apod_mask=mask_mars, ell_min=lmin, 
                                           ell_max=lmax, delta_ell=dl, 
                                           res=res_mars)

            # correct for size of Mars, pixel window function, transfer function
            Bm = ba.Bl_planet(ell_mars, size=15.)
            spec = cls / Bm**2. / pw_mars_binned**2. / f_l[f]

            if f not in mars_xspecs:
                mars_xspecs[f] = spec
            else:
                mars_xspecs[f] = np.vstack([mars_xspecs[f], spec])
    np.save(os.path.join(args.planet_map_dir, 'mars_specs.npy'), mars_xspecs)

fig0, ax0 = plt.subplots(1, 1)
fig_highl, ax_highl = plt.subplots(1, 1)

headers.append('ell')
all_ells = np.arange(lmax)
Bl_all_ells = all_ells

jitter_beam = {}
mean_mars_beam = {}
err_mars_beam = {}

for i, f in enumerate(freqs):
    # normalization as in Henning et al
    l600_idx = np.argmin(np.abs(ells-600))
    l1000_idx = np.argmin(np.abs(ells-1000))
    for j, spec in enumerate(mars_xspecs[f]):
        # beam is sqrt of power spectrum
        mars_xspecs[f][j] = np.sqrt(spec)
        mars_xspecs[f][j] /= np.nanmean(mars_xspecs[f][j][l600_idx:l1000_idx])
    mean_mars_beam[f] = np.mean(mars_xspecs[f], axis=0)
    err_mars_beam[f] = np.std(mars_xspecs[f], axis=0) / \
                       np.sqrt(len(mars_xspecs[f]))

    # also high ell normalization to match point sources
    l2000_idx = np.argmin(np.abs(ells-2000))
    l4000_idx = np.argmin(np.abs(ells-4000))
    mean_mars_highl = np.nanmean(mean_mars_beam[f][l2000_idx:l4000_idx])

    mean_mars_beam[f][mean_mars_beam[f]==0] = np.nan
    good_inds = ~np.isnan(mean_mars_beam[f])
    popt, pcov = ba.fit_gauss_lspace(ells[good_inds], 
                                     mean_mars_beam[f][good_inds])
    mars_width = np.deg2rad(np.abs(popt[0]))/(2.*np.sqrt(2*np.log(2.)))

    if args.do_ps:
        l2000_idx = np.argmin(np.abs(ells-2000))
        l4000_idx = np.argmin(np.abs(ells-4000))

        ax_highl.errorbar(ells[l2000_idx:], 
                          mean_mars_beam[f][l2000_idx:], 
                          yerr=err_mars_beam[f][l2000_idx:], 
                          linestyle='-', color=colors[f][0],
                          label='{} Mars fwhm={:.1f} arcmin'.format(
                              f, np.abs(popt[0]*60)))

        for j, spec in enumerate(ps_xspecs[f]):
            # beam is sqrt of power spectrum
            ps_xspecs[f][j] = np.sqrt(spec)
            # normalize each cross spectrum
            high_l_ps = np.nanmean(ps_xspecs[f][j][l2000_idx:l4000_idx])
            ps_xspecs[f][j] *= mean_mars_highl / high_l_ps
        mean_ps_beam = np.mean(ps_xspecs[f], axis=0)
        err_ps_beam = np.std(ps_xspecs[f], axis=0)/np.sqrt(len(ps_xspecs[f]))
        good_inds = ~np.isnan(mean_ps_beam[l2000_idx:])
        
        popt_ps, pcov_ps = ba.fit_gauss_lspace(
            ells[l2000_idx:][good_inds], mean_ps_beam[l2000_idx:][good_inds])
        ps_width = np.deg2rad(np.abs(popt_ps[0]))/(2.*np.sqrt(2*np.log(2.)))
        jitter_width = np.sqrt(ps_width**2-mars_width**2)
        jitter_width_deg = np.rad2deg(jitter_width*2*np.sqrt(2*np.log(2.)))
        jitter_beam[f] = ba.Bl(ells, jitter_width_deg, 1)
        
        ax0.errorbar(ells[l2000_idx:], mean_ps_beam[l2000_idx:], 
                     yerr=err_ps_beam[l2000_idx:], alpha=0.8,
                     color=colors[f][2],
                     label='{} PS fwhm={:.1f} arcmin'.format(
                         f, np.abs(popt_ps[0]*60)))
        ax_highl.errorbar(ells[l2000_idx:], mean_ps_beam[l2000_idx:], 
                          yerr=err_ps_beam[l2000_idx:], alpha=0.8, 
                          color=colors[f][2],
                          label='{} PS fwhm={:.1f} arcmin'.format(
                              f, np.abs(popt_ps[0]*60)))

if args.do_ps:
    # call jitter beam the average over all three frequencies
    mean_jitter = np.mean([jb for jb in jitter_beam.values()], axis=0)
    good_inds = ~np.isnan(mean_jitter[l2000_idx:])
    popt_j, pcov_j = ba.fit_gauss_lspace(ells[l2000_idx:][good_inds], 
                                         mean_jitter[l2000_idx:][good_inds])
    print('jitter fwhm: {} arcsec'.format(popt_j[0]*60*60))
    for i, f in enumerate(freqs):
        # correct mars beam for it
        mean_mars_beam[f] *= mean_jitter
        err_mars_beam[f] *= mean_jitter
        ax_highl.errorbar(ells[l2000_idx:], 
                          mean_mars_beam[f][l2000_idx:], 
                          yerr=err_mars_beam[f][l2000_idx:], 
                          linestyle='--', color=colors[f][1],
                          label='{} Mars with jitter'.format(f))
    
for i, f in enumerate(freqs):
    fig_lowl, ax_lowl = plt.subplots(1, 1)
    fig_lowl_lo, ax_lowl_lo = plt.subplots(1, 1)
    fig_lowl_hi, ax_lowl_hi = plt.subplots(1, 1)

    # do a cubic interpolation to get ell by ell beam
    good_inds = ~np.isnan(mean_mars_beam[f])
    B_interp = interpolate.interp1d(
        ells[good_inds], mean_mars_beam[f][good_inds], 
        fill_value=(mean_mars_beam[f][good_inds][0], 0),
        bounds_error=False, kind='cubic')
    Bl0 = B_interp(all_ells)
    Bl_all_ells = np.vstack([Bl_all_ells, Bl0])
    headers.append(str(f))
    ax0.plot(Bl0, color=colors[f][0])
    ax_lowl.plot(Bl0[:1500], color='C4')
    ax_lowl_lo.plot(Bl0[:1500], color='C4')
    ax_lowl_hi.plot(Bl0[:1500], color='C4')

    if args.do_planck:
        # now get the various ways to get field map beam
        # (HMxS)/(HMxHM)*(Bl_p/Bl_s)=norm_fac (over some l range)
        # once have norm_fac, can then instead solve for Bl_s

        norm_bins = np.logical_and(ells > 600, 
                                   ells < 1000)
        norm_bins_mars = np.logical_and(ells > 600, 
                                        ells < 1000)

        sxs = field_xspecs['sxs'][f]
        sxp = np.mean([field_xspecs['sxhm1'][f], field_xspecs['sxhm2'][f]], 
                      axis=0)
        pxp = field_xspecs['hm1xhm2'][f]

        norm_fac_sxs = sxs[norm_bins] / \
                       (sxp[norm_bins] * mean_mars_beam[f][norm_bins_mars])
        norm_fac_sxs = np.mean(norm_fac_sxs)
        norm_fac_sxp = sxp[norm_bins] / \
                       (pxp[norm_bins] * mean_mars_beam[f][norm_bins_mars])
        norm_fac_sxp = np.mean(norm_fac_sxp)
        norm_fac_pxp = np.sqrt(sxs[norm_bins] / \
                               (pxp[norm_bins] * 
                                mean_mars_beam[f][norm_bins_mars]**2))
        norm_fac_pxp = np.mean(norm_fac_pxp)
        l1500_idxp = np.argmin(np.abs(ells-1500))
        l1500_idx = np.argmin(np.abs(ells-1500))
        ax_lowl.plot(ells[:l1500_idxp], 
                     (sxs / sxp / norm_fac_sxs)[:l1500_idxp], 
                     label=r'$\frac{\mathrm{SxS}}{\mathrm{SxP}}B_\ell^P, $'
                     r'$\epsilon_T=$'+'{:0.3f}'.format(norm_fac_sxs))
        ax_lowl.plot(
            ells[:l1500_idxp], 
            (np.sqrt(sxs / pxp) / norm_fac_pxp)[:l1500_idxp],
            label=r'$\sqrt{\frac{\mathrm{SxS}}{\mathrm{PxP}}}B_\ell^P$, ' 
            r'$\epsilon_T=$'+'{:0.3f}'.format(norm_fac_pxp))
        ax_lowl.plot(ells[:l1500_idxp], 
                     (sxp / pxp / norm_fac_sxp)[:l1500_idxp], 
                     label=r'$\frac{\mathrm{SxP}}{\mathrm{PxP}}B_\ell^P$, '
                     r'$\epsilon_T=$'+'{:0.3f}'.format(norm_fac_sxp))
        ax_lowl.errorbar(ells[:l1500_idx], mean_mars_beam[f][:l1500_idx], 
                         yerr=err_mars_beam[f][:l1500_idx], alpha=0.8, 
                         color='C4', linestyle='', label='Mars')
        ax_lowl.set_title(r'{}'.format(f))
        ax_lowl.set_ylabel(r'$B_\ell$')
        ax_lowl.set_xlabel(r'$\ell$')
        ax_lowl.legend()
        ax_lowl.set_ylim(0.85,1.1)
        fig_lowl.savefig(
            os.path.join(args.fig_dir, 
                         'planck_cross_{}.png'.format(f)))
        plt.close(fig_lowl)

        # do the lower half of the field
        sxs = field_xspecs_lo['sxs'][f]
        sxp = np.mean([field_xspecs_lo['sxhm1'][f], 
                       field_xspecs_lo['sxhm2'][f]], 
                      axis=0)
        pxp = field_xspecs_lo['hm1xhm2'][f]

        norm_fac_sxs = sxs[norm_bins] / \
                       (sxp[norm_bins] * mean_mars_beam[f][norm_bins_mars])
        norm_fac_sxs = np.mean(norm_fac_sxs)
        norm_fac_sxp = sxp[norm_bins] / \
                       (pxp[norm_bins] * mean_mars_beam[f][norm_bins_mars])
        norm_fac_sxp = np.mean(norm_fac_sxp)
        norm_fac_pxp = np.sqrt(sxs[norm_bins] / \
                               (pxp[norm_bins] * 
                                mean_mars_beam[f][norm_bins_mars]**2))
        norm_fac_pxp = np.mean(norm_fac_pxp)
        l1500_idxp = np.argmin(np.abs(ells-1500))
        l1500_idx = np.argmin(np.abs(ells-1500))
        ax_lowl_lo.plot(ells[:l1500_idxp], 
                        (sxs / sxp / norm_fac_sxs)[:l1500_idxp], 
                        label=r'$\frac{\mathrm{SxS}}{\mathrm{SxP}}B_\ell^P, $'
                        r'$\epsilon_T=$'+'{:0.3f}'.format(norm_fac_sxs))
        ax_lowl_lo.plot(
            ells[:l1500_idxp], 
            (np.sqrt(sxs / pxp) / norm_fac_pxp)[:l1500_idxp],
            label=r'$\sqrt{\frac{\mathrm{SxS}}{\mathrm{PxP}}}B_\ell^P$, ' 
            r'$\epsilon_T=$'+'{:0.3f}'.format(norm_fac_pxp))
        ax_lowl_lo.plot(ells[:l1500_idxp], 
                     (sxp / pxp / norm_fac_sxp)[:l1500_idxp], 
                     label=r'$\frac{\mathrm{SxP}}{\mathrm{PxP}}B_\ell^P$, '
                     r'$\epsilon_T=$'+'{:0.3f}'.format(norm_fac_sxp))
        ax_lowl_lo.errorbar(ells[:l1500_idx], mean_mars_beam[f][:l1500_idx], 
                         yerr=err_mars_beam[f][:l1500_idx], alpha=0.8, 
                         color='C4', linestyle='', label='Mars')
        ax_lowl_lo.set_title('{} Dec < -46 deg'.format(f))
        ax_lowl_lo.set_ylabel(r'$B_\ell$')
        ax_lowl_lo.set_xlabel(r'$\ell$')
        ax_lowl_lo.legend()
        ax_lowl_lo.set_ylim(0.85,1.1)
        fig_lowl_lo.savefig(
            os.path.join(args.fig_dir, 
                         'planck_cross_{}_low_dec.png'.format(f)))
        plt.close(fig_lowl_lo)

        # do the lower half of the field
        sxs = field_xspecs_hi['sxs'][f]
        sxp = np.mean([field_xspecs_hi['sxhm1'][f], 
                       field_xspecs_hi['sxhm2'][f]], 
                      axis=0)
        pxp = field_xspecs_hi['hm1xhm2'][f]

        norm_fac_sxs = sxs[norm_bins] / \
                       (sxp[norm_bins] * mean_mars_beam[f][norm_bins_mars])
        norm_fac_sxs = np.mean(norm_fac_sxs)
        norm_fac_sxp = sxp[norm_bins] / \
                       (pxp[norm_bins] * mean_mars_beam[f][norm_bins_mars])
        norm_fac_sxp = np.mean(norm_fac_sxp)
        norm_fac_pxp = np.sqrt(sxs[norm_bins] / \
                               (pxp[norm_bins] * 
                                mean_mars_beam[f][norm_bins_mars]**2))
        norm_fac_pxp = np.mean(norm_fac_pxp)
        l1500_idxp = np.argmin(np.abs(ells-1500))
        l1500_idx = np.argmin(np.abs(ells-1500))
        ax_lowl_hi.plot(ells[:l1500_idxp], 
                     (sxs / sxp / norm_fac_sxs)[:l1500_idxp], 
                     label=r'$\frac{\mathrm{SxS}}{\mathrm{SxP}}B_\ell^P, $'
                     r'$\epsilon_T=$'+'{:0.3f}'.format(norm_fac_sxs))
        ax_lowl_hi.plot(
            ells[:l1500_idxp], 
            (np.sqrt(sxs / pxp) / norm_fac_pxp)[:l1500_idxp],
            label=r'$\sqrt{\frac{\mathrm{SxS}}{\mathrm{PxP}}}B_\ell^P$, ' 
            r'$\epsilon_T=$'+'{:0.3f}'.format(norm_fac_pxp))
        ax_lowl_hi.plot(ells[:l1500_idxp], 
                     (sxp / pxp / norm_fac_sxp)[:l1500_idxp], 
                     label=r'$\frac{\mathrm{SxP}}{\mathrm{PxP}}B_\ell^P$, '
                     r'$\epsilon_T=$'+'{:0.3f}'.format(norm_fac_sxp))
        ax_lowl_hi.errorbar(ells[:l1500_idx], mean_mars_beam[f][:l1500_idx], 
                         yerr=err_mars_beam[f][:l1500_idx], alpha=0.8, 
                         color='C4', linestyle='', label='Mars')
        ax_lowl_hi.set_title('{} Dec > -46 deg'.format(f))
        ax_lowl_hi.set_ylabel(r'$B_\ell$')
        ax_lowl_hi.set_xlabel(r'$\ell$')
        ax_lowl_hi.legend()
        ax_lowl_hi.set_ylim(0.85,1.1)
        fig_lowl_hi.savefig(
            os.path.join(args.fig_dir, 
                         'planck_cross_{}_high_dec.png'.format(f)))
        plt.close(fig_lowl_hi)

    popt, pcov = ba.fit_gauss_lspace(ells[good_inds], 
                                     mean_mars_beam[f][good_inds])
    print('{} final fwhm: {} arcmin'.format(f, np.abs(popt[0]*60)))

    ax0.errorbar(ells, mean_mars_beam[f], linestyle='',
                 yerr=err_mars_beam[f], color=colors[f][0],
                 label='{} Mars fwhm={:.1f} arcmin'.format(f, 
                                                      np.abs(popt[0]*60)))

ax0.set_xlabel(r'$\ell$')
ax_highl.set_xlabel(r'$\ell$')
ax0.set_ylabel(r'$B_\ell$')
ax_highl.set_ylabel(r'$B_\ell$')
ax0.legend()
fig0.savefig(os.path.join(args.fig_dir, 'Bl_allfreq.png'))
ax_highl.legend()
fig_highl.savefig(os.path.join(args.fig_dir, 'Bl_allfreq_highl.png'))

np.savetxt(os.path.join(args.output_path, 'beam00.txt'), Bl_all_ells.T, 
           header=' '.join(headers))
plt.show()
