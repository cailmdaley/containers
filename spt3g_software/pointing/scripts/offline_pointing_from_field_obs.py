#!/usr/bin/env python

#    Wrapper to pointing model code intended to be used
#    for calculating offline pointing model corrections
#    based on source locations in CMB field maps (or
#    source-centered cutout maps from CMB field observation
#    data).

from __future__ import print_function

import numpy as np
from spt3g import core, std_processing, pointing, mapmaker, coordinateutils
from spt3g.mapmaker import mapmakerutils as mm
from spt3g.pointing import pointing_model_fit_tools as pmf
from spt3g.pointing.astrometry import check_astrometry_at20
from spt3g.pointing import offline_pointing as op
from spt3g.sources.source_utils import parse_point_source_file
import argparse as ap
import sys, glob
from scipy import ndimage
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord

# Usage: offline_pointing_from_field_obs.py <files.g3> [-o output.g3]

P = ap.ArgumentParser(
    description='Get pointing model corrections from source locations in CMB field observations',
    formatter_class=ap.ArgumentDefaultsHelpFormatter,
)
P.add_argument('input', action='store', nargs='+', default=[], help='Input files')
P.add_argument(
    '-o',
    '--output',
    action='store',
    default=None,
    help='Output filename, if none provided will use obsid',
)
P.add_argument(
    '-n',
    '--n-sigma',
    action='store',
    default='10',
    help='Signal-to-noise threshold for including sources in fit',
)
P.add_argument(
    '--fixed-params',
    nargs='+',
    choices=[
        'boom_terms',
        'az_tilts',
        'el_tilt',
        'x_collimation',
        'y_collimation',
        'az_encoder_offset',
    ],
    default=['boom_terms', 'az_tilts', 'x_collimation', 'az_encoder_offset'],
)
P.add_argument(
    '-f',
    '--maps-used-offline-pointing',
    action='store_true',
    default='True',
    help='Input maps were made using offline pointing corrections',
)
P.add_argument(
    '-s',
    '--sources-to-ignore-file',
    action='store',
    default='',
    help='File with sources to ignore in fit (for instance if we want to save these for measuring the beam)',
)
args = P.parse_args()

# parse some args
nsigma = np.float(args.n_sigma)
if args.sources_to_ignore_file != '':
    indices, ras2ignore, decs2ignore, source_radii = parse_point_source_file(
        args.sources_to_ignore_file
    )
    coords2ignore = SkyCoord(
        ras2ignore / core.G3Units.deg, decs2ignore / core.G3Units.deg, unit="deg"
    )
else:
    coords2ignore = SkyCoord(0.0, 90.0, unit="deg")

# Parse input files, grab pointing model parameters used and maps,
# hand maps to astrometry checker.
input_files = args.input
obsid = None
# assume if there is only one file, it's a full-field map; otherwise
# assume there are many source-centered thumbnail maps. No, actually,
# the thumbnail maps are many frames in a single file. Anyway, also
# grab calframe (for offline pointing corrections already applied) and
# the online pointing model (if it's there).
maplist = []
f1 = core.G3File(input_files)
scanframe = None
calframe = None
for thisframe in f1:
    if thisframe.type is core.G3FrameType.Observation:
        obsid = thisframe['ObservationID']
    if thisframe.type is core.G3FrameType.Calibration:
        calframe = thisframe
    if thisframe.type is core.G3FrameType.Map:
        if 'GHz' in thisframe['Id']:
            frame = thisframe
            thiscoord = SkyCoord(
                frame['T'].alpha_center / core.G3Units.deg,
                frame['T'].delta_center / core.G3Units.deg,
                unit="deg",
            )
            # only include in maplist if source is not within 10 arcmin
            # of a source in the "ignore" list
            if np.min(coords2ignore.separation(thiscoord).arcmin > 10.0):
                if 'OnlinePointingModel' in frame.keys():
                    model_from_frame = frame['OnlinePointingModel']
                mm.RemoveWeightModule(frame)
                maplist.append(frame)

if len(maplist) == 1:
    # for full-field map, get full map, smooth, make pixel mask, hand to
    # AT20G comparing code. also grab offline pointing corrections (if
    # any were used).
    fwhm_pix = 1.8 / (frame['T'].res / core.G3Units.arcmin)
    map_sm = ndimage.gaussian_filter(np.asarray(frame['T']), fwhm_pix * 0.42)
    map_new = coordinateutils.FlatSkyMap(
        map_sm,
        res=frame['T'].res,
        is_weighted=False,
        alpha_center=frame['T'].alpha_center,
        delta_center=frame['T'].delta_center,
        proj=frame['T'].proj,
        pol_type=core.MapPolType.T,
    )
    frame_new = core.G3Frame(core.G3FrameType.Map)
    frame_new['Wunpol'] = frame['Wunpol']
    frame_new['T'] = map_new
    wtemp = np.asarray(frame_new['Wunpol'].TT)
    wtemp_sm = ndimage.gaussian_filter(wtemp, 16.0)
    medwt = np.median(wtemp_sm[np.where(wtemp_sm > np.max(wtemp_sm) / 2.0)])
    whok = np.where(wtemp_sm > 0.8 * medwt)
    pixel_mask = np.zeros(wtemp_sm.shape)
    pixel_mask[whok] = 1
    pdict = check_astrometry_at20(
        frame_new,
        check_beam_size=False,
        close_cut_arcsec=360.0,
        pixel_mask=pixel_mask,
        nsigma=nsigma,
        plot=False,
    )
    # also try and find the raw data from this obs to get
    # online pointing model
    data_dirs = [
        '/spt_data/bolodata/downsampled/',
        '/spt/data/bolodata/downsampled/',
    ]  # Pole vs. Chicago
    for dd in data_dirs:
        scanfile = glob.glob(dd + '/*/' + str(obsid) + '/0000.g3')
        if len(scanfile) > 0:
            f2 = core.G3File(scanfile[0])
            for thisframe in f2:
                if thisframe.type is core.G3FrameType.Scan:
                    scanframe = thisframe
                    break
else:
    # concantenate all thumbnail map frames, hand them to astrometry
    # checker. also get obsid and online pointing (and any offline
    # corrections) from one of the frames.
    if scanframe is None:
        scanframe = core.G3Frame(core.G3FrameType.Scan)
        scanframe['OnlinePointingModel'] = model_from_frame
        pdict = check_astrometry_at20(
            maplist,
            check_beam_size=False,
            close_cut_arcsec=360.0,
            nsigma=nsigma,
            plot=False,
            gauss_fit_switch=True,
        )

# get measured offsets and pointing model used into format
# expected by pointing model fitter. assume el = -dec,
# az = ra - lst (need to calculate lst), delta_el = -delta_dec,
# and delta_az = delta_ra. both check_astrometry_at20 and
# pointing_model_fit_tools are terrible about units, so you just
# have to know that both work in degrees.
model_nocorr = pmf.read_model(scanframe)
if args.maps_used_offline_pointing == 'True':
    corrector = op.ApplyPointingCorrection()
    corrector(calframe)
    corrector(scanframe)
    model = pmf.read_model(scanframe, use_offline_model=True)
else:
    model = model_nocorr
az_offset = pdict['dra']
el_offset = -pdict['ddec']
t1 = Time(
    std_processing.obsid_to_g3time(obsid).mjd,
    format='mjd',
    location=('-44.65d', '-89.999d'),
)
dt1 = TimeDelta(3600.0, format='sec')
t2 = t1 + dt1
mean_lst = t2.sidereal_time('apparent').value * 15.0
az = np.mod(pdict['ra'] - mean_lst, 360.0)
el = -pdict['dec']

# make comma-separated file in the format pmf expects
datafile = str(obsid) + '_source_offset_data.txt'
f = open(datafile, 'w')
f.write('#Written by offline_pointing_from_field_obs.py. \n')
f.write(
    '#tilts '
    + str(model['a2'])
    + ', '
    + str(model['a3'])
    + ', '
    + str(model['a4'])
    + ' \n'
)
f.write('#flexure radio, ' + str(model['a0']) + ', ' + str(model['a1']) + ' \n')
f.write('#collimate_fixed radio, ' + str(model['a5']) + ', ' + str(model['a6']) + ' \n')
f.write(
    '# Name, R.A. [deg], dec. [deg], MJD, Source az [deg], Source el [deg], Measured Az [deg], dAz [arcmin], Measured El [deg], dEl [arcmin], tilt0 [deg], tilt1 [deg], tilt2 [deg] \n'
)
for q in np.arange(len(az)):
    f.write(
        'source'
        + str(q)
        + ', '
        + str(pdict['ra'][q])
        + ', '
        + str(pdict['dec'][q])
        + ', '
        + str(std_processing.obsid_to_g3time(obsid).mjd)
        + ', '
        + str(az[q])
        + ', '
        + str(el[q])
        + ', '
        + str(az[q] + az_offset[q])
        + ', '
        + '0.001, '
        + str(el[q] + el_offset[q])
        + ', '
        + '0.001, '
        + str(model['a2'])
        + ', '
        + str(model['a3'])
        + ', '
        + str(model['a4'])
        + ' \n'
    )
f.close()

# point fitter at file
fit_result = pmf.run(
    str(obsid), spt_files=[datafile], fixed_params=args.fixed_params, process=True
)

# gather output
fit_params = fit_result[1]['fit_params']
dpmodel = core.G3MapVectorDouble()
dpmodel['fixedCollimation'] = (
    np.asarray([fit_params['a5_spt'], fit_params['a6_spt']])
    - np.asarray([model_nocorr['a5'], model_nocorr['a6']])
) * core.G3Units.deg
dpmodel['flexure'] = (
    np.asarray([fit_params['a0'], fit_params['a1']])
    - np.asarray([model_nocorr['a0'], model_nocorr['a1']])
) * core.G3Units.deg
dpmodel['tilts'] = (
    np.asarray([fit_params['a2'], fit_params['a3'], fit_params['a4']])
    - np.asarray([model_nocorr['a2'], model_nocorr['a3'], model_nocorr['a4']])
) * core.G3Units.deg

# write new offline pointing model corrections to file
calframe['FieldScanPointingModelCorrection'] = dpmodel

if args.output is None:
    output_file = str(obsid) + '.g3'
else:
    output_file = args.output
writer = core.G3Writer(output_file)
writer(calframe)
writer(core.G3Frame(core.G3FrameType.EndProcessing))  # Make sure data gets flushed
