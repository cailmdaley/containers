'''
Script to make planet maps in source-centered coordinates, either with
data or a simulated input map. Defaults are set for beam determination,
including flagging scans that saturate on the planet.

Because filtering can cause striping features if the source is not masked,
a point source mask should be used. This mask is constructed by finding the 
brightest point in the map specified by the mask_map argument. Because
the planet tends to have a constant offset from the center of the map 
independent of observation, one can make a single first order poly-filtered
map for the mask_map, and feed that map to all ensuing higher order
poly-filtered maps.

If running locally, can set input_files argument to the observation
directory containing TODs and calibration files. If submitting to the
grid, use submit_planet_maps.py, which points to this script.
'''

#Note, this script has been disassembled and reassembled by Neil
#We want to implement a number of things, first of which is a larger flag radius
#We also want to make it run with the current version of the 3g pipeline 

import argparse as ap
import numpy as np
import os
from glob import glob
from spt3g.coordinateutils.coordsysmodules import FillCoordTransRotations
from spt3g import core, std_processing, mapmaker, calibration, \
    coordinateutils, timestreamflagging
from spt3g.pointing import CalculateCoordTransRotations, \
    CalculateLocalOffsetPointing

#Note from past Neil, this should be killed once we are done debugging
#HACK, to be removed before this script is committed
from matplotlib.pyplot import *

P = ap.ArgumentParser(description='Make planet maps with special options',
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input_files', action='store', nargs='+', default=[],
               help='Input files')
P.add_argument('--cal-file', action='store', default=None,
               help='If input_files is a directory and the calibration file'
               ' in that directory is not wanted, set cal_file to desired'
               ' calibration g3 file.')
P.add_argument('-o', '--output', action='store', default='output.g3',
               help='Output filename')
P.add_argument('-s', '--source', action='store', 
               default='mars', help='name of source')
P.add_argument('--poly-order', action='store', type=int,
               default=3, help='polynomial filter order')
P.add_argument('-k', '--source-relative', action='store_true',
               default=False, help='calibrate in source-relative units rather '
               'than K_cmb')
P.add_argument('-r', '--res', action='store', 
               default=0.2, help='resolution [arcmin]')
P.add_argument('-x', '--xlen', action='store', type=float,
               default=4, help='map width [deg]')
P.add_argument('-y', '--ylen', action='store', type=float,
               default=4, help='map height [deg]')
P.add_argument('--bands', action='store', nargs='+', 
               default=['90', '150', '220'], 
               help='Frequency bands to make maps for')
P.add_argument('--wafer-split', action='store_true', default=False,
               help='Split maps by wafer')
P.add_argument('--lr', action='store_true', default=False,
               help='Split maps by scan direction')
P.add_argument('--mask-map', action='store', default=None,
               help='Map to use to get point source mask')
P.add_argument('--no-flag-saturated', dest='flag_saturated', 
               action='store_false', default=True,
               help='Do not flag detectors that saturate on Mars.')
P.add_argument('--produce-simstub', action='store_true', default=False,
               help='Produce a simstub for mock observations.')
P.add_argument('--sim', action='store_true', default=False,
               help='Mock observe a 1 arcmin gaussian point source. Input '
               'files should be simstubs.')
P.add_argument('--mask-radius', type=float, default=30.,
                help='Radius in arcmin around the center pixel to mask, default is 30')
P.add_argument('--saturate-threshold', type=float, default=1.15,
                help='Flagging threshold for Annes figure of merit on Mars saturation, voided by --no-flag-saturated')

args = P.parse_args()

res = float(args.res) * core.G3Units.arcmin
xlen_pix = int(args.xlen * core.G3Units.degrees / res)
ylen_pix = int(args.ylen * core.G3Units.degrees / res)

bands = []
for band in args.bands:
    bands.append(band+'GHz')

#Note, this looks as expected...
#print(bands)

# If input_files is a directory, search it for data
# Use cal file in that directory if cal_file arg isn't specified.
if os.path.isdir(args.input_files[0]):
    if args.cal_file is None:
        cal = glob(os.path.join(args.input_files[0],'*offline_calibration.g3'))
    else:
        cal = [args.cal_file]
    dat = sorted(glob(os.path.join(args.input_files[0],'[0-9]*.g3')))
    # calibration must be first argument!
    args.input_files = cal + dat

# If splitting by wafer, get a list of wafers
if args.wafer_split:
    wafer_list = []
    for fname in args.input_files:
        for frame in core.G3File(fname):
            if 'BolometerProperties' in frame:
                for bolo, props in frame['BolometerProperties'].items():
                    if props.wafer_id not in [None,'']:
                        waf = props.wafer_id.capitalize()
                        if waf not in wafer_list:
                            wafer_list.append(waf)
                break
        if len(wafer_list) > 0:
            break

# Generate point source mask
# this is used for 1) fitting the polynomial coefficients in the filtering, 
# and 2) computing the inverse variance of the filtered timestream for 
# weighting
if args.mask_map is not None:
    # Find mask_map's brightest point to mask it
    mask_map = core.G3File(os.path.join(args.mask_map))
    print('got file')
    xvals = []
    yvals = []
    xshapes=[]
    yshapes=[]
    for fr in mask_map:
        if fr.type == core.G3FrameType.Map:
            tmap = np.asarray(fr['T'])
            max_x, max_y = np.unravel_index(np.abs(tmap).argmax(), tmap.shape)
            xvals.append(max_x)
            yvals.append(max_y)
            xshapes.append(tmap.shape[0])
            yshapes.append(tmap.shape[1])
    
    #This is just using argmax and assuming that the brightest thing is Mars
    #empirically this is fine
    #For some reason this breaks with the sparse_maps
    #std::vector<double> FlatSkyMap::pixel_to_angle(size_t x, size_t y) const {
    #    return proj_info.pixel_to_angle(y*xpix_ + x, false);
    #}

    # call the best location the median of the three frequency maps' centers
    xloc = int(np.median(xvals))
    yloc = int(np.median(yvals))
    xlen = int(np.median(xshapes))
    ylen = int(np.median(yshapes))
    print('Masking pixel [{}, {}]'.format(xloc,yloc))
    # numpy array indexing-- need to switch x and y here
    mask_center = fr['T'].pixel_to_angle(yloc * xlen + xloc)
    radius = args.mask_radius *core.G3Units.arcmin

    psmask = coordinateutils.FlatSkyMap(
        x_len=xlen_pix, y_len=ylen_pix, alpha_center=0, delta_center=0, res=res,
        proj=coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
        pol_type=coordinateutils.MapPolType.T, coord_ref=coordinateutils.MapCoordReference.Local)
    ra, dec = list(mask_center) 

    #Replacement line of code
    mapmaker.make_point_source_mask_cpp([ra,], [dec,], [radius,], False, psmask)

    #Anne's original method to do point source masking
    #mapmaker.fill_point_source_mask_flat_map([ra], [dec], [radius], False, 
    #                                         False, psmask)
    #Note to future Neil, this seems to be doing the correct thing now
    
    ps_map_id = 'PointSourceMap'
else:
    ps_map_id = None

def hack_y_offset(frame):
    """
    Flip sign of dector y-offset from boresight. For some reason, this
    is needed for source-centered maps.
    """
    if 'BolometerProperties' in frame.keys():
        for det in frame['BolometerProperties'].keys():
            frame['BolometerProperties'][det].y_offset *= -1

def SimStub(fr, ts_key, valid_ids_key='valid_ids', flag_key='Flags', 
            weight_key='TodWeights'):
    """
    Create the file containing the keys that sims require listed in to_keep
    """
    to_keep = [valid_ids_key, flag_key, weight_key,
               'DfMuxHousekeeping', 'RawBoresightAz', 'RawBoresightEl',
               'OnlineBoresightAz', 'OnlineBoresightEl', 'OnlineBoresightRa',
               'OnlineBoresightDec', 'OnlineRaDecRotation',
               'OnlinePointingModel', 'OffsetRotation']
   
    if fr.type == core.G3FrameType.Scan:
        assert(ts_key in fr)
        fr[valid_ids_key] = core.G3VectorString(fr[ts_key].keys())
        for k in fr:
            if k not in to_keep:
                del fr[k]

def check_for_sim_keys(fr, valid_ids_key = 'valid_ids',
                       out_tsm_key = 'valid_ids_tsm'):
    '''
    CalculatePointing expects a TimestreamMap, so make sure there's one
    in the frame when doing mock-observations.
    valid_ids_key should already exist in simstub and point to a list.
    '''
    if fr.type == core.G3FrameType.Scan:
        if not valid_ids_key in fr:
            raise KeyError(
                "%s not found in frame. "%valid_ids_key +
                "List of valid bolometer ids required for mock-observing.")
        tsm = core.G3TimestreamMap()
        for bid in fr[valid_ids_key]:
            tsm[bid]=core.G3Timestream()
        fr[out_tsm_key] = tsm

# Generate map stub
smstub = coordinateutils.FlatSkyMap(
    x_len=xlen_pix, y_len=ylen_pix, alpha_center=0, delta_center=0, res=res,
    proj=coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
    pol_type=coordinateutils.MapPolType.T, coord_ref=coordinateutils.MapCoordReference.Local)

# Begin pipeline
pipe = core.G3Pipeline()

pipe.Add(core.G3Reader, filename=args.input_files)
    
if not args.sim:
    pipe.Add(hack_y_offset)

    pipe.Add(core.DeduplicateMetadata)

    # Cut turnarounds
    pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])

    # Clean up pre-existing timestreams
    pipe.Add(
        core.Delete,
        keys=['OnlineBoresightAz', 'OnlineBoresightEl',
              'OnlineBoresightRa', 'OnlineBoresightDec', 'OnlineRaDecRotation']
    )

    pipe.Add(
        CalculateCoordTransRotations,
        raw_az_key='RawBoresightAz',
        raw_el_key='RawBoresightEl',
        output='OnlineBoresight',
        transform_store_key='OnlineRaDecRotation',
        model='OnlinePointingModel',
        flags=['az_tilts', 'el_tilts', 'flexure', 'collimation', 'refraction']
    )

    pipe.Add(CalculateLocalOffsetPointing, source=args.source,
             x_offset_key='AzOffset', y_offset_key='ElOffset',
             trans_key='OnlineRaDecRotation', ts_ref_key='OnlineBoresightAz',
             max_throw=args.xlen * core.G3Units.deg * 1.66)

    pipe.Add(FillCoordTransRotations,
             transform_store_key='OffsetRotation',
             bs_az_key='RawBoresightAz', bs_el_key='RawBoresightEl',
             bs_ra_key='AzOffset', bs_dec_key='ElOffset',
             do_bad_transform=True)

    # do some pre-calibration flagging
    
    #Looks for NaN in the data and negative DAN gain channels, among other things. Probably no need to touch this
    pipe.Add(std_processing.flagsegments.FlagInvalidData, flag_key='Flags', 
             ts_key='RawTimestreams_I')
    #Checks the chopped thermal source and elnod signal to noise and cuts if these are below some threshold
    pipe.Add(std_processing.flagsegments.FlagNonResponsive, flag_key='Flags')



    pipe.Add(std_processing.flagsegments.FlagUncalibratable, 
             ts_key='RawTimestreams_I', flag_key='Flags')
    if args.flag_saturated:
        pipe.Add(timestreamflagging.RemoveFlagged, 
                 input_ts_key='RawTimestreams_I', input_flag_key='Flags', 
                 output_ts_key='DeflaggedTimestreamsFirst')
        #Note from Neil, this eventually goes to timestreamflagging/python/miscflagmodules.py:FlagSaturatedBolos
        #TODO run this code sweeping ver the threshold parameter


        pipe.Add(std_processing.flagsegments.FlagSaturatedBolosMars, 
                 ts_key='RawTimestreams_I', flag_key='Flags',thresh=args.saturate_threshold)
        
        #Now we need to check the focal plane distribution of what is being cure from these maps!!
        #Will have to cache this data somehow

        pipe.Add(timestreamflagging.RemoveFlagged, 
                 input_ts_key='DeflaggedTimestreamsFirst',
                 input_flag_key='Flags', output_ts_key='DeflaggedTimestreams')
        pipe.Add(core.Delete, keys=['DeflaggedTimestreamsFirst'])

    else:
        pipe.Add(timestreamflagging.RemoveFlagged, 
                 input_ts_key='RawTimestreams_I', input_flag_key='Flags', 
                 output_ts_key='DeflaggedTimestreams')

    # Collect statistics on flagged detectors
    stats = timestreamflagging.GenerateFlagStats(flag_key='Flags')
    pipe.Add(stats)

    # Calibrate, the timestreams
    pipe.Add(std_processing.CalibrateRawTimestreams, 
             units=core.G3TimestreamUnits.Power, output='TimestreamsWatts')
    
    #It would be interesting do dig at this a bit and see what the map looks like in units of watts
    pipe.Add(calibration.ApplyTCalibration, InKCMB=not args.source_relative,
             OpacityCorrection=not args.source_relative, 
             Input='TimestreamsWatts', Output='DeflaggedCalTimestreams')

else:
    # Do unique to sim stuff
    if args.mask_map is None:
        raise ValueError('For sims, must specify a mask_map argument to '
                         'determine center pixel of sim map')

    # Make sure a list of valid ids exists and create dummy timestream map
    pipe.Add(check_for_sim_keys, valid_ids_key='valid_ids',
             out_tsm_key='valid_ids_tsm')

    # make sim map to reobserve
    res_s = 0.1 * core.G3Units.arcmin
    xlen_s = int(7 * core.G3Units.deg / res_s)
    ylen_s = int(7 * core.G3Units.deg / res_s)
    if xlen_s % 2 == 0:
        xlen_s += 1
        ylen_s += 1

    mars_amp = 5000.
    # make a 2-d 1 arcmin gaussian 
    # first, get distance from center in arcmin for where the source is in mask
    # map 
    xcen_am = (xloc - xlen_pix / 2) * res
    ycen_am = (yloc - ylen_pix / 2) * res
    # now put it in source pixels
    xcen_s = int((xcen_am + xlen_s*res_s / 2) / res_s)
    ycen_s = int((ycen_am + ylen_s*res_s / 2) / res_s)
    width_s = 1. * core.G3Units.arcmin / (2*np.sqrt(2*np.log(2))) / res_s
    xp_s, yp_s = np.indices((xlen_s, ylen_s))
    # build the 2-d gaussian array
    g = mars_amp * np.exp(-(((xcen_s - xp_s)/width_s)**2. + 
                            ((ycen_s - yp_s)/width_s)**2.)/2.)

    #Note to future Neil, I don't think this is quite right, actually
    #We probably want to put in a Gaussian + theta^-3 tail to better approximate the large scale structure of the beam profile
    #Then we can simulate running that through the TOD scanning pipeline and compute the biased B_ell
    #Alternately we can work in the no-filtering regime and compute cross spectra as a consistency check \
    #that we can debias, but I should think about that later


    m_sim = coordinateutils.FlatSkyMap(
        x_len=xlen_s, y_len=ylen_s, res=res_s, alpha_center=0, delta_center=0,
        proj=coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea, 
        coord_ref=coordinateutils.MapCoordReference.Local, pol_type=coordinateutils.MapPolType.T)

    np.asarray(m_sim)[:] = g
    pipe.Add(mapmaker.MapInjector, map_id='sim_map',
             maps_lst=[m_sim], is_stub=False)

    pipe.Add(mapmaker.CalculatePointing, map_id='sim_map',
        pointing_store_key='SimPixelPointing',
        ts_map_key='valid_ids_tsm',
        trans_key='OffsetRotation')

    # Generate simulated TOD
    pipe.Add(mapmaker.mapmakerutils.FillSimTodSegment,
             out_ts_key='DeflaggedCalTimestreams',
             ts_to_get_sample_rate_key='OnlineBoresightRa',
             interp_sim=False,
             map_is_healpix=False,
             trans_key='OffsetRotation',
             sim_map_id='sim_map',
             valid_ids='valid_ids',
            sim_pointing_key='SimPixelPointing')

if args.mask_map is not None:
    pipe.Add(mapmaker.MapInjector, map_id=ps_map_id, maps_lst=[psmask], 
             is_stub=False)

pipe.Add(mapmaker.MapInjector, map_id='bsmap', maps_lst=[smstub,], 
         is_stub=False)
pipe.Add(mapmaker.CalculatePointing, map_id='bsmap',
         pointing_store_key='PixelPointing', 
         ts_map_key='DeflaggedCalTimestreams', 
         trans_key='OffsetRotation')

# Basic timestream filtering
#This appears to go to 

pipe.Add(mapmaker.TodFiltering, ts_in_key='DeflaggedCalTimestreams',
         ts_out_key='PolyFilteredTimestreams', use_dynamic_source_filter=False, 
         poly_order=args.poly_order, point_source_mask_id=ps_map_id, 
         point_source_pointing_store_key='PixelPointing', 
         filter_mask_key='FilterMask')

if not args.sim:
    # Standard bolometer weighting
    if args.mask_map is not None:
        # Mask the point source when calculating weight
        pipe.Add(std_processing.weighting.AddMaskedVarWeight, 
                 input='PolyFilteredTimestreams', output='TodWeights')
    else:
        pipe.Add(std_processing.weighting.AddSigmaClippedWeight,
                 input='PolyFilteredTimestreams', output='TodWeights')

    #TODO note to future Neil, figure out why this is killing all of the data...
    #pipe.Add(timestreamflagging.noiseflagging.FlagUnphysicallyLowVariance, 
    #         ts_key='PolyFilteredTimestreams') 

pipe.Add(timestreamflagging.RemoveFlagged,  
         input_ts_key='PolyFilteredTimestreams', input_flag_key='Flags',
         output_ts_key='DeflaggedPFTimestreams')

pipe.Add(lambda frame: print(frame))

if args.sim:
    stats = timestreamflagging.GenerateFlagStats(flag_key='Flags')
    pipe.Add(stats)
else:
    pipe.Add(stats.update)

pipe.Add(std_processing.AddMetaData, stats=stats, args=vars(args))

# Clean up detritus
pipe.Add(core.Delete, keys=['RawTimestreams_I', 'RawTimestreams_Q', 
                            'TimestreamsWatts', 'CalTimestreams',
                            'PolyFilteredTimestreams',
                            'DeflaggedTimestreams'])
pipe.Add(core.Dump)

# split left/right scans (optional)
if args.lr:
    pipe.Add(std_processing.ScanDirections.split_left_right_scans,
             ts_in_key='DeflaggedPFTimestreams',
             boresight_az_key='OnlineBoresightAz')
    pipe.Add(core.Delete, keys=['DeflaggedPFTimestreams'])
    directions = ['Left', 'Right']
else:
    directions = ['']

# split by wafer (optional)
if args.wafer_split:
    wafers = wafer_list
    for direction in directions:
        pipe.Add(calibration.SplitByWafer, 
                 input='DeflaggedPFTimestreams'+direction,
                 output_root='DeflaggedPFTimestreams'+direction)
    pipe.Add(core.Delete, keys=['DeflaggedPFTimestreams'+direction])
else:
    wafers = ['']


# split by band (always)
for direction in directions:
    for wafer in wafers:
        print('direction/wafer')
        print(direction)
        print(wafer)
        pipe.Add(calibration.SplitByBand,
                 input='DeflaggedPFTimestreams'+direction+wafer,
                 output_root='DeflaggedPFTimestreams'+direction+wafer)
        pipe.Add(core.Delete, keys=['DeflaggedPFTimestreams'+direction+wafer])
mapids = ['%s%s%s' % (direction, wafer, band)
          for direction in directions for wafer in wafers for band in bands]

for mapid in mapids:
    pipe.Add(mapmaker.MapInjector, map_id=mapid, maps_lst=[smstub], 
             is_stub=True, make_polarized=False, do_weight=True)
    pipe.Add(mapmaker.mapmakerutils.BinMap, map_id=mapid,
             ts_map_key='DeflaggedPFTimestreams%s' % (mapid),
             trans_key='OffsetRotation', pointing_store_key='PixelPointing',
             timestream_weight_key='TodWeights')

if args.produce_simstub:
    # Write sim stub
    pipe.Add(SimStub, ts_key='DeflaggedCalTimestreams')
    pipe.Add(core.G3Writer,
             filename = args.output.replace(
                 args.output.split('/')[-1],
                 'simstub_'+args.output.split('/')[-1]),
             streams=[core.G3FrameType.Observation, core.G3FrameType.Wiring,
                      core.G3FrameType.Calibration, core.G3FrameType.Scan,
                      core.G3FrameType.PipelineInfo,
                      core.G3FrameType.EndProcessing])

# Final cleanup
pipe.Add(lambda fr: fr.type != core.G3FrameType.Scan)

# Write map to disk
pipe.Add(core.G3Writer, filename=args.output)

pipe.Run()


