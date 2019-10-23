from __future__ import print_function
from spt3g import core, calibration, coordinateutils
from spt3g.mapmaker import MapBinner, MapPointingCalculator
from spt3g.mapmaker import SimulatedTimestreamFiller, remove_weight
from spt3g.mapmaker import remove_weight_t_cpp, make_determinant_map_cpp
from spt3g.mapmaker import MapPointingCalculator, MapBinner
from spt3g.mapmaker import SkyMapInterpTodFiller
from spt3g.mapmaker import SimulatedTimestreamFiller

from spt3g.coordinateutils import ang_to_quat, quat_to_ang
from spt3g import todfilter

import numpy as np
import copy


'''
To understand how things work it's really important to understand the format of the Map frames.

This is documented in the frames.rst file, but just as a quick bit of information, each
Map frame has a field called Id that stores a string identifying which map this is.

This allows us to pass point source masks, individual bolometer maps and polarized coadds around
and have a variety of modules work on them.  In the case of instances where we are storing information in a map frame that doesn't map to a physical signal (like the point source mask) store the map in the "T" field of the frame.  In this case, T stands for "The map."

'''

def find_timestreams_with_pointing(ra, dec, ts_ra_map, ts_dec_map, slop):
    '''
    A utility function to find list of timestreams in ts_ra/dec_map within distance, slop, of (ra,dec).
    '''
    out_ks = []
    for k in ts_ra_map.keys():
        ras = np.asarray(ts_ra_map[k]) % ( 2.0 * np.pi )
        decs = np.asarray(ts_dec_map[k])

        ra_delts = np.asarray(np.abs(ras - ra))
        dec_delts = np.asarray(np.abs(decs - dec))
        print(np.min(ra_delts), np.min(dec_delts))
        #print ra_delts, ras, ra
        good_inds = np.where( (ra_delts <= slop) & (dec_delts <= slop) )[0]

        if len(good_inds) > 0:
            out_ks.append(k)
    return out_ks


@core.indexmod
class CheckForCalFrames(object):
    '''
    Raises an exception if we get scan frames before cal frames
    '''

    def __init__(self):
        self.has_cal_frame = False
        self.bolo_props_keys = None
    def __call__(self, frame):
        if frame.type == core.G3FrameType.Scan:
            if not self.has_cal_frame:
                core.log_fatal("Calibration frame not present before scan frames.")
        if frame.type == core.G3FrameType.Calibration:
            self.has_cal_frame = True

@core.indexmod
def ValidateMapFrames(frame):
    '''
    Raises a RuntimeError if any map frames do not fit the proper structure.
    '''
    if (not isinstance(frame, core.G3Frame)) or (not isinstance(frame, dict)):
        return
    if isinstance(frame, core.G3Frame):
        if frame.type != core.G3FrameType.Map:
            return
    if 'T' not in frame:
        raise RuntimeError("Right now the code assumes T will always be in a map frame.")
    if 'Q' in frame and not 'U' in frame:
        raise RuntimeError("Q is in map frame, but U is not")
    if 'U' in frame and not 'Q' in frame:
        raise RuntimeError("U is in map frame, but Q is not")
    if 'Wpol' in frame and 'Wunpol' in frame:
        raise RuntimeError("Both kinds of weight are present in the map")

def BoresightPointingOffsetCalc(frame,
                                boresight_alpha_key,
                                boresight_delta_key,
                                out_detector_alpha_key,
                                out_detector_delta_key):
    '''
    Used to change the type of frame[boresight_delta/alpha_key] from G3MapDouble to
    G3MapVectorDouble where the data is stored in the key "Boresight"
    '''

    if frame.type != core.G3FrameType.Scan:
        return
    out_alpha = core.G3MapVectorDouble()
    out_delta = core.G3MapVectorDouble()
    out_alpha['Boresight'] = core.G3VectorDouble(frame[boresight_alpha_key])
    out_delta['Boresight'] = core.G3VectorDouble(frame[boresight_delta_key])
    frame[out_detector_alpha_key] = out_alpha
    frame[out_detector_delta_key] = out_delta

@core.indexmod
def AddBoresightBoloProps(frame, new_bp_key="BoresightBoloProps"):
    if frame.type != core.G3FrameType.Calibration:
        return
    bpm = calibration.BolometerPropertiesMap()
    bp = calibration.BolometerProperties()
    bp.x_offset = 0
    bp.y_offset = 0
    bp.pol_efficiency = 0
    bp.pol_angle = 0
    bpm["Boresight"] = bp
    frame[new_bp_key]=bpm

def AddBsTimestream(frame, tsm_key = 'BoresightTimestreams'):
    if frame.type != core.G3FrameType.Scan:
        return
    tsm = core.G3TimestreamMap()
    tsm['Boresight'] = core.G3Timestream()
    frame[tsm_key] = tsm

def ExpandBoresightPointing(frame,
                            ts_map_key,
                            pointing_key):
    '''
    Expands the G3MapVectorInt pointed to by pointing key.  It takes the pointing information mapped to by
    the key 'Boresight' and copies it to all the detectors listed in the keys of frame[ts_map_key]
    '''
    if frame.type != core.G3FrameType.Scan:
        return
    old_pointing = frame[pointing_key]
    new_pointing = core.G3MapVectorInt()
    for k in frame[ts_map_key].keys():
        new_pointing[k] = old_pointing['Boresight']
    del frame[pointing_key]
    frame[pointing_key] = new_pointing
    return

def FlattenBoresightPointing(frame, pointing_key):
    '''
    Converts the G3MapVectorInt with a single 'Boresight' entry pointing to by
    frame[pointing_key] into a single G3VectorInt of the same name.
    '''
    if pointing_key not in fr:
        return
    pointing = frame[pointing_key]['Boresight']
    del frame[pointing_key]
    frame[pointing_key] = pointing
    return

def AddUnityWeights(frame, ts_key, out_wgt_key):
    '''
    Creates a G3MapDouble and maps every timestream key in frame[ts_key] to 1.0.

    Input:
    ts_key -> G3TimestreamMap

    Ouput:
    out_wgt_key -> G3MapDouble
    '''
    if frame.type != core.G3FrameType.Scan:
        return
    frame[out_wgt_key] = core.G3MapDouble()
    if ts_key in frame:
        ks = frame[ts_key].keys()
        for k in ks:
            frame[out_wgt_key][k] = 1.0

def AddTodNoiseWeights(frame,
                       averaged_bin_num_to_use,
                       psd_key,
                       weight_key = 'TimestreamWeights',
                   ):
    '''
    This code has not been tested.  If you are using it because of the name and just think
    it will work, you should read documentation more.


    -This needs to be investigated...

    Sets the weight to be the inverse of the psd, unless it's 0, then it sets it to 0.
    Input:
      psd_key ->G3MapDouble
    Output:
      weight_key -> G3MapDouble
    '''
    if frame.type != core.G3FrameType.Scan:
        return
    ts_weights = core.G3MapDouble()
    psds = frame[psd_key]
    for k in psds.keys():
        if psds[k] == 0:
            ts_weights[k] = 0.0
        else:
            ts_weights[k] = 1.0/(psds[k][averaged_bin_num_to_use])
    frame[weight_key] = ts_weights

@core.indexmod
def AddEllBasedFilterSampleRate(frame,
                                boresight_az_key, boresight_el_key,
                                ell_filter_key):
    '''
    Estimates the scan speed from the angular distance between the starting and ending position of the boresight pointing
     and uses that to estimate the sample frequency of the TOD in units of "ell"
    '''
    if frame.type != core.G3FrameType.Scan:
        return
    if ell_filter_key in frame:
        return

    az = frame[boresight_az_key]
    el = frame[boresight_el_key]

    az_0 = az[0]/core.G3Units.rad
    az_1 = az[len(az)-1]/core.G3Units.rad

    el_0 = el[0]/core.G3Units.rad
    el_1 = el[len(az)-1]/core.G3Units.rad

    cart_0 = np.array( [np.cos(az_0) * np.cos(el_0), np.sin(az_0) * np.cos(el_0), np.sin(el_0)])
    cart_1 = np.array( [np.cos(az_1) * np.cos(el_1), np.sin(az_1) * np.cos(el_1), np.sin(el_1)])

    angular_dist = np.arccos(np.sum(cart_0 * cart_1)) * core.G3Units.rad
    frame[ell_filter_key] = 2.0*np.pi * len(az)/angular_dist

@core.usefulfunc
def covers_same_patch_of_sky(m0, m1):
    '''
    returns true if m0 and m1 cover the same patch of sky (checks that the shapes and coordinates and such are the same)
    '''
    if isinstance(m0, coordinateutils.G3SkyMapWeights):
        return covers_same_patch_of_sky(m0.TT, m1)
    if isinstance(m1, coordinateutils.G3SkyMapWeights):
        return covers_same_patch_of_sky(m0, m1.TT)
    if isinstance(m0, coordinateutils.G3SkyMap) and isinstance(m1, coordinateutils.G3SkyMap):
        return m0.IsCompatible(m1)

    raise RuntimeError("Map type not recognized")

@core.indexmod
class ExtractTheMaps(object):
    '''
    For any map frames coming through the pipe, caches them in self.maps.
    '''
    def __init__(self):
        self.maps = {}
    def __call__(self, frame):
        if frame.type == core.G3FrameType.Map:
            self.maps[frame["Id"]] = {}
            for m_id in ["T", "Q", "U", "Wpol", "Wunpol"]:
                if m_id in frame:
                    self.maps[frame["Id"]][m_id] = frame[m_id]

@core.indexmod
class MapInjector(object):
    '''
    MapInjector takes a list of maps (FlatSkyMap or HealpixSkyMap)
    (maps_lst) and adds them to a map frame with the proper format.  The Id of this map frame is map_id.

    Each MapInjector only adds one Map frame to the pipeline.

    is_stub should be set to False if the map you are injecting contains data like a point source
    mask or a simulated sky.  If it is the map being binned you should set it to True.

    If is_stub is False.  It stores the maps in maps_lst under the appropriate key in frame.
    If it is true, maps_lst, should be a list that only contains 1 element.
    That element is a G3SkyMap that used to specify the shape/location/resolution of the
    map for the frame. Empty maps are stored in the appropriate keys with the
    shape information specified by maps_lst[0].

    If make_polarized is True, T,Q, and U are stored.
    If do_weight is true, also adds weight maps

    For additional documentation see doc/mapmaking.rst
    '''
    def __init__(self, map_id, maps_lst,
                 is_stub = False,
                 make_polarized = True,
                 do_weight = True
             ):
        self.map_id = map_id
        self.has_emitted = False

        if is_stub:
            assert(len(maps_lst) == 1)
            map_in = maps_lst[0]
            if (make_polarized):
                T = map_in.Clone(False)
                T.pol_type = coordinateutils.MapPolType.T
                Q = map_in.Clone(False)
                Q.pol_type = coordinateutils.MapPolType.Q
                U = map_in.Clone(False)
                U.pol_type = coordinateutils.MapPolType.U
                self.maps_lst = [T,Q,U]
            else:
                T = map_in.Clone(False)
                T.pol_type = coordinateutils.MapPolType.T
                self.maps_lst = [T]
            if do_weight:
                wtype = coordinateutils.WeightType.Wpol if make_polarized else coordinateutils.WeightType.Wunpol
                wmap = coordinateutils.G3SkyMapWeights(map_in, wtype)
                self.maps_lst.append(wmap)
        else:
            assert(len(maps_lst) > 0)
            m0 = maps_lst[0]
            for m in maps_lst[1:]:
                if not covers_same_patch_of_sky(m0, m):
                    core.log_fatal("Maps passed to map injector do not cover the same patch of sky")
            self.maps_lst = maps_lst

    def __call__(self, frame):
        if not self.has_emitted:
            self.has_emitted = True
            new_frame = core.G3Frame(core.G3FrameType.Map)
            new_frame["Id"] = self.map_id
            for m in self.maps_lst:
                assert(isinstance( m, coordinateutils.G3SkyMap) or
                       isinstance( m, coordinateutils.G3SkyMapWeights))
                if (isinstance( m, coordinateutils.G3SkyMap)):
                    new_frame[str(m.pol_type)] = m
                else:
                    new_frame[str(m.weight_type)] = m
            return [ new_frame, frame ]
        else:
            return [ frame ]


@core.indexmod
def MapSplitter(frame, input_map_id, output_map_ids, also_split_weights):
    '''
    Takes an input map frame with Id input_map_id and outputs N map frames with the ids
    in output_map_ids.  The data stored in the maps in the input frame are not copied over
    just the patch of sky information.

    This is a helper module for making individual bolometer maps
    '''

    if frame.type != core.G3FrameType.Map:
        return
    if frame['Id'] == input_map_id:
        output_frames = [frame]

        is_first = True
        for oid in output_map_ids:
            o_frame = core.G3Frame(core.G3FrameType.Map)
            o_frame['Id'] = oid
            if also_split_weights or is_first:
                map_types = ['T', 'Q', 'U', 'Wpol', 'Wunpol']
                is_first = False
            else:
                map_types = ['T', 'Q', 'U']
            for mtype in map_types:
                if mtype in frame:
                    o_frame[mtype] = frame[mtype].Clone(False)
            output_frames.append(o_frame)
        return output_frames


@core.pipesegment_nodoc
def CalculatePointing(pipe,
                      map_id,
                      pointing_store_key, #output
                      trans_key,
                      ts_map_key="",
                      use_boresight_pointing = False,
                      bolo_props_key = 'BolometerProperties'):
    '''
    Calculates which pixel in the map the detectors are pointing at information for one map.

    Processing Arguments:
      map_id [string] : the id of the map to make maps for
      use_boresight_pointing [bool] : if true calculates boresight pointing and uses that for
          each detectors maps.
    Frame Input data:
      ts_map_key -> G3TimestreamMap: If ts_map_key =='', calculates the detector pointing
          for the detectors in the BolometerPropertiesMap at bolo_props_key
          If the ts_map_key is specified, calculates the pointing for
          the detectors in the timestream map
      bolo_props_key : BolometerProperties
      trans_key: The G3VectorQuat that maps (1,0,0) to our coordinate system of interest.
          You only need to specify this if you are also doing polarization angle rotation

    Frame Output Data:
      pointing_store_key -> G3VectorInt  where we want to store the G3MapVectorInt that record
          which pixel of the map each sample corresponds to.

    '''

    if use_boresight_pointing:
        return CalculateBoresightPointing(pipe, map_id=map_id,
                                          pointing_store_key=pointing_store_key,
                                          trans_key=trans_key,
                                          expand=True,
                                          bolo_map_key=ts_map_key or bolo_props_key)

    pipe.Add(CheckForCalFrames)

    pipe.Add(MapPointingCalculator, map_id = map_id,
             ts_map_key = ts_map_key,
             trans_key = trans_key,
             bolo_props_key=bolo_props_key,
             detector_alpha_key="",
             detector_delta_key="",
             pointing_store_key = pointing_store_key )


@core.pipesegment
def CalculateBoresightPointing(pipe,
                               map_id,
                               pointing_store_key, #output
                               trans_key,
                               expand=False,
                               bolo_map_key='BolometerProperties'):
    '''
    Calculates which pixel in the map the telescope boresight is pointing at
    information for one map.

    Processing Arguments
    --------------------
    map_id : [string]
        the id of the map to make maps for
    expand : [bool]
        if True, copy the Boresight pointing information to all channels
        mapped in the input G3Map object stored in frame[bolo_map_key].  Otherwise,
        store the boresight pointing information as a vector directly in
        frame[pointing_store_key].

    Frame Input data
    ----------------
    trans_key: The G3VectorQuat that maps (1,0,0) to our coordinate system of interest.
          You only need to specify this if you are also doing polarization angle rotation
    bolo_map_key: A G3Map object whose keys define which channels to use when copying the
          boresight pointing information.  Only used if expand is True.

    Frame Output Data
    -----------------
    pointing_store_key -> G3VectorInt  where we want to store the G3VectorInt that record
        which pixel of the map each sample corresponds to.
    '''

    pipe.Add(CheckForCalFrames)

    bolo_props_key = 'BoresightBoloProps'
    pointing_ts_key = 'BoresightTimestreams'
    pipe.Add(AddBoresightBoloProps, new_bp_key=bolo_props_key)
    pipe.Add(AddBsTimestream, tsm_key=pointing_ts_key)

    pipe.Add(MapPointingCalculator, map_id=map_id,
             ts_map_key=pointing_ts_key,
             trans_key=trans_key,
             bolo_props_key=bolo_props_key,
             detector_alpha_key="",
             detector_delta_key="",
             pointing_store_key=pointing_store_key)

    if expand:
        pipe.Add(ExpandBoresightPointing,
                 ts_map_key=bolo_map_key,
                 pointing_key=pointing_store_key)
    else:
        pipe.Add(FlattenBoresightPointing,
                 pointing_key=pointing_store_key)


@core.pipesegment_nodoc
def BinMap(pipe, map_id, ts_map_key,
           pointing_store_key,
           trans_key = '',
           timestream_weight_key = None,
           bolo_props_key = 'BolometerProperties',
           use_boresight_pointing = False,
           individual_bolos_to_map = None,
           force_making_sum_map = False,
           include_pol_rotation = False,
           use_unity_weights = False):
    '''
    Bins Tod into a Map.  When the map frame with ID map_id arrives it stores it in the frame.
    Whenever a scan frame passes by it bins the TOD into the map.  It will emit the map frame
    when it receives an end processing frame or if there is a Scan with the
    key "EmitMap" -> bool which is True

    Processing Arguments:

      map_id: The string Id of the map we are binning.  see the mapmaking.rst file

      use_boresight_pointing:  Whether we are using boresight pointing.  If we are
          we only need to store one weight
          map when making individual bolo maps

      individual_bolos_to_map:  A list of bolometer ids to generate individual maps for.
          If this is not None, we take the input map frame specified by map_id
          and make identical map frames with Id's that are the
          bolometer ids in the individual_bolos_to_map.  We then bin the
          respective timestreams in those maps.

    force_making_sum_map: if individual_bolos_to_map is not None it will
        skip making the coadd map. If you also want that coadd set this to True.
        Just to be clear this will be a coadd of all the detectors
        in the timestream map passed to this not just the ones in individual_bolos_to_map

    include_pol_rotation: boolean, if true includes the rotation from the coordinate transform
        in general this only needs to be true if you are making maps in galactic coordinates
    trans_key: The G3VectorQuat that maps (1,0,0) to our cooridnate system of interest.
        You only need to specify this if you are also doing polarization angle rotation

    Input Frame Keys:

    pointing_store_key -> G3MapVectorInt calculated in the CalculatePointing pipesegment

    timestream_weight_key ->G3MapDouble the weights assigned to each TOD for the map
        binning  This is an input to the code.

    '''

    if (individual_bolos_to_map == None):
        bolo_ids_to_map = core.G3VectorString()

    if use_unity_weights:
        timestream_weight_key = 'UnityTsWeights'
        pipe.Add(AddUnityWeights,
                 ts_key = ts_map_key,
                 out_wgt_key = timestream_weight_key)
    assert(not timestream_weight_key is None)


    if (not individual_bolos_to_map is None):
        pipe.Add(MapSplitter,
                 input_map_id = map_id,
                 output_map_ids = individual_bolos_to_map,
                 also_split_weights = not use_boresight_pointing)

    if ((individual_bolos_to_map is None) or force_making_sum_map):
        pipe.Add(MapBinner,
                 map_id = map_id,
                 ts_map_key = ts_map_key,
                 pointing_key = pointing_store_key,
                 weight_key = timestream_weight_key,
                 bolo_props_key = bolo_props_key,
                 trans_key = trans_key,
                 bolo_ids_to_map = core.G3VectorString(),
                 include_pol_rotation=include_pol_rotation )

    if (not individual_bolos_to_map is None):
        for map_id in individual_bolos_to_map:
            pipe.Add(MapBinner,
                     map_id = map_id,
                     ts_map_key = ts_map_key,
                     pointing_key = pointing_store_key,
                     weight_key = timestream_weight_key,
                     bolo_props_key = bolo_props_key,
                     trans_key = trans_key,
                     bolo_ids_to_map = core.G3VectorString([map_id]),
                     include_pol_rotation=include_pol_rotation )

@core.pipesegment_nodoc
def TodFiltering(pipe,
                 ts_in_key,
                 ts_out_key,

                 mhpf_cutoff = -1,
                 poly_order = -1,
                 lpf_filter_frequency = -1,
                 filters_are_ell_based = False,

                 point_source_mask_id = None,
                 point_source_pointing_store_key = None,

                 delete_input_ts = False,
                 keep_all_data = False,

                 use_dynamic_source_filter = False,
                 dynamic_source_filter_order = 1,
                 dynamic_source_filter_thresh = 4.0,
                 point_source_mask_debug = False,

                 boresight_az_key = 'BoresightAz',
                 boresight_el_key = 'BoresightEl',
                 filter_key = 'LowPassFilter',
                 ell_filter_effective_sample_rate_key = 'EllFilterSampleRate',
                 filter_mask_key = 'FilterMask',
                 fft_padding_key = 'FftPadding'):
    '''
    Filters the Tod.  Bam, comment done.

    But in all seriousness, this code applies a myriad of time ordered data filters to the data, specifically,
    a masked linear least squares (LLS) filter and FFT based low pass filter.  Remember, UNITS

    Processing Arguments:

    filters_are_ell_based:  If true, instead of using the units system, mhpf_cutoff and
      lpf_filter_frequency are specified in 'Ell' or spatial space.  This is estimated from the
      boresight pointing.

    mhpf_cutoff:  The high pass filter frequency cutoff for the LLS masked high pass filter, -1 means do not mhpf the data

    poly_order:  the order of polynomial to fit, -1 means don't do this

    lpf_filter_frequency: LPF cutoff frequency, -1 means don't lpf the data


    point_source_mask_id:  The Map Id of the point source mask. For the point source mask: 1 = masked, 0 = not masked.  This map should have been injected with the MapInjector


    point_source_pointing_store_key:  Where the pixel pointing (G3VectorInt) for
       the point source mask is stored.


    delete_input_ts:  Whether we should delete the input timestream.

    keep_all_data:  If True we keep a hold of all the intermediate data products, By that I mean, all of the copies of the timestreams from each filtering step.  This drastically increases memory usage but is super useful for debugging.


    use_dynamic_source_filter:  If true use dynamic filtering.
      Dynamic filtering applies a poly filter of dynamic_source_filter_order to the data
      Estimates the standard deviation of the signal with a MAD estimator and finds all the
      samples above dynamic_source_filter_thresh std devs and flags them.  This is not
      a linear operation so don't use it for mapmaking.

    Input Frame Data Keys:
      ts_in_key: Input key for the G3TimestreamMap

      boresight_az_key -> G3Timestream

      boresight_el_key -> G3Timestream

    Output Frame Data Keys:
      ts_out_key: Output key for the filtered G3TimestreamMap

    Frame Keys Used in the Mid Processing:
      filter_key -> G3VectorComplexDouble the frequency based filter

      ell_filter_effective_sample_rate_key -> G3Double the ell based sample rate for filtering

      filter_mask_key -> FilterMask  stores the masked pixels

      fft_padding_key -> G3Int the length of the padded fft for filtering.


    '''

    def delete_prev_ts_map(pipe, prev_key):
        if keep_all_data:
            return
        if delete_input_ts or prev_key != ts_in_key:
            pipe.Add(core.Delete, keys = [prev_key])
    prev_ts_map_key = ts_in_key
    run_fft_filter = lpf_filter_frequency > 0
    pipe.Add(lambda frame: (frame.type != core.G3FrameType.Scan or
                            len(frame[ts_in_key].keys()) > 0))

    if filters_are_ell_based:
        filter_sample_rate_key = ell_filter_effective_sample_rate_key
        pipe.Add(AddEllBasedFilterSampleRate,
                 boresight_az_key = boresight_az_key,
                 boresight_el_key = boresight_el_key,
                 ell_filter_key = ell_filter_effective_sample_rate_key)
    else:
        filter_sample_rate_key = None

    if (mhpf_cutoff > 0 or poly_order >= 0):
        is_masked = False
        if point_source_mask_id != None:
            is_masked = True
            #pipe.Add(core.InjectDebug, type = core.G3FrameType.Map)
            pipe.Add(todfilter.FilterMaskInjector,
                     point_src_mask_id = point_source_mask_id,
                     filter_mask_key = filter_mask_key,
                     pointing_key = point_source_pointing_store_key)
        elif use_dynamic_source_filter:
            core.log_debug("Using Dynamic Source")
            is_masked = True
            pipe.Add( todfilter.polyutils.DynamicSourceFiltering,
                      ts_map_key = prev_ts_map_key,
                      out_mask_key = filter_mask_key,
                      pre_filter_poly_order = dynamic_source_filter_order,
                      threshold = dynamic_source_filter_thresh)
        if point_source_mask_debug:
            core.log_warn("DOING DYNAMIC MASK DEBUGGING")
            poly_order = 0
            curr_ts_map_key = 'PntSrcMaskTod'
            pipe.Add(todfilter.polyutils.ReplaceTsMapWithDynamicMask,
                     in_dynamic_mask_key = filter_mask_key,
                     in_timestream_key  = prev_ts_map_key,
                     out_ts_map_key = curr_ts_map_key)
            delete_prev_ts_map(pipe, prev_ts_map_key)
            prev_ts_map_key = curr_ts_map_key
        if poly_order >= 0:
            curr_ts_map_key = 'PolyTsData'
            pipe.Add(todfilter.MaskedPolyHpf,
                     in_ts_map_key = prev_ts_map_key,
                     out_ts_map_key = curr_ts_map_key,
                     poly_order = poly_order,
                     high_pass_freq_cutoff = -1,
                     is_masked = is_masked,
                     mask_key = filter_mask_key,
                     sample_rate_override_key = ("" if filter_sample_rate_key is None
                                                 else filter_sample_rate_key)
                 )
            delete_prev_ts_map(pipe, prev_ts_map_key)
            prev_ts_map_key = curr_ts_map_key
        if mhpf_cutoff > 0:
            curr_ts_map_key = 'HpfTsData'
            pipe.Add(todfilter.MaskedPolyHpf,
                     in_ts_map_key = prev_ts_map_key,
                     out_ts_map_key = curr_ts_map_key,
                     poly_order = 0,
                     high_pass_freq_cutoff = mhpf_cutoff,
                     is_masked = is_masked,
                     mask_key = filter_mask_key,
                     sample_rate_override_key = ("" if filter_sample_rate_key is None
                                                 else filter_sample_rate_key)
                 )
            delete_prev_ts_map(pipe, prev_ts_map_key)
            prev_ts_map_key = curr_ts_map_key
    #Run the FFT filter (I know these comments are obvious, i mainly use the color coding to find things fast when looking)
    if (lpf_filter_frequency >= 0):
        pipe.Add(todfilter.dftutils.AddFftPaddingLengthKey, ts_map_key = prev_ts_map_key,
                 padding_length_key = fft_padding_key)

    if (lpf_filter_frequency >= 0):
        pipe.Add(todfilter.dftutils.LowPassFilterSpecifier,
                 ts_map_key = prev_ts_map_key,
                 sample_rate_override_key = filter_sample_rate_key,
                 input_filter_field = None,
                 output_filter_field = filter_key,
                 cutoff_freq = lpf_filter_frequency,
                 already_specified_key = 'LowPassFilterSpecified',
                 padding_key = fft_padding_key,
             )

    if (run_fft_filter):
        curr_ts_map_key = 'FFTFilteredTS'
        pipe.Add(todfilter.dftutils.FftFilter,
                 in_ts_key = prev_ts_map_key,
                 filter_path = filter_key,
                 out_ts_key = curr_ts_map_key
        )
        delete_prev_ts_map(pipe, prev_ts_map_key)
        prev_ts_map_key = curr_ts_map_key

    pipe.Add(core.Rename,
             keys = {prev_ts_map_key: ts_out_key},
             type = core.G3FrameType.Scan)

@core.pipesegment_nodoc
def FillSimTodSegment(pipe, sim_map_id, out_ts_key,
                      trans_key = 'RaDecRotation',
                      sim_pointing_key = 'SimPointing',
                      valid_ids = '',
                      interp_sim = False,
                      include_pol_rotation = False,
                      ts_to_get_sample_rate_key = 'OnlineBoresightRa',
                      bolo_props_key = 'BolometerProperties',
                      detector_alpha_pointing_key = '',
                      detector_delta_pointing_key = ''):
    '''
    '''
    pipe.Add(CheckForCalFrames)
    pipe.Add(core.Delete, type = core.G3FrameType.Scan, keys = [out_ts_key])

    if interp_sim:
        pipe.Add(SkyMapInterpTodFiller,
                 map_id = str(sim_map_id),
                 bolo_props_key = str(bolo_props_key),
                 ts_to_get_sample_rate_key = str(ts_to_get_sample_rate_key),
                 trans_key=str(trans_key),
                 include_pol_rot = bool(include_pol_rotation),
                 valid_ids = valid_ids,
                 detector_alpha_key = detector_alpha_pointing_key,
                 detector_delta_key = detector_alpha_pointing_key,
                 out_ts_key = str(out_ts_key))
    else:
        pipe.Add(SimulatedTimestreamFiller,
                 map_id = str(sim_map_id),
                 sim_pointing_key = str(sim_pointing_key),
                 bolo_props_key = str(bolo_props_key),
                 ts_to_get_sample_rate_key = str(ts_to_get_sample_rate_key),
                 ts_lst_key = str(valid_ids),
                 trans_key=str(trans_key),
                 include_pol_rot = bool(include_pol_rotation),
                 out_ts_key = str(out_ts_key))





def RejectTurnArounds(frame):
    return not (frame.type == core.G3FrameType.Scan and
        'Turnaround' in frame and frame['Turnaround'])

@core.usefulfunc
def remove_weight_t(Ti,W):
    '''
    Returns T,Q,U with the weight removed
    '''
    T = Ti.Clone(False)
    remove_weight_t_cpp(Ti,W,T)
    return T

@core.usefulfunc
def make_determinant_map(weight):
    f = weight.TT.Clone(False)
    make_determinant_map_cpp( weight, f)
    return f

def RemoveWeightModule(frame):
    ValidateMapFrames(frame)
    if 'Wpol' in frame:
        if frame['T'].is_weighted:
            T,Q,U = remove_weight(
                frame['T'], frame['Q'], frame['U'], frame['Wpol']
            )
            del frame['T']
            del frame['Q']
            del frame['U']
            frame['T'] = T
            frame['Q'] = Q
            frame['U'] = U
    elif 'Wunpol' in frame:
        if frame['T'].is_weighted:
            T = remove_weight_t(frame['T'], frame['Wunpol'])
            del frame['T']
            frame['T'] = T
    else:
        if frame['T'].is_weighted:
            raise TypeError('Map is weighted, but no weight is found')

@core.usefulfunc
def weights_cond(weight, mode=None):
    """
    Lifted with minor modification from Sasha's qpoint library

    Hits-normalized projection matrix condition number for
    each pixel.

    Arguments
    ---------
    weight : [G3SkyMapWeights]
        Weights object from which to calculate condition
    mode : [None, 1, -1, 2, -2, inf, -inf, 'fro'], optional
        condition number order.  See `numpy.linalg.cond`.
        Default: None (2-norm from SVD)

    Returns
    -------
    cond : array_like
        Condition number of each pixel.
    """

    # return if unpolarized
    if not weight.weight_type == coordinateutils.WeightType.Wpol:
        return weight.TT.Clone(True)

    cond = weight.TT.Clone(False)
    cond_array = np.asarray(cond)
    npix = weight.TT.shape[0] * weight.TT.shape[1]
    proj = np.zeros([6, npix])
    for i, m0 in enumerate(['TT', 'TQ', 'TU', 'QQ', 'QU', 'UU']):
        proj[i] = np.ravel(np.asarray(getattr(weight, m0)))

    nmap = 3
    nproj = len(proj)

    # normalize
    m = proj[0].astype(bool)
    proj[:, m] /= proj[0, m]
    proj[:, ~m] = np.inf

    # projection matrix indices
    idx = np.zeros((3, 3), dtype=int)
    rtri, ctri = np.triu_indices(nmap)
    idx[rtri, ctri] = idx[ctri, rtri] = np.arange(nproj)

    # calculate for each pixel
    # faster if numpy.linalg handles broadcasting
    npv = np.__version__.split('.')
    if len(npv) > 3:
        if 'dev' in npv[-1]:
            npv = [0,0,0]
        else:
            npv = np.version.short_version.split('.')
    npv = [int(x) for x in npv]
    if npv >= [1,10,0]:
        proj[:, ~m] = 0
        c0 = np.linalg.cond(proj[idx].transpose(2,0,1), p=mode)
        c0[~m] = np.inf
        # threshold at machine precision
        c0[c0 > 1./np.finfo(float).eps] = np.inf
        c0 = np.reshape(c0, weight.TT.shape)
        cond_array[:] = c0
        return cond

    # slow method, loop over pixels
    def func(x):
        if not np.isfinite(x[0]):
            return np.inf
        c = np.linalg.cond(x[idx], p=mode)
        # threshold at machine precision
        if c > 1./np.finfo(float).eps:
            return np.inf
        return c
    c0 = np.apply_along_axis(func, 0, proj)
    c0 = np.reshape(c0, weight.TT.shape)
    cond_array[:] = c0
    return cond


def ZeroMapNans(frame):
    '''
    Turns NaNs (and infs) in maps in a Map frame to 0
    '''
    np.asarray(frame['T'])[np.invert(np.isfinite(frame['T']))] = 0
    if 'Wpol' in frame:
        np.asarray(frame['Q'])[np.invert(np.isfinite(frame['Q']))] = 0
        np.asarray(frame['U'])[np.invert(np.isfinite(frame['U']))] = 0
