import numpy as np
import scipy.stats
from copy import deepcopy
from spt3g import core
from spt3g.util.genericutils import sort_two_lists
from spt3g.calibration.template_groups import get_template_groups

'''
Flags are stored in a G3MapVectorString that has the form: 
  {bolometer_id: ['BadCalibrator', 'Glitchy', 'ActuallyAKitten']}

The bolometer id maps to a list of reasons why the bolometer is flagged.  

When developing code to flag timestreams we are trying to not force 
flagging decisions on people using intermediate data products.
By that I mean, if we need to compute something in order to flag the timestreams, 
we should compute the thing, store it in the frame, and then use a separate, 
later, Module to actually set the flags. 

As an example, with the glitch finder, we have a module that finds the number of 
glitches above a threshold, and a seperate module for flagging based off of those stored values.

Remember that data products in the frames are supposed to be immutable.  
Do *NOT* do something like frame['Flags'][bolo_id].append("NoisyTimestream").  
BAD BAD BAD, NO. SHAME. BAD

You will need to copy the flags, append the flag you want, 
delete the previous flags in the frame and  then add the new flags.  
Because having you do that every time would be insane please use the add_flag
convenience function.
'''

@core.usefulfunc
def add_flag(frame, flag_key, flag_id, timestreams_to_flag):
    '''
 Adds a flag to the timestreams listed in timestreams_to_flag.

    frame : [G3Frame]
        The frame data
    flag_key : [G3MapVectorString]
        The key in the frame where the flags object is 
        stored, eg "Flags"
    flag_id : [string]
        The string labelling the flag we are 
        adding, eg 'Glitchy', 'BadCalibratorResponse'
    timestreams_to_flag : array-like
        A   list of timestreams to flag.

    '''
    if flag_key not in frame:
        flags = core.G3MapVectorString()
    else:
        flags = deepcopy(frame.pop(flag_key))
    for k in timestreams_to_flag:
        if not k in flags:
            flags[k] = core.G3VectorString()
        flags[k].append(flag_id)
    frame[flag_key] = flags

@core.indexmod
def RemoveFlaggedTimestreams(frame, 
                             input_ts_key, 
                             input_flag_key, 
                             output_ts_key):
    '''
    Flags are stored in a G3MapVectorString that has the form: {bolometer_id: ['BadCalibrator', 'Glitchy', 'ActuallyAKitten']}
    
    The bolometer id maps to a list of reasons why the bolometer is flagged.  

    input_ts_key [->G3TimestreamMap] input timestreams
    input_flag_key [->G3MapVectorString] specified flags
    output_ts_key [->G3TimestreamMap] output timestream
    
    '''

    if frame.type != core.G3FrameType.Scan:
        return
    ts_map = core.G3TimestreamMap(frame[input_ts_key])
    flags = frame[input_flag_key]

    core.log_debug('Flagging: %d'%len(flags.keys()))

    for k in flags.keys():
        assert(len(flags[k]) > 0)
        ts_map.pop(k, None)
    frame[output_ts_key] = ts_map

@core.scan_func_cache_data(m_key=None)
def FlagBadG3MapValue(frame, m_key=None, 
                      flag_key='Flags', flag_reason='NoReasonGiven', 
                      min_val = None, max_val = None):
    '''
    Flags the values of a G3MapDouble mapped to by m_key that are below the min val
     or above the max val.  This will grab the flag_key information from scan frames
    or, if the key is present in other frames, will cache the value.

    m_key [->G3MapDouble] the map we are doing flagging on.  Keys are bolo ids
    min_val [float] if frame[m_key][bolo_id] < min_val it's flagged
    max_val [float] if frame[m_key][bolo_id] > min_val it's flagged
    frame_type: the type of frame to run on

    flag_key [G3MapVectorString] where the flags are stored
    flag_reason [string] the string identifying why these detectors are being flagged

    '''
    if m_key is None:
        core.log_error("m_key is None", flag_reason)
    m = m_key
    bad_lst = []
    for k in m.keys():
        if (((not min_val is None) and (np.isnan(m[k]) or m[k] < min_val)) or
            ((not max_val is None) and (np.isnan(m[k]) or m[k] > max_val))):
            bad_lst.append(k)
    add_flag(frame, flag_key, flag_reason, bad_lst)

FlagBadG3MapValueCal = FlagBadG3MapValue

@core.scan_func_cache_data(m_key = None)
def SigmaclipFlagG3MapValue(frame, m_key = None, low = 4.0, high = 4.0,
                            flag_key = 'Flags', flag_reason = 'NoReasonGiven'):
    '''
    Flag values in a G3MapDouble mapped by m_key that are more than
    the supplied number of standard deviations from the mean.  
    Standard deviation and mean are iteratively calculated using
    scipy.stats.sigmaclip.
    
    INPUTS
    ------
    m_key: str
        the key on which to calculate statistics and apply the flagging
    low: float
    high: float
        The number of sigma at which to flag values.  Any values between `low`
        and `high` will NOT be flagged.
    flag_key: str ['Flags']
        The key under which flags are stored
    flag_reason: str ['NoReasonGiven']
        A string identifier indicating why this flag has been applied.
    '''
    if m_key is None:
        core.log_error("m_key is None", flag_reason)
    m = m_key
    bolos = np.array(m.keys())
    values = m.values()
    clipper = scipy.stats.sigmaclip(values, low, high)
    bad = np.where((values > clipper.upper) | (values < clipper.lower))[0]
    add_flag(frame, flag_key, flag_reason, bolos[bad])

@core.scan_func_cache_data(bolo_props = 'BolometerProperties', 
                           wiring_map = 'WiringMap', m_key = None)
def SigmaclipFlagGroupG3MapValue(frame, m_key = None, low = 4.0, high = 4.0,
                                 per_band = True, per_wafer = False, 
                                 per_squid = False, 
                                 flag_reason = 'NoReasonGiven',
                                 flag_key = 'Flags',
                                 wiring_map = None, bolo_props= None):
    '''
    Flag values by group in a G3MapDouble mapped by m_key that are more than
    the supplied number of standard deviations from the mean.  
    Standard deviation and mean are iteratively calculated using
    scipy.stats.sigmaclip.  Groups are either bands, squids or wafers.
    
    INPUTS
    ------
    m_key: str
        the key on which to calculate statistics and apply the flagging
    low: float
    high: float
        The number of sigma at which to flag values.  Any values between `low`
        and `high` will NOT be flagged.
    flag_key: str ['Flags']
        The key under which flags are stored
    flag_reason: str ['NoReasonGiven']
        A string identifier indicating why this flag has been applied.
    per_band, pre_wafer, per_squid: bool
        Split by band (default), wafer, or SQUID.
    bolo_props, wiring_map are values cached from Calibration frames.  
        In general, there is no need to change these.
    '''
    if m_key is None:
        core.log_error("m_key is None", flag_reason)
    m = m_key
    template_groups = get_template_groups(bolo_props, wiring_map = wiring_map,
                                          per_band = per_band, 
                                          per_wafer = per_wafer,
                                          per_squid = per_squid, 
                                          include_keys = True )
    for dev_id, bolos in template_groups.items():
        # Build a list of bolos in this group and in the key, in case 
        # bolos have been cut already
        bolos = set(bolos).intersection(set(m.keys()))
        bolos = np.array(list(bolos))
        values = [m[bolo] for bolo in bolos]
        clipper = scipy.stats.sigmaclip(values, low, high)
        bad = np.where((values > clipper.upper) | (values < clipper.lower))[0]
        add_flag(frame, flag_key, flag_reason, bolos[bad])
        
@core.scan_func_cache_data(bolo_props = 'BolometerProperties', wiring_map = 'WiringMap')
def GroupAverageG3MapValue(frame, input_g3_map_key, output_g3_map_key,
                           average_func = np.mean,
                           per_band = True, per_wafer = False,  per_squid = False, 
                           cache_device_settings_key = None,

                           wiring_map = None, bolo_props= None):
    '''
    Averages a G3Map with bolometers as keys over specified groups, where the groups are
    done per band, wafer or squid (or any combination of them.
    
    input_g3_map_key ->[G3MapDouble] The G3MapDouble with keys of bolometer ids
    output_g3_map_key ->[G3MapDouble] The averaged values.   The keys are stringified
                    versions of the groups we are splitting on.

    per_band, per_wafer, per_squid : They are all bools.  If they are true the division
                        is split on those values.  If multiple are true, it splits on the 
                        combination of the values.  So if per_band and per_squid are true
                        it will divide the bolometers into groups where only bolometers
                        with the sames squid and band are grouped.

    if cache_device_settings_key is a string it will store the per_band, per_wafer and per_squid
    values in the frame

    '''

    assert(not bolo_props is None)
    mp = frame[input_g3_map_key]
    out_mp = core.G3MapDouble()
    template_groups = get_template_groups(bolo_props, wiring_map = wiring_map,
                                          per_band = per_band, per_wafer = per_wafer,
                                          per_squid = per_squid, include_keys = True )
    for dev_id, glst in template_groups.iteritems():
        av_val = average_func([mp[k] for k in glst])
        out_mp[dev_id] = av_val

    if not cache_device_settings_key is None:
        frame[cache_device_settings_key] = core.G3VectorInt([per_band, per_wafer,  per_squid])

    frame[output_g3_map_key] = out_mp

@core.scan_func_cache_data(bolo_props = 'BolometerProperties', wiring_map = 'WiringMap')
def FlagGroupAverageG3MapValue(frame, g3_map_key,
                               min_val = None, max_val = None,
                               flag_reason = 'Unspecified', flag_key = 'Flags',
                               wiring_map = None, bolo_props = None,
                               per_band = True, per_wafer = False, per_squid = False,
                               cache_device_settings_key = None):

    '''
    Flags timestreams based off the values averaged with GroupAverageG3MapValue
    
    
    g3_map_key ->[G3MapDouble] the output map from GroupAverageG3MapValue
    per_band, per_wafer, per_squid : They are all bools.  If they are true the division
                        is split on those values.  If multiple are true, it splits on the 
                        combination of the values.  So if per_band and per_squid are true
                        it will divide the bolometers into groups where only bolometers
                        with the sames squid and band are grouped.

    cache_device_settings_key: Used the match the cache_device_setting_key in 
         GroupAverageG3MapValue.  If provided overrides the per_bad, per_wafer, per_squid
         settings.

    '''

    assert(not bolo_props is None)


    if not cache_device_settings_key is None:
        dev_settings = frame[cache_device_settings_key]
        assert(len(dev_settings) == 3)
        per_band = dev_settings[0]
        per_wafer = dev_settings[1]
        per_squid = dev_settings[2]

    bad_bolos = []
    mp = frame[g3_map_key]
    out_mp = core.G3MapDouble()
    template_groups = get_template_groups(bolo_props,
                                          wiring_map = wiring_map,
                                          per_band = per_band,
                                          per_wafer = per_wafer,
                                          per_squid = per_squid,
                                          include_keys = True  )
    for dev_id, glst in template_groups.iteritems():
        is_bad = False
        if (not min_val is None) and mp[dev_id] < min_val:
            is_bad = True
        if (not max_val is None) and mp[dev_id] > max_val:
            is_bad = True
        if is_bad:
            bad_bolos += list(glst)
    add_flag(frame, flag_key, flag_reason, bad_bolos)


class GenerateFlagStats(object):
    def __init__(self,flag_key='Flags', verbose=True):
        self.flag_key = flag_key
        self.verbose = verbose
        self.scan_stats = []
        self.scan_bolos = []
        self.total_stats = {}
    def __call__(self, frame):
        if frame.type == core.G3FrameType.EndProcessing:
            if self.verbose:
                print(self.get_string_stats())
        if frame.type != core.G3FrameType.Scan:
            return
        flags = frame[self.flag_key]
        cut_stat = {}
        cut_bolo = {}
        for k, v in flags.iteritems():
            for reason in v:
                cut_stat[reason] = cut_stat.get(reason, 0) + 1
                if reason not in cut_bolo:
                    cut_bolo[reason] = [k]
                else:
                    cut_bolo[reason].append(k)
        self.scan_stats.append(cut_stat)
        self.scan_bolos.append(cut_bolo)
        for k in cut_stat:
            self.total_stats[k] = self.total_stats.get(k,0)+cut_stat[k]
    def get_string_stats(self, frame_num = None):
        o_s = '******By Scan******\n'
        for i,d in enumerate(self.scan_stats):
            ks = d.keys()
            vs = d.values()
            ks, vs = sort_two_lists(ks,vs)
            o_s += 'Scan %d:\t'%i
            for k,v in zip(ks,vs):
                o_s += '%s\t%d   '%(k,v)
            o_s += '\n'
        o_s += '\n\n******Total Flag Reasons ****** \n'
        d = self.total_stats
        ks = d.keys()
        vs = d.values()
        vs, ks = sort_two_lists(vs,ks)
        for k,v in zip(ks,vs):
            o_s += '%s\t%d\tAvg per Scan:\t%d\n'%(k,v,v/(i+1))
        return o_s
    def update(self, frame):
        if frame.type != core.G3FrameType.Scan:
            return
        flags = frame[self.flag_key]
        cut_stat = {}
        cut_bolo = {}
        old_count = {}
        for k, v in flags.iteritems():
            for reason in v:
                cut_stat[reason] = cut_stat.get(reason, 0) + 1

                if reason not in self.scan_bolos[-1]:
                    self.scan_bolos[-1][reason] = [k]
                else:
                    if k not in self.scan_bolos[-1][reason]:
                        self.scan_bolos[-1][reason].append(k)

        for reason, val in cut_stat.items():
            if reason not in self.scan_stats[-1]:
                old_count[reason] = 0
                self.scan_stats[-1][reason] = val
            if self.scan_stats[-1][reason] < val:
                old_count[reason] = self.scan_stats[-1][reason]
                self.scan_stats[-1][reason] = val
            else:
                old_count[reason] = val
        for k in cut_stat:
            self.total_stats[k] = self.total_stats.get(k,0)+cut_stat[k]-old_count[k]







