import numpy as np
from spt3g import core, calibration

def RemoveAzGlitches(frame, LoCutoff = 0.00004*core.G3Units.rad,
    HiCutoff = 6*core.G3Units.rad):
    '''
    Drops frames containing az-following glitches from the mapmaking pipeline.
    The default values for how big of a discrepancy qualifies as a glitch have
    been determined by inspection, and can be changed as an input argument
    (though this is not recommended). The search for glitches begins at count
    300 of each scan in order to not pick up any ringdown effects from the
    turnaround.
    '''
    if 'TrackerStatus' not in frame:
        return
    status = frame['TrackerStatus']
    diff = np.abs(status.az_command - status.az_pos)
    dif = diff[300:len(diff)]
    if np.any((dif > LoCutoff) & (dif < HiCutoff)):
        return False

def CutOnScanSpeed(frame, boresight_az_key = 'OnlineBoresightAz', 
                   min_speed = 0 * core.G3Units.deg/core.G3Units.sec,
                   max_speed = 1.95 * core.G3Units.deg/core.G3Units.sec):
    '''
    Drop frames during which the on-bearing azimuth speed drops below
    min_speed or exceeds max_speed.
    Mostly useful for dropping scans containing an az wrap/unwrap
    '''
    if boresight_az_key not in frame:
        return
    az_speed = (np.abs(np.gradient(frame[boresight_az_key]))
                * frame[boresight_az_key].sample_rate)
    # using indices of sorted array instead of min, max for some wiggle room
    # and to avoid e.g. jump going from 0 deg to 360 deg
    if sorted(az_speed)[5]<min_speed or sorted(az_speed)[-5]>max_speed:
        return False

def RemoveShortScans(frame, nsamples = 1000):
    '''
    Drops anomalously short Scan frames.
    Default cut somewhat arbitrarily set to 10 seconds.
    '''
    if 'Turnaround' in frame or 'TrackerStatus' not in frame:
        return
    if len(frame['TrackerStatus'].az_pos) < nsamples:
        return False

def RemoveSparseScans(frame, ts_key, min_num = 1):
    '''
    Drops frames in which the specified ts_key is not in the frame
    or contains fewer than min_num elements.
    '''
    if frame.type != core.G3FrameType.Scan:
        return
    if ts_key not in frame or len(frame[ts_key].keys()) < min_num:
        return False
