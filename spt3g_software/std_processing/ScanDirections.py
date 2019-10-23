from spt3g import core
import numpy


@core.indexmod
def SelectScanDirection(fr, left=False, az='OnlineBoresightAz'):
    '''
    Select only scans moving in a particular direction (left, to lower az, or
    right, to increasing az). Internally drops turnarounds, since they have
    no sensible interpretation here. If 'left' is True, keeps only
    left-going scans; otherwise, only right-going.
    '''

    if fr.type != core.G3FrameType.Scan:
        return None
    if fr.get('Turnaround', False):
        return False

    azchange = numpy.diff(fr[az])
    if left and sum(azchange > 0) > sum(azchange < 0):
        return False
    if (not left) and sum(azchange > 0) < sum(azchange < 0):
        return False
    
@core.indexmod
def split_left_right_scans(fr, ts_in_key = 'DeflaggedTimestreams',
                           boresight_az_key ='OnlineBoresightAz',
                           ts_left_key = None, ts_right_key = None):
    '''
    Splits scans into left-going (decreasing Az) and right-going (increasing Az).
    `ts_in_key` will be renamed to `ts_in_key`Left for left-going scans,
    or `ts_in_key`Right for right-going scans if `ts_left_key` or
    `ts_right_key` is None.  Otherwise, `ts_left(right)_key` will contain
    left(right)-going scans.
    '''
    if fr.type != core.G3FrameType.Scan:
        return
    if ts_left_key is None or ts_right_key is None:
        ts_left_key = ts_in_key + 'Left'
        ts_right_key = ts_in_key + 'Right'
    daz = numpy.gradient(fr[boresight_az_key])
    ts_out_key = ts_left_key if daz[0] < 0 else ts_right_key
    fr[ts_out_key] = fr[ts_in_key]
