from spt3g import core

from spt3g.coordinateutils.azel import convert_azel_to_radec
from spt3g.pointing.offline_pointing import CorrectBoresightPointing

from spt3g.coordinateutils import create_det_az_el_trans, create_lazy_det_ra_dec_trans
from spt3g.coordinateutils import create_det_ra_dec_trans, convert_ra_dec_trans_to_gal

from copy import copy
import numpy as np

@core.indexmod
def CalcAzElToRaDecOffsetFromBoresight(
        frame, 
        correct_bs_pointing_flags = ['az_tilts', 'el_tilts', 'flexure', 'collimation',
                                     'refraction'],
        bs_az_key = 'RawBoresightAz', bs_el_key='RawBoresightEl',
        offset_az_key='OffsetBoresightAz', offset_el_key='OffsetBoresightEl', 
        offset_ra_key='OnlineOffsetRa', offset_dec_key='OnlineOffsetDec',
        offset_val = 5 * core.G3Units.arcmin):
    '''
    Repeats the offline pointing model calculation for an az/el coordinate slightly offset
    from the boresight pointing.

    Inputs:
    bs_az_key, bs_el_key : raw boresight az and el

    offset_val : the angular distance of the offset coordinate, this should be small

    Outputs:
    offset_az_key, offset_el_key : our offset values in local coordinates 
    offset_ra_key, offset_dec_key : our calculated offline pointing ra/dec


    ::::::WARNING::::::
    Right now this code repeats the calculation of the boresight pointing on the offset 
    points.  This is *I think* true with the way the code is implemented now, but it is a 
    huge abstraction violation since the boresight pointing code is only meant to calculate
    things for boresight.

    If the boresight pointing code is modified this code will need to be modified.
    '''


    if frame.type != core.G3FrameType.Scan:
        return
    az = frame[bs_az_key]
    el = copy(frame[bs_el_key])

    #with some slop, if we are near the pole the ell offset can do some weird things
    np.asarray(el)[np.where(el > (np.pi - 2*offset_val))] -= 2 * offset_val
    el += offset_val
    frame[offset_az_key] = az
    frame[offset_el_key] = el

    #so we need to calculate the Real Az El for a detector offset from boresight pointing.
    #This bit of code is violating the abstraction barrier of the CorrectBoresightPointing
    #in that, CorrectBoresightPointing is only meant for applying boresight corrections not 
    #corrections for detectors off of boresight.
    #We have taken a peak under the hood
    cbp = CorrectBoresightPointing(raw_az_key=offset_az_key, raw_el_key=offset_el_key,
                                   flags = correct_bs_pointing_flags,
                                   output='OffsetReal')
    cbp(frame)
    ra, dec = convert_azel_to_radec(frame['OffsetRealAz'], frame['OffsetRealEl'])
    frame[offset_ra_key] = ra
    frame[offset_dec_key] = dec
