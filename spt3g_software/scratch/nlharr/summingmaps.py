from spt3g import todfilter, core
import numpy as np
from copy import copy


'''
This code mostly existed to be hacky fixes before the pointing code got implemented.  Some of the 
functions in this file will probably be useful in general.
'''



def find_brightest_pixel(m):
    '''
    returns the brightest pixel in a map...
    '''
    return np.asarray(map(lambda x: x[0], np.where(m == np.max(m))))

def shift_map(m, offset):
    '''
    Shifts a 2d sky map by the value specified in offset.
    Offset has the form [y,x]
    '''
    m_shift = np.roll(np.roll(m, offset[0], 0), offset[1], 1)
    if offset[0] > 0:
        m_shift[ : offset[0], : ] = 0
    elif offset[0] < 0:
        m_shift[ offset[0]:, : ] = 0
    if offset[1] > 0:
        m_shift[:, : offset[1]] = 0
    elif offset[1] < 0:
        m_shift[:, offset[1]:] = 0
    out_m = copy(m)
    np.asarray(out_m)[:] = m_shift

    return out_m

def get_brightest_pixel_offset(m, desired_pos):
    '''
    returns the offset usable by shift_map for moving the map to the brightest pixel
    '''
    return desired_pos - find_brightest_pixel(m)

def DeleteWeights(frame):
    '''
    Removes the weights keys from the map frames com
    '''

    if frame.type != core.G3FrameType.Map:
        return
    ks = ['Wpol', 'Wunpol']
    for k in ks:
        if k in frame:
            del frame[k]

#removed weight first
class OffsetTMapFromPreviousOffset(object):
    '''
    This module exists for the hack cena analysis with bad pointing.
    This module aligns the individual detector maps with the coadd and
    centers all of them at the brightest pixel.

    When the frame order looks like:
    [Coadd Map Frame, individual map frame 0, individual map frame 1, ...]
    and you pass the offset_map_id as the id of the coadd map 
    to this module, it will find the offset needed to
    shift the map to the desired_center_pixel from the coadd, and then shift all
    the individual detector maps following.
    
    I don't know, you try documenting hack fix code.  Jerk.
    '''
    def __init__(self, offset_map_id, desired_center_pixel, 
                 center_pixel_ra_key = 'BrightestPixelRaOrig',
                 center_pixel_dec_key = 'BrightestPixelDecOrig',
             ):
        self.offset_map_id = offset_map_id
        self.cached_offset = None
        self.desired_center_pixel = desired_center_pixel
        self.center_pixel_ra_key = center_pixel_ra_key
        self.center_pixel_dec_key = center_pixel_dec_key

    def __call__(self, frame):
        if frame.type != core.G3FrameType.Map:
            return
        if frame['Id'] == self.offset_map_id:
            brightest_pixel = find_brightest_pixel(frame['T'])
            ra_dec = frame['T'].pixel_to_angle(brightest_pixel[1],
                                               brightest_pixel[0])
            self.cached_offset = (np.asarray(self.desired_center_pixel) 
                                  - brightest_pixel)
            frame[self.center_pixel_ra_key] = ra_dec[0]
            frame[self.center_pixel_dec_key] = ra_dec[1]
        if not self.cached_offset is None:
            T = shift_map(frame['T'], self.cached_offset)
            del frame['T']
            frame['T'] = T

def estimate_map_noise_variance(m, mask):
    '''
    estimates the variance using todfilter.get_mad_std of the map m wherever the mask is 1.
    '''
    mm = np.nan_to_num( m * (1-mask))
    return todfilter.get_mad_std(core.G3Timestream(mm[np.where(mm != 0)]))**2.0

def get_map_weight_from_map(m, mask):
    '''
    creates a weight map from our map with noise estimated from the map
    '''
    scalar_variance = estimate_map_noise_variance(m, mask)
    if scalar_variance == 0:
        return m * 0, 0
    else:
        return np.abs(np.sign(m)) * (1.0/scalar_variance), scalar_variance

def MakeMapsWeightedFromMap(frame, mask):
    '''
    This is the module form of get_map_weight_from_map.

    Creates the weight maps from the maps coming through.
    '''
    if frame.type != core.G3FrameType.Map:
        return
    mp = copy(frame['T'])
    del frame['T']
    if 'Wunpol' in frame:
        del frame['Wunpol']
    wunpol, w_scale = get_map_weight_from_map(mp, mask)
    weight = core.G3SkyMapWeights(mp, weight_type=core.WeightType.Wunpol)
    np.asarray(weight.TT)[:] = wunpol
    frame['Wunpol'] = weight
    mp.is_weighted = True
    np.asarray(mp)[:] = mp * wunpol
    frame['T'] = mp


class AddUnpolarizedMapsComingThrough(object):
    def __init__(self):
        self.map_dic = {}
        self.map_types_to_add = ['T', 'Wunpol']

    def __call__(self, frame):
        if frame.type == core.G3FrameType.EndProcessing:
            frames_to_return = []
            for k, v in self.map_dic.items():
                fr = core.G3Frame(core.G3FrameType.Map)
                fr['Id'] = k
                for i, t_key in enumerate(self.map_types_to_add):
                    fr[t_key] = v[i]
                frames_to_return.append(fr)
            frames_to_return.append(frame)
            return frames_to_return
        if frame.type == core.G3FrameType.Map:
            mid = frame['Id']
            if not mid in self.map_dic:
                self.map_dic[mid] = []
                for t_key in self.map_types_to_add:
                    self.map_dic[mid].append(frame[t_key])
            else:
                for i, t_key in enumerate(self.map_types_to_add):
                    self.map_dic[mid][i] += frame[t_key]
            return []

@core.cache_frame_data(type = core.G3FrameType.Map, bolo_props = 'BolometerProperties')
def FilterMapsNotInBoloProps(frame, bolo_props = None):
    if not frame['Id'] in bolo_props:
        #print frame['Id'], not frame['Id'] in bolo_props
        return []
    else:
        return

@core.cache_frame_data(type = core.G3FrameType.Map, bolo_props = 'BolometerProperties')
def FilterMapToBand(frame, band = 150 * core.G3Units.GHz, bolo_props = None):
    if frame.type != core.G3FrameType.Map:
        return
    if bolo_props[frame['Id']].band != band:
        return []
    else:
        return

class AddAllTheMaps(object):
    def __init__(self):
        self.t = None
        self.q = None
        self.u = None
        self.w = None
    def __call__(self, frame):
        if frame.type == core.G3FrameType.Map:
            if self.t is None:
                self.t = frame['T']
                self.q = frame['Q']
                self.u = frame['U']
                self.w = frame['Wpol']
            else:
                self.t += frame['T']
                self.q += frame['Q']
                self.u += frame['U']
                self.w += frame['Wpol']
            return []
        if frame.type == core.G3FrameType.EndProcessing:
            fr = core.G3Frame(core.G3FrameType.Map)
            fr['Id'] = 'SumMap'
            fr['T'] = self.t
            fr['Q'] = self.q
            fr['U'] = self.u
            fr['Wpol'] = self.w
            return [fr, frame]
