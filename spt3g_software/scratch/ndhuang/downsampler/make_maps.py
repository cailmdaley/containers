import os
import numpy as np
from spt3g import core, mapmaker, sptpol
from downsample import Downsampler, fancy_downsample

def make_sim_map_np(center, like_el = True):
    '''
    500 d obs, proj 0 --> 20.0 d x 47.0 d
    iceboard samplerate 152 Hz
    2 arcmin pixels in the map
    for now filter comparisons, make pixels of the same size
    for powspec, small pixels
    '''
    m_el = np.zeros((600 * 20, 600 * 47))
    el = np.linspace(center[1] - 10, center[1] + 10, np.shape(m_el)[0])
    ra = np.linspace(center[0] - 47./2, center[0] + 47./2, np.shape(m_el)[1])
    if like_el:
        for row in xrange(np.shape(m_el)[0]):
            m_el[row] = el[row]
        return m_el
    for row in xrange(np.shape(m_el)[0]):
        m_el[row] = 1 / np.sin(np.deg2rad(el[row]))
    m_ra = np.zeros_like(m_el)
    for col in xrange(np.shape(m_ra)[1]):
        m_ra[:, col] = (ra[col] - ra[0]) / (ra[0] - ra[-1]) * .01 + .99
    return m_el * m_ra

def make_sim_maps(center, like_el = True):
    m = make_sim_map_np(center, like_el)
    simmap = coordinateutils.FlatSkyMap(m,
                                 res = .1 * core.G3Units.arcmin,
                                 is_weighted = False,
                                 proj = mapmaker.MapProjection.Proj0,
                                 alpha_center = center[0] * core.G3Units.deg, 
                                 delta_center = center[1] * core.G3Units.deg,
                                 pol_type = core.MapPolType.T)
    return simmap

def interpolate_ptng(fr):
    if not fr.type == core.G3FrameType.Scan:
        return
    az = fr['BoresightAz']
    el = fr['BoresightEl']
    ra = fr['BoresightRa']
    dec = fr['BoresightDec']
    for key in ['BoresightAz', 'BoresightEl', 'BoresightRa', 'BoresightDec']:
        vec = fr[key]
        t = np.arange(len(vec))
        tslow = np.linspace(0, len(vec), 20./25 * len(vec))
        vecslow = np.interp(tslow, t, vec)
        fr[key + 'Slow'] = core.G3Timestream(vecslow)
    return fr

def uniform_ts_weights(fr):
    if not fr.type == core.G3FrameType.Scan:
        return
    ts_weight = fr.pop('TimestreamWeights')
    for bolo, weight in ts_weight:
        if weight != 0:
            ts_weight[bolo] = 1.0
        else:
            ts_weight[bolo] = 0.0
    fr['TimestreamWeights'] = ts_weight
    

def debugger(fr, debug = True):
    if debug:
        import IPython
        IPython.embed()

def periodic_dump(fr):
    if fr.type == core.G3FrameType.Scan:
        if fr['ScanNumber'] % 10 == 0:
            print fr

def delete_mapping_cache(fr):
    for k in ['MapPointing', 'DetectorAlphaPointing', 'DetectorDeltaPointing']:
        # also stokes
        fr.pop(k, None)

def dosim():
    savedir = '/home/ndhuang/spt/data/dec_fancyconvolve_elmap_3x'
    # savedir = '/data/ndhuang/test/downsampler/dec_convolve'
    try:
        os.makedirs(savedir)
    except OSError:
        pass
    # set up the maps, make the sim atmosphere
    center = [0., -57.5]
    bs_ra = 'BoresightRa'
    bs_dec = 'BoresightDec'
    suffix = '_downsampled'
    sim_map = make_sim_maps(center)
    ds_factor = 3
    # sim_map = np.load('/data/ndhuan/test/downsampler/sim_map_big.npy')
    map_params = coordinateutils.FlatSkyMap(x_len = 47 * 30, y_len = 20 * 30,
                                     res = 2 * core.G3Units.arcmin,
                                     alpha_center = center[0] * core.G3Units.deg,
                                     delta_center = center[1] * core.G3Units.deg,
                                     proj = mapmaker.MapProjection.Proj0)
    # Set up modules
    map_extractor = mapmaker.mapmakerutils.ExtractTheMaps()
    ds = Downsampler(ds_factor, key_suffix = suffix, keys = ('sim_ts',), 
                     decimate_keys = (bs_ra, bs_dec),
                     debug = True)
    #pipeit
    pipe = core.G3Pipeline()
    pipe.add = pipe.Add
    pipe.Add(sptpol.DirectIdfReader, 
             # filename = '/data/ndhuang/test/ds_idf_20150717_005656_150ghz.h5',
             filename = '/home/ndhuang/spt/data/ds_idf_20150717_005656_150ghz.h5',
             load_bolo_data = False,
             include_turn_arounds = True)
    pipe.Add(mapmaker.mapmakerutils.FillSimTimestreams,
             sim_map_id = 'sim_maps',
             out_ts_key = 'sim_ts',
             T = sim_map,
             interp_sim = False)
    pipe.Add(delete_mapping_cache)
    pipe.Add(ds)
    pipe.Add(fancy_downsample, ds_factor = ds_factor, debug = True, 
             key_suffix = '_fancy', keys = ('sim_ts',),
             decimate_keys = None)
    pipe.Add(core.G3Writer, 
             filename = os.path.join(savedir, 'sim_frames.g3'))
    pipe.Add(lambda fr: not ('Turnaround' in fr and fr['Turnaround']))
    # pipe.Add(core.Dump, added_message = 'Before mapmaking')
    # pipe.Add(periodic_dump)
    pipe.Add(mapmaker.mapmakerutils.MakeMap, map_in = map_params,
             map_id = 'full_res_sim', ts_in_key = 'sim_ts',
             do_weight = True, make_polarized = False,
             boresight_x_key = bs_ra,
             boresight_y_key = bs_dec)
    pipe.Add(delete_mapping_cache)
    pipe.Add(mapmaker.mapmakerutils.MakeMap, map_in = map_params,
             map_id = 'ds_sim', ts_in_key = 'sim_ts' + suffix,
             do_weight = True, make_polarized = False,
             boresight_x_key = bs_ra + suffix,
             boresight_y_key = bs_dec + suffix)
    pipe.Add(delete_mapping_cache)
    pipe.Add(mapmaker.mapmakerutils.MakeMap, map_in = map_params,
             map_id = 'fancy_sim', ts_in_key = 'sim_ts' + '_fancy',
             do_weight = True, make_polarized = False,
             boresight_x_key = bs_ra + suffix,
             boresight_y_key = bs_dec + suffix)
    pipe.Add(map_extractor)

    pipe.Run(profile = True)

    # Save the maps separately
    fullres = map_extractor.maps['full_res_sim']
    mapmaker.mapmakerutils.save_spt3g_map(os.path.join(savedir, 'full_res.fits'),
                                          fullres['T'], W = fullres['Wunpol'],
                                          overwrite = True)
    basic = map_extractor.maps['ds_sim']
    mapmaker.mapmakerutils.save_spt3g_map(os.path.join(savedir, 'basic.fits'),
                                          basic['T'], W = basic['Wunpol'],
                                          overwrite = True)
    fancy = map_extractor.maps['fancy_sim']
    mapmaker.mapmakerutils.save_spt3g_map(os.path.join(savedir, 'fancy.fits'),
                                          fancy['T'], W = fancy['Wunpol'],
                                          overwrite = True)

if __name__ == '__main__':
    dosim()
