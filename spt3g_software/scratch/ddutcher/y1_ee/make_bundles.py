"""
Functions for combining 1500d subfield observations
and mock-observations into bundles.
"""
import os
import gc
import argparse
import numpy as np
from glob import glob
from spt3g import core, calibration, coordinateutils, util
from spt3g.util import files

spt3g_software = os.environ['SPT3G_SOFTWARE_PATH']
map_dir = '/spt/user/ddutcher/'


def define_bundles(num_bundles=30,
                   output_dir='/spt/user/ddutcher/bundles',
                   map_stats='/home/ddutcher/data/1500d_map_stats.pkl',
                   return_data=False,
                   save=True):
    """
    Create coadds of observations ("bundles") that cover the full
    1500d field and have approximately equal weights

    Parameters
    ---------
    num_bundles : int
        The number of bundles to define
    output_dir : str
        Where to save the .pkl of bundle definitions
    map_stats : str
        The location of a pickle file containing the
        outputs of map_stats.py (or get_map_stats.py) for each observation,
        grouped by field, e.g.{'ra0hdec-44.75':{'47429579':{'med_w':,...}}}
        This contains the median weight for each observation.
    save : bool
        Save the resulting bundle_obsids.pkl to output_dir
    return_data : bool
        Return bundle_obsids.pkl

    Returns
    -------
    Dictionary mapping bundle index to list of ObsIDs
    """
    # Strategy is to coadd with each subfield chronologically
    # until the desired average weight is achieved, then coadd
    # the subfield coadds.
    
    bad_obs = [48085553, 53588441, 57262064, 57282273, 57248970]

    map_stats = files.load_pickle(map_stats)

    med_w = dict()
    for field in map_stats.keys():
        med_w[field]={'obsid':[], 'w':[]}
        for obsid, d in map_stats[field].items():
            if int(obsid) in bad_obs:
                continue
            if np.isnan(d['med_w']):
                continue
            med_w[field]['obsid'].append(obsid)
            med_w[field]['w'].append(d['med_w'])

    bundle_dict= dict()
    for field, d in med_w.items():
        # keeping list of sums like this is useful for checking
        # the uniformity of the subfield bundle weights
        bsums = np.zeros(num_bundles)
        bidx = 0
        tot = sum(d['w'])
        # In case the obsids are not in numerical (chronological) order:
        ind = np.argsort(d['obsid'])
        for i in ind:
            curr_diff = np.abs(bsums[bidx] - tot/num_bundles)
            # figure out if new obs helps, or if we were better before
            if np.abs(d['w'][i]+bsums[bidx] - tot/num_bundles) > curr_diff:
                if bidx == len(bsums)-1:
                    # don't create extra lame bundle
                    pass
                else:
                    bidx += 1
            bsums[bidx] += d['w'][i]

            obsids = bundle_dict.get(bidx, [])
            bundle_dict[bidx] = obsids + [d['obsid'][i]]
    if save:
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        files.save_pickle(
            bundle_dict, os.path.join(output_dir, 'bundle_obsids.pkl'))
    if return_data:
        return bundle_dict


def define_az_bundles(num_bundles=30,
                      output_dir='/spt/user/ddutcher/az_bundles',
                      map_stats='/home/ddutcher/data/1500d_map_stats.pkl',
                      config='/spt/user/ddutcher/null_configs/azimuth.pkl',
                      return_data=False,
                      save=True):
    """
    Create coadds of observations ("bundles") that cover the full
    1500d field and have approximately equal weights.

    Paramters
    ---------
    num_bundles : int
        The number of bundles to define
    output_dir : str
        Where to save the .pkl of bundle definitions
    map_stats : str
        The location of a pickle file containing the
        outputs of map_stats.py (or get_map_stats.py) for each observation,
        grouped by field, e.g.{'ra0hdec-44.75':{'47429579':{'med_w':,...}}}
        This contains the median weight for each observation.
    save : bool
        Save the resulting az_bundle_obsids.pkl to output_dir
    return_data : bool
        Return az_bundle_obsids.pkl

    Returns
    -------
    Dictionary mapping bundle index to list of ObsIDs
    """
    # Strategy:
    # 1. Define target weight per subfield based on num_bundles
    # 2. Sort subfield observations by azimuthal distance from 153
    # 3. Coadd observations moving out from 153 until we hit target weight
    # 4. Coadd subfield coadds

    bad_obs = [48085553, 53588441, 57262064, 57282273, 57248970]

    map_stats = files.load_pickle(map_stats)
    med_w = dict()
    for field in map_stats.keys():
        med_w[field]={'obsid':[], 'w':[]}
        for obsid, d in map_stats[field].items():
            if int(obsid) in bad_obs:
                continue
            if np.isnan(d['med_w']):
                continue
            med_w[field]['obsid'].append(obsid)
            med_w[field]['w'].append(d['med_w'])

    config = files.load_pickle(config)

    bundle_dict= dict()
    for field, d in med_w.items():
        # keeping list of sums like this is useful for checking
        # the uniformity of the subfield bundle weights
        bsums = np.zeros(num_bundles)
        bidx = 0
        tot = sum(d['w'])

        # Sort the obs by azimuthal distance from 153
        dists = []
        for obsid in d['obsid']:
            az = config[int(obsid)]
            dist = np.abs(153 - az)
            if dist > 180:
                dist = 180 - (dist % 180)
            dists.append(dist)
        ind = np.argsort(dists)

        for i in ind:
            curr_diff = np.abs(bsums[bidx] - tot/num_bundles)
            # figure out if new obs helps, or if we were better before
            if np.abs(d['w'][i]+bsums[bidx] - tot/num_bundles) > curr_diff:
                if bidx == len(bsums)-1:
                    # don't create extra lame bundle
                    pass
                else:
                    bidx += 1
            bsums[bidx] += d['w'][i]
            obsids = bundle_dict.get(bidx, [])
            bundle_dict[bidx] = obsids + [d['obsid'][i]]
    if save:
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        files.save_pickle(
            bundle_dict, os.path.join(output_dir, 'az_bundle_obsids.pkl'))
    if return_data:
        return bundle_dict


def define_sat_bundles(num_bundles=30,
                       output_dir='/spt/user/ddutcher/sat_bundles',
                       map_stats='/home/ddutcher/data/1500d_map_stats.pkl',
                       config='/spt/user/ddutcher/null_configs/saturation.pkl',
                       return_data=False,
                       save=True):
    """
    Create coadds of observations ("bundles") that cover the full
    1500d field and have approximately equal weights.

    Parameters
    ----------
    num_bundles : int
        The number of bundles to define
    output_dir : str
        Where to save the .pkl of bundle definitions
    map_stats : str
        The location of a pickle file containing the
        outputs of map_stats.py (or get_map_stats.py) for each observation,
        grouped by field, e.g.{'ra0hdec-44.75':{'47429579':{'med_w':,...}}}
        This contains the median weight for each observation.
    save : bool
        Save the resulting sat_bundle_obsids.pkl to output_dir
    return_data : bool
        Return sat_bundle_obsids.pkl

    Returns
    -------
    Dictionary mapping bundle index to list of ObsIDs

    Notes
    -----
    The existing bundles have a decent spread of saturation values,
    so I have chosen not to use these sat_bundles.
    """
    # Strategy:
    # 1. Define target weight per subfield based on num_bundles
    # 2. Sort subfield observations by avg. number of saturated bolos per scan.
    # 3. Coadd observations until we hit target weight
    # 4. Coadd subfield coadds

    bad_obs = [48085553, 53588441, 57262064, 57282273, 57248970]

    map_stats = files.load_pickle(map_stats)
    med_w = dict()
    for field in map_stats.keys():
        med_w[field]={'obsid':[], 'w':[]}
        for obsid, d in map_stats[field].items():
            if int(obsid) in bad_obs:
                continue
            if np.isnan(d['med_w']):
                continue
            med_w[field]['obsid'].append(obsid)
            med_w[field]['w'].append(d['med_w'])

    config = files.load_pickle(config)

    bundle_dict= dict()
    for field, d in med_w.items():
        # keeping list of sums like this is useful for checking
        # the uniformity of the subfield bundle weights
        bsums = np.zeros(num_bundles)
        bidx = 0
        tot = sum(d['w'])

        # Sort the obs by saturation
        num_sats = []
        for obsid in d['obsid']:
            num_sats.append(config[int(obsid)])
        ind = np.argsort(num_sats)[::-1]

        for i in ind:
            curr_diff = np.abs(bsums[bidx] - tot/num_bundles)
            # figure out if new obs helps, or if we were better before
            if np.abs(d['w'][i]+bsums[bidx] - tot/num_bundles) > curr_diff:
                if bidx == len(bsums)-1:
                    # don't create extra lame bundle
                    pass
                else:
                    bidx += 1
            bsums[bidx] += d['w'][i]
            obsids = bundle_dict.get(bidx, [])
            bundle_dict[bidx] = obsids + [d['obsid'][i]]
    if save:
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        files.save_pickle(
            bundle_dict, os.path.join(output_dir, 'sat_bundle_obsids.pkl'))
    if return_data:
        return bundle_dict, bsums


class GetInfo(object):
    def __init__(self, info_to_drop= ['SurvivingBolos','DroppedBolos',
                                      'DroppedBoloStats',
                                      'input_files', 'output',
                                      'split_left_right',
                                      'wafers_to_include']):
        self.drop_list = info_to_drop
    def __call__(self, fr):
        if fr.type == core.G3FrameType.PipelineInfo:
            if len(fr.keys()) < 3:
                return False
            for k in self.drop_list:
                fr.pop(k, None)
            return
        if fr.type == core.G3FrameType.Observation:
            print(fr)
        if fr.type == core.G3FrameType.EndProcessing:
            return
        if not ('Id' in fr and 'combined' in fr['Id']):
            return False


def make_bundles(bundle_def='/spt/user/ddutcher/bundles/bundle_obsids.pkl',
                 output_dir='/spt/user/ddutcher/bundles',
                 out_tag=None,
                 map_runs=['high_150GHz_left_maps'],
                 map_ids=['Left150GHz'],
                 compress=False,
                 sim_map_index=None):
    """
    Coadd obervations according to the bundle definitions in bundle_def
    """
    # This is slooow, so I recommending running multiple copies    
    bundle_def = files.load_pickle(bundle_def)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    if not isinstance(map_ids, list):
        map_ids = [map_ids]
        
    if not isinstance(map_runs, list):
        map_runs = [map_runs]

    if out_tag is None:
        out_tag = map_runs[0]

    if sim_map_index is not None:
        sim_map_index = int(sim_map_index)
        sim = 'simmap%04d*' % sim_map_index
    else:
        sim = ''

    for bidx, obsids in bundle_def.items():
        if os.path.isfile(os.path.join(
            output_dir, 'bun%02d_%s.g3*' % (bidx, out_tag))):
            continue
        print(bidx)
        obs_to_add = []
        for obsid in obsids:
            for map_run in map_runs:
                pth = os.path.join(
                    map_dir,
                    'ra0hdec-??.??',
                    'y1_ee_20190811',
                    map_run,
                    obsid,
                    sim+map_run+'*.g3*')
                mp = glob(pth)
                if len(mp) != 1:
                    raise FileNotFoundError(pth)
                obs_to_add += mp

        if compress:
            suffix='.gz'
        else:
            suffix=''

        pipe = core.G3Pipeline()
        pipe.Add(core.G3Reader,filename = obs_to_add)
        for map_id in map_ids:
            pipe.Add(util.framecombiner.MapFrameCombiner, fr_id = map_id)
        pipe.Add(GetInfo)
        pipe.Add(core.DeduplicateMetadata,
                 dataframetype = [core.G3FrameType.Map])
        pipe.Add(core.G3Writer, filename = os.path.join(
            output_dir, 'bun%02d_%s.g3%s' % (bidx, out_tag, suffix)))
        pipe.Run()

        gc.collect()

def make_all_az_bundles(bands=['90GHz', '150GHz', '220GHz'],
                        compress=False, sim_map_index=None):
    '''
    Coadds individual observations into azimuth bundles.
    '''
    if not isinstance(bands, list):
        bands = [bands]
    for band in bands:
        assert band in ['90GHz', '150GHz', '220GHz']
        if sim_map_index is not None:
            sim_map_index=int(sim_map_index)
            output_dir = os.path.join(
                '/spt/user/ddutcher/sim_bundles/simmap%04d' % sim_map_index,
                band,
                'azimuth')
        else:
            output_dir = os.path.join(
                '/spt/user/ddutcher/bundles', band, 'azimuth')
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        make_bundles(bundle_def='/spt/user/ddutcher/bundles/az_bundle_obsids.pkl',
                     output_dir=output_dir,
                     map_runs=['low_'+band+'_left_maps',
                               'high_'+band+'_left_maps',
                               'high_'+band+'_right_maps',
                               'low_'+band+'_right_maps'],
                     map_ids='*'+band,
                     out_tag='azimuth_'+band,
                     compress=compress,
                     sim_map_index=sim_map_index)


def make_all_bundles(bands=['90GHz', '150GHz', '220GHz'],
                     compress=False, sim_map_index=None):
    '''
    Coadds individual observations into bundle components that
    can be then summed into total, left-right, or wafer bundles.
    '''
    if not isinstance(bands, list):
        bands = [bands]
    for band in bands:
        assert band in ['90GHz', '150GHz', '220GHz']
        if sim_map_index is not None:
            sim_map_index=int(sim_map_index)
            output_dir = os.path.join(
                '/spt/user/ddutcher/sim_bundles/simmap%04d' % sim_map_index,
                band)
        else:
            output_dir = '/spt/user/ddutcher/bundles/'+band
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        for direction in ['left', 'right']:
            for resp in ['high', 'low']:
                map_run = '_'.join([resp, band, direction, 'maps'])
                make_bundles(bundle_def='/spt/user/ddutcher/bundles/bundle_obsids.pkl',
                             output_dir=output_dir,
                             map_runs=[map_run],
                             map_ids='*'+band,
                             compress=compress,
                             sim_map_index=sim_map_index)


def _adder(input, output):
    '''
    Simple map-coadding helper function.
    '''
    outdir = os.path.dirname(output)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename = input)
    pipe.Add(util.framecombiner.MapFrameCombiner)
    pipe.Add(lambda frame: frame.type == core.G3FrameType.Map)
    pipe.Add(core.G3Writer, filename=output)
    pipe.Run()
    
    gc.collect()

def sum_bundles(bands=['90GHz', '150GHz', '220GHz'],
                sim_map_index=None):
    '''
    Combines bundle components into high-low, left-right, and total bundles.
    '''
    if not isinstance(bands, list):
        bands = [bands]
    if sim_map_index is not None:
        sim = 'simmap%04d' % int(sim_map_index)
        bun_dir = os.path.join(
            '/spt/user/ddutcher/sim_bundles',
            sim)
    else:
        bun_dir = '/spt/user/ddutcher/bundles'

    for band in bands:
        assert band in ['90GHz', '150GHz', '220GHz']
        bundles = sorted(glob(os.path.join(
            bun_dir, band, '*.g3*')))
        assert len(bundles) == 120
        print('Found %s bundles' % len(bundles))
        bdict = dict()

        for bun in bundles:
            bun_label = os.path.basename(bun).split('_')[0]
            to_add = bdict.get(bun_label, [])
            to_add.append(bun)
            bdict[bun_label] = to_add

        for bun_label, blist in bdict.items():
            print(bun_label, blist)
            high = []
            low = []
            left = []
            right = []
            for pth in blist:
                if 'high' in pth: #combines left,right 'high's
                    high.append(pth)
                if 'low' in pth: # combines left, right 'low's
                    low.append(pth)
                if 'left' in pth: # combines high and low 'left's
                    left.append(pth)
                if 'right' in pth: # combines high and low 'right's
                    right.append(pth)                    

            _adder(high, os.path.join(
                bun_dir, '%s/high_low/%s_high_%s.g3.gz' % (band, bun_label, band)))

            _adder(low, os.path.join(
                bun_dir, '%s/high_low/%s_low_%s.g3.gz' % (band, bun_label, band)))

            _adder(left, os.path.join(
                bun_dir, '%s/left_right/%s_left_%s.g3.gz' % (band, bun_label, band)))

            _adder(right, os.path.join(
                bun_dir, '%s/left_right/%s_right_%s.g3.gz' % (band, bun_label, band)))

            sum1 = os.path.join(
                bun_dir, '%s/left_right/%s_right_%s.g3.gz' % (band, bun_label, band))
            sum2 = os.path.join(
                bun_dir, '%s/left_right/%s_left_%s.g3.gz' % (band, bun_label, band))

            _adder([sum1, sum2], os.path.join(
                bun_dir, '%s/total/%s_%s.g3.gz' % (band, bun_label, band)))

            

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["make", "az", "sum"])
    parser.add_argument("--sim-map-index", type=int, default=None)
    parser.add_argument("--num-maps", type=int, default=6828)
    args = parser.parse_args()

    if args.sim_map_index is not None:
        filestr = 'simmap%04d*.g3*' % args.sim_map_index
    else:
        filestr = '[hl]*.g3*'

    if args.mode != "sum":
        test = glob(
            '/spt/user/ddutcher/ra0hdec-??.??/y1_ee_20190811/*/*/' + filestr)
        assert len(test) == args.num_maps
        if args.mode=="az":
            make_all_az_bundles(sim_map_index=args.sim_map_index, compress=True)
        else:
            make_all_bundles(sim_map_index= args.sim_map_index, compress=True)
    else:
        sum_bundles(sim_map_index=args.sim_map_index)
