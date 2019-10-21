"""
Functions for combining Planck mock-observations into bundles.
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
                 field='ra0hdec-??.??',
                 map_runs=['high_150GHz_left_maps'],
                 map_ids=['Left150GHz'],
                 compress=False,
                 sim_prefix=None):
    """
    Coadd obervations according to the bundle definitions in bundle_def
    
    Note: This splits coadds up by subfield
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

    if sim_prefix is not None:
        if 'planck' in str(sim_prefix):
            sim = sim_prefix
        else:
            sim = 'simmap%04d' % int(sim_map_index)
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
                    field,
                    'y1_ee_20190811',
                    map_run,
                    obsid,
                    sim+'*'+map_run+'*.g3*')
                mp = glob(pth)
                if len(mp) != 1:
                    continue
#                     raise FileNotFoundError(pth)
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
        pipe.Add(core.G3Writer,
                 filename = os.path.join(
                     output_dir,
                     'bun%02d_%s.g3%s' % (bidx, out_tag, suffix)))
        pipe.Run()
        
        gc.collect()


def make_all_az_bundles(bands=['90GHz', '150GHz', '220GHz'],
                        compress=False, sim_prefix=None):
    '''
    Coadds individual observations into azimuth bundles.
    '''
    if not isinstance(bands, list):
        bands = [bands]
    for band in bands:
        assert band in ['90GHz', '150GHz', '220GHz']
        if sim_prefix is not None:
            if 'planck' in str(sim_prefix):
                sim = sim_prefix
            else:
                sim = 'simmap%04d' % int(sim_prefix)
            output_dir = os.path.join(
                '/spt/user/ddutcher/sim_bundles', sim, band, 'azimuth')
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
                     sim_prefix=sim_prefix)


def make_all_bundles(bands=['90GHz', '150GHz', '220GHz'],
                     field='ra0hdec-??.??',
                     compress=False,
                     sim_prefix=None):
    '''
    Coadds individual observations into bundle components that
    can be then summed into total, left-right, or wafer bundles.
    '''
    if not isinstance(bands, list):
        bands = [bands]
    for band in bands:
        assert band in ['90GHz', '150GHz', '220GHz']
        if sim_prefix is not None:
            if 'planck' in str(sim_prefix):
                sim = sim_prefix
            else:
                sim = 'simmap%04d' % int(sim_prefix)
            if field != 'ra0hdec-??.??':
                output_dir = os.path.join(
                    '/spt/user/ddutcher/sim_bundles', sim, band, field)
            else:
                output_dir = os.path.join(
                    '/spt/user/ddutcher/sim_bundles', sim, band)
        else:
            output_dir = '/spt/user/ddutcher/bundles/'+band
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        for direction in ['left', 'right']:
            for resp in ['high', 'low']:
                map_run = '_'.join([resp, band, direction, 'maps'])
                make_bundles(bundle_def='/spt/user/ddutcher/bundles/bundle_obsids.pkl',
                             output_dir=output_dir,
                             field=field,
                             map_runs=[map_run],
                             map_ids='*'+band,
                             compress=compress,
                             sim_prefix=sim_prefix)


def _adder(input, output):
    outdir = os.path.dirname(output)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader,filename = input)
    pipe.Add(util.framecombiner.MapFrameCombiner)
    pipe.Add(lambda frame: frame.type == core.G3FrameType.Map)
    pipe.Add(core.G3Writer, filename=output)
    pipe.Run()
    
    gc.collect()


def sum_bundles(bands=['90GHz', '150GHz', '220GHz'],
#                 fields=['ra0hdec-44.75', 'ra0hdec-52.25'],
#                         'ra0hdec-59.75', 'ra0hdec-67.25'],
                sim_prefix=None):
    '''
    Combines bundle components into total bundles.
    '''
    if not isinstance(bands, list):
        bands = [bands]
    if sim_prefix is not None:
        if 'planck' in str(sim_prefix):
            sim = sim_prefix
        else:
            sim = 'simmap%04d' % int(sim_prefix)
        bun_dir = os.path.join(
            '/spt/user/ddutcher/sim_bundles', sim)
    else:
        bun_dir = '/spt/user/ddutcher/bundles'

    for band in bands:
        assert band in ['90GHz', '150GHz', '220GHz']
        bundles = sorted(
            glob(os.path.join(bun_dir, band, 'high_el', '*.g3*'))
            + glob(os.path.join(bun_dir, band, 'low_el', '*.g3*')))
        assert len(bundles) == 60
        print('Found %s bundles' % len(bundles))
        bdict = dict()

        for bun in bundles:
            bun_label = os.path.basename(bun).split('_')[0]
            to_add = bdict.get(bun_label, [])
            to_add.append(bun)
            bdict[bun_label] = to_add

        for bun_label, blist in bdict.items():
            print(bun_label, blist)

            _adder(blist, os.path.join(
                bun_dir, '%s/full/%s_%s.g3.gz' % (band, bun_label, band)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["make", "az", "sum"])
    parser.add_argument("sim_prefix")
    parser.add_argument("--field", default='ra0hdec-??.??', type=str)
    args = parser.parse_args()

    if args.mode != "sum":
        if args.sim_prefix is not None:
            if 'planck' in str(args.sim_prefix):
                sim = args.sim_prefix
            else:
                sim = 'simmap%04d' % int(args.sim_prefix)
        test = glob(
            '/spt/user/ddutcher/ra0hdec-??.??/y1_ee_20190811/*/*/%s*' % args.sim_prefix)
        assert len(test) == 6828
        if args.mode=="az":
            make_all_az_bundles(sim_prefix=args.sim_prefix, compress=True)
        else:
            make_all_bundles(sim_prefix= args.sim_prefix, compress=True, field=args.field)
    else:
        sum_bundles(sim_prefix=args.sim_prefix)
