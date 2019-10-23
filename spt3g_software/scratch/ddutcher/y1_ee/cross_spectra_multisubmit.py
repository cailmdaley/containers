"""
cross_spectra_multisubmit.py

This is a command-line script for computing many, many cross-spectra.
Given a directory containing N maps, it will compute all
N-choose-2 cross-spectra.

For help, run
    python cross_spectra_multisubmit.py -h
"""
import os, sys
import argparse as ap
import numpy as np
from glob import glob
from spt3g import core
from spt3g.mapspectra import map_analysis
from spt3g.util import files
from spt3g.cluster.condor_tools import condor_submit


def _get_all_map_pairs(map_dir):
    """
    Given a map directory containing N maps, this function will produce
    a list of pairs of absolute paths for all N-choose-2 map pairs.
    """
    maps = sorted(glob(os.path.join(map_dir, '*.g3*')))
    if len(maps) == 0:
        raise FileNotFoundError("No g3 files in %s"%map_dir)
    map_pairs = []
    for i in maps:
        for j in maps:
            if i != j and (j,i) not in map_pairs:
                map_pairs.append( ( os.path.abspath(i),
                                    os.path.abspath(j) ) )
    return map_pairs


def _create_submission_files(remote_script,
                             map_pairs=None,
                             overwrite=False,
                             remote_args=[],
                             aux_input_files=[],
                             dag_script='job.dag',
                             log_dir='',
                             output_dir='',
                             grid_proxy=None,
                             create_only=False,
                             request_disk=4*core.G3Units.GB,
                             request_memory=4*core.G3Units.GB):
    """
    This function creates the executables, submit files, and dag file
    necessary for submitting jobs to the Open Science Grid via condor.

    `map_pairs` is a list of 2-tuples of absolute paths to map .g3(.gz) files,
    as produced by _get_all_map_pairs()
    """
    assert map_pairs is not None
    
    if not os.path.isfile(remote_script):
        raise FileNotFoundError(remote_script)
    
    if not isinstance(aux_input_files, list):
        aux_input_files = [aux_input_files]
    if not isinstance(remote_args, list):
        remote_args = [remote_args]

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    dag_script = os.path.join(log_dir, dag_script)
    f = open(dag_script,'w')
    for i,j in map_pairs:
        bun_a = os.path.basename(i).partition('.')[0] + '_x_'
        bun_b = os.path.basename(j).partition('.')[0]
        root = bun_a + bun_b
        if not overwrite:
            if os.path.isfile(os.path.join(output_dir, root + '_xspectra.pkl')):
                continue
        in_files = [os.path.basename(i), os.path.basename(j)]
        args = '{infiles} -o {output} {extra}'.format(
            infiles=' '.join(in_files),
            output=root + '_xspectra.pkl',
            extra=' '.join(remote_args))
        jobname = root + '_xspectra'

        condor_submit(remote_script,
                      create_only=True,
                      args=[args],
                      log_root=log_dir,
                      output_root=output_dir,
                      retry=False,
                      jobname=jobname,
                      input_files=[i, j],
                      aux_input_files=aux_input_files,
                      output_files=[root + '_xspectra.pkl'],
                      request_disk=request_disk,
                      request_memory=request_memory,
                      grid_proxy=grid_proxy,
                      verbose=False,
                      globus_uri = \
                      'gsiftp://ceph-gridftp1.grid.uchicago.edu:2811/cephfs')

        f.write('JOB %s %s\n' % (jobname, os.path.join(log_dir, jobname + '.submit')))
        f.write('RETRY %s 5\n' % jobname)

    f.close()
    if not create_only:
        os.system('condor_submit_dag -f -maxidle 100 %s' % dag_script)

# =============================================================================

if __name__ == "__main__":
    # Parent parser for input/output data locations
    # These are used in 'local' and 'submit' modes
    data_args = ap.ArgumentParser(add_help=False)
    data_args.add_argument('--script',
                           default='/home/ddutcher/code/spt3g_software/scratch/ddutcher/y1_ee/y1_ee_cross_spectra.py',
                           help='Script used on maps to compute cross spectra.')
    data_args.add_argument('--map-dir', required=True,
                           help='The directory containing all the maps to cross.')
    data_args.add_argument('--output-dir', default = None,
                           help='The local directory for storing outputs.')
    data_args.add_argument('--mask-file',
                           default='/spt/user/ddutcher/masks/3band_res2_bundle_total_mask.pkl',
                           help='Location of apodization mask pkl file')
    data_args.add_argument('--overwrite', default = False, action='store_true',
                           help='Recompute data already on disk.')

    # Main argument parser
    P0 = ap.ArgumentParser(
        description='Compute the cross-spectra of many maps.',
        formatter_class=ap.ArgumentDefaultsHelpFormatter)
    S = P0.add_subparsers(dest='mode', metavar='MODE', title='subcommands',
                          help='Mode in which to run script. For help, call'
                          ' python %(prog)s %(metavar)s -h')

    # Subparser for calculating all cross-spectra locally.
    # This mode is not recommended for large numbers of cross-spectra,
    # but is useful for testing.
    P1 = S.add_parser('local', parents=[data_args],
                      help='Compute all cross-spectra locally.')

    # Subparser for creating and submitting condor jobs to calculate
    # cross-spectra on the grid
    P2 = S.add_parser('submit', parents=[data_args],
                      help = 'Write executables and submit condor '
                      'dag job to OSG.')
    P2.add_argument('--request-memory', default=4*core.G3Units.GB,
                    help='RAM requirement', type=float)
    P2.add_argument('--request-disk', default=4*core.G3Units.GB,
                    help='Disk space requirement', type=float)
    P2.add_argument('--proxy', default='/home/ddutcher/ddutcher_proxy',
                    help="Location of user's grid proxy cert")
    P2.add_argument('--create-only', default=False, action='store_true',
                    help='Create submission files but do not submit them')
    P2.add_argument('--log-dir', default='.', type=str,
                    help='Location to store condor logging outputs')
    P2.add_argument('--dag-script', default='out.dag', type=str,
                    help='Name of dag script to be submitted')

    args, extra_args = P0.parse_known_args()

    if args.mode == 'submit':
        for arg, val in args.__dict__.items():
            if val == 'None':
                args.__dict__[arg] = None
        if args.mask_file is not None:
            extra_args += ['--apod %s' % os.path.basename(args.mask_file)]

        map_pairs = _get_all_map_pairs(args.map_dir)
        _create_submission_files(
            args.script,
            map_pairs=map_pairs,
            aux_input_files=args.mask_file,
            dag_script=args.dag_script,
            log_dir = args.log_dir,
            output_dir=args.output_dir,
            grid_proxy=args.proxy,
            create_only=args.create_only,
            request_disk=args.request_disk,
            request_memory=args.request_memory,
            remote_args=extra_args,
            overwrite=args.overwrite)

    if args.mode == 'local':
        for arg,val in args.__dict__.items():
            if val == 'None':
                args.__dict__[arg] = None
        if args.mask_file is not None:
            extra_args += ['--apod %s' % args.mask_file]
        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)

        map_pairs = _get_all_map_pairs(args.map_dir)
        for map1, map2 in map_pairs:
            filename = (os.path.basename(map1).partition('.')[0]
                        + '_x_'
                        + os.path.basename(map2).partition('.')[0]
                        + '_xspectra.pkl')
            if not args.overwrite:
                if os.path.isfile(os.path.join(args.output_dir, filename)):
                    continue
            print('Making '+filename)
            user_cmd = 'python {script} {infiles} -o {output} {extra}'.format(
                script=args.script,
                infiles=' '.join([map1, map2]),
                output=os.path.join(args.output_dir, filename),
                extra=' '.join(extra_args))
            
#             print('Running user command: '+user_cmd)
            os.system(user_cmd)


# python cross_spectra_multisubmit.py submit --map-dir /spt/user/ddutcher/null_maps/data/150GHz/azimuth --output-dir /spt/user/ddutcher/null_xspectra/data/150GHz/azimuth --log-dir /scratch/ddutcher/condor_logs/xspectra/data/150GHz/azimuth --dag-script /scratch/ddutcher/condor_logs/xspectra/data/150GHz/azimuth.dag --create-only

# python cross_spectra_multisubmit.py submit --map-dir /spt/user/ddutcher/sim_bundles/150GHz/simmap0000/total --output-dir /spt/user/ddutcher/xspectra/sim/simmap0000/150GHz --log-dir /scratch/ddutcher/condor_logs/xspectra/sim/150GHz --dag-script /scratch/ddutcher/condor_logs/xspectra/sim/150GHz/simmap0000.dag

# python cross_spectra_multisubmit.py submit --map-dir /spt/user/ddutcher/null_maps/sim/simmap0001/150GHz/wafer --output-dir /spt/user/ddutcher/null_xspectra/sim/simmap0001/150GHz/wafer --log-dir /scratch/ddutcher/condor_logs/xspectra/sim/150GHz/wafer --dag-script /scratch/ddutcher/condor_logs/xspectra/sim/150GHz/simmap0001_wafer.dag