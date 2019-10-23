import os
from glob import glob
from spt3g.cluster.condor_tools import condor_submit
from spt3g import core
import argparse as ap
import time

user = os.getenv('LOGNAME')
default_batch_name = 'condor_{}'.format(time.strftime('%Y%m%d_%H%M%S'))
default_log_root = os.path.join('/scratch', user)
default_output_root = os.path.join('/spt/user', user)

P = ap.ArgumentParser()
P.add_argument('script', help='Script to run.')
P.add_argument('obstype', help='Observation type.')
P.add_argument('obsids', nargs='+',
               help='Globbable list of observation ids of type obstype to '
               'process.')
P.add_argument('--args', nargs='+', default=[],
               help='Arguments passed to the user script')
P.add_argument('--submit-jobs', action='store_true',
               help='Submit the jobs to the grid. When not set, submission '
               'files are created but jobs are not actually submitted. This '
               'is useful for debugging.')
P.add_argument('--batch-name', default=default_batch_name,
               help='Name of the submit job.  This is used to create the log '
               'and output directories, and to tag the input software.')
P.add_argument('--raw-data-root', default='/spt/data/bolodata/fullrate/',
               help='Path to raw bolometer data.')
P.add_argument('--log-root', default=default_log_root,
               help='Local directory where the log files are stored. This '
               'should be a location on /scratch, since these are typically '
               'lots of small files.')
P.add_argument('--output-root', default=default_output_root,
               help='Local directory where the output files are stored. This '
               'should be a location on cephfs (/spt), since output files '
               'must be transfered using GridFTP')
P.add_argument('--request-cpus', type=int, default=1,
               help='Number of CPUs to request for the job')
P.add_argument('--request-memory', type=int, default=2*core.G3Units.GB,
               help='Memory required, in GB')
P.add_argument('--request-disk', type=int, default=4*core.G3Units.GB,
               help='Disk space required, in GB')
P.add_argument('--clustertools-version', default=None,
               choices=['py2-v1', 'py3-v1', 'py3-v2'],
               help='CVMFS environment to use. If not supplied, check user '
               'environment.')

args = P.parse_args()

# extract the paths to the observation data to process
obsid_paths = []
for obsid_pattern in args.obsids:
    obsid_paths.extend([path for path in glob(os.path.join(args.raw_data_root,
                                                           args.obstype,
                                                           obsid_pattern))])

# submit the jobs, one per observation
for obs_path in obsid_paths:
    print(obs_path)
    obsid = os.path.basename(obs_path)
    infiles = glob(os.path.join(obs_path, '*'))
    job_name = '{}_{}_{}'.format(args.batch_name,
                                 args.obstype,
                                 obsid)

    print('Processing observation {}'.format(obsid))
    cluster, f_submit, f_script = condor_submit(args.script,
                                                create_only = (not args.submit_jobs),
                                                args = args.args,
                                                log_root = args.log_root, 
                                                output_root = args.output_root,
                                                jobname = job_name,
                                                grid_proxy = '/home/adama/.globus/grid_proxy',
                                                input_files = infiles,
                                                output_files = '{}.g3'.format(job_name),
                                                request_cpus = args.request_cpus,
                                                request_disk = args.request_disk,
                                                request_memory = args.request_memory,
                                                clustertools_version = args.clustertools_version)
