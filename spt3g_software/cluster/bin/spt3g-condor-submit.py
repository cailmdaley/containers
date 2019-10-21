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
default_grid_proxy = os.getenv('X509_USER_PROXY')

P = ap.ArgumentParser(description='Submit jobs to grid from amundsen/scott.',
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
# global condor settings that could never depend on choice of workflow
P.add_argument('--caller', default='python',
               help='Program that calls the user script')
P.add_argument('--submit-jobs', action='store_true',
               help='Submit the jobs to the grid. When not set, submission '
               'files are created but jobs are not actually submitted. This '
               'is useful for debugging.')
P.add_argument('--batch-name', default=default_batch_name,
               help='Name of the submit job.  This is used to create the log '
               'and output directories, and to tag the input software.')
P.add_argument('--grid-proxy', default=default_grid_proxy,
               help='Path to a valid grid proxy file.  Required if any '
               'input or output files are supplied.')
P.add_argument('--no-retry', dest='retry', action='store_false', default=True,
               help='If supplied, do not mark failed jobs as held to periodically retry.')
P.add_argument('--requirements', default=None,
               help='Computing requirements specification')
P.add_argument('--extra-requirements', default=None,
               help='Additional computing requirements to add to the defaults.')
P.add_argument('--bigmem', default=False, action='store_true',
               help='Add the SPT_BIGMEM requirement to the requirements list.')
P.add_argument('--request-cpus', type=int, default=1,
               help='Number of CPUs to request for the job')
P.add_argument('--request-memory', type=int, default=2*core.G3Units.GB,
               help='Memory required, in bytes')
P.add_argument('--request-disk', type=int, default=4*core.G3Units.GB,
               help='Disk space required, in bytes')
P.add_argument('--when-to-transfer-output', default='ON_EXIT',
               help='When the `aux_output_files` are to be transfered '
               'back to the local `log_root`.')
P.add_argument('--globus-uri', default='gsiftp://gridftp.grid.uchicago.edu:2811/cephfs',
               help='URI prefix to use for accessing files on the Ceph pool')
P.add_argument('--clustertools-version', default=None, choices=['py2-v1', 'py3-v1', 'py3-v2'],
               help='CVMFS environment to use. If not supplied, check user environment.')
subparsers = P.add_subparsers(title='workflow',
                              description='Type of processing workflow to use '
                              'for executing script on the grid.')

P_obs = subparsers.add_parser('process_obs',
                              help='Workflow in which a script is run on '
                              'many observations, one per job submitted '
                              'to the grid. Useful for processing many '
                              'observations in a standardized way that '
                              'depends only on each observation and its '
                              'calibration data.',
                              formatter_class=ap.ArgumentDefaultsHelpFormatter)
# Arguments related to the workflow for processing observations. Many of
# these are duplicated in the `simple` workflow, but we define them here
# because they way that they are handled could, in principle, be different for
# the two workflows (even though they typically are not).
P_obs.add_argument('script', help='Script to run.')
P_obs.add_argument('obstype', help='Observation type.')
P_obs.add_argument('obsids', nargs='+',
                   help='Globbable list of observation ids of type obstype to '
                   'process.')
P_obs.add_argument('--args', nargs='+', default=[],
               help='Arguments passed to the user script')
P_obs.add_argument('--raw-data-root', default='/spt/data/bolodata/fullrate/',
               help='Path to raw bolometer data.')
P_obs.add_argument('--log-root', default=default_log_root,
               help='Local directory where the log files are stored. This '
               'should be a location on /scratch, since these are typically '
               'lots of small files.')
P_obs.add_argument('--output-root', default=default_output_root,
               help='Local directory where the output files are stored. This '
               'should be a location on cephfs (/spt), since output files '
               'must be transfered using GridFTP')
P_obs.add_argument('--aux-input-files', nargs='+', default=[],
               help='Small files to transfer with the submit script.  '
               'The command script and grid proxy are automatically '
               'added to this list.')
P_obs.add_argument('--aux-output-files', nargs='+', default=[],
               help='Small files to transfer to `log_root` on job completion.  '
               'Paths must be relative to the remote working directory.')
P_obs.add_argument('--user-code', default='',
               help='User code to add to the run script, before the main script is called, '
               'but after the environment is configured and input files are transfered.')
P_obs.add_argument('--user-code-post', default='',
               help='User code to add to the run script, after the main script is called, '
               'but before output files are transfered.')

P_simple = subparsers.add_parser('simple',
                                 help='Workflow in which a script is run on '
                                 'the same arbitrary, user-supplied data, '
                                 'which is copied to every grid node. '
                                 'Possibly useful for simulations in which '
                                 'experimental data is not required for '
                                 'each job.',
                                 formatter_class=ap.ArgumentDefaultsHelpFormatter)
# Arguments related to the `simple` workflow. See note about repeated arguments
# in the comment near the definition of the `process_obs` mode.
P_simple.add_argument('script', help='Script to run.')
P_simple.add_argument('--args', nargs='+', default=[],
                      help='Arguments passed to the user script')
P_simple.add_argument('--log-root', default=default_log_root,
                      help='Local directory where the log files are stored. This '
                      'should be a location on /scratch, since these are typically '
                      'lots of small files.')
P_simple.add_argument('--output-root', default=default_output_root,
                      help='Local directory where the output files are stored. This '
                      'should be a location on cephfs (/spt), since output files '
                      'must be transfered using GridFTP')
P_simple.add_argument('--aux-input-files', nargs='+', default=[],
                      help='Small files to transfer with the submit script.  '
                      'The command script and grid proxy are automatically '
                      'added to this list.')
P_simple.add_argument('--aux-output-files', nargs='+', default=[],
                      help='Small files to transfer to `log_root` on job completion.  '
                      'Paths must be relative to the remote working directory.')
P_simple.add_argument('--user-code', default='',
                      help='User code to add to the run script, before the main script is called, '
                      'but after the environment is configured and input files are transfered.')
P_simple.add_argument('--user-code-post', default='',
                      help='User code to add to the run script, after the main script is called, '
                      'but before output files are transfered.')
P_simple.add_argument('--queue', default='1',
                      help='The job queue specification, either an integer '
                      'number of job copies or a specification string -- see '
                      'Condor documentation.')

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

    # organize log files by observation type and ID in file tree
    obs_log_dirname = os.path.join(args.log_root, args.batch_name,
                                   args.obstype, obsid)
    os.makedirs(obs_log_dirname)

    # organize output files by observation type and ID in file tree
    obs_out_dirname = os.path.join(args.output_root, args.batch_name,
                                   args.obstype, obsid)
    os.makedirs(obs_out_dirname)

    print('Processing observation {}'.format(obsid))
    cluster, f_submit, f_script = condor_submit(args.script,
                                                create_only = (not args.submit_jobs),
                                                args = args.args,
                                                log_root = obs_log_dirname, 
                                                output_root = obs_out_dirname,
                                                jobname = job_name,
                                                grid_proxy = '/home/adama/.globus/grid_proxy',
                                                input_files = infiles,
                                                output_files = '{}.g3'.format(job_name),
                                                request_cpus = args.request_cpus,
                                                request_disk = args.request_disk,
                                                request_memory = args.request_memory,
                                                clustertools_version = args.clustertools_version)
