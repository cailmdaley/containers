'''
submit_ps_maps.py

This script submits map jobs for each of the N brightest point sources
using a file that has been pre-generated to list the observations containing
the point source of interest.
'''
import os
import numpy as np
from glob import glob
import argparse
import getpass
from spt3g import core
from spt3g.cluster.condor_tools import condor_submit

parser = argparse.ArgumentParser()
parser.add_argument('--submit', action='store_true', default=False,
                    help='Submit jobs. Otherwise, just create files.')
parser.add_argument('--log-dir', action='store',
                    default=os.path.join('/scratch', getpass.getuser(),
                                         'condor_logs'),
                    help='Absolute path to condor log directory')
parser.add_argument('--retries', action='store', type=int, default=0,
                    help='Number of times to retry a job before giving up')
parser.add_argument('--grid-proxy', action='store',
                    default=os.path.join(os.getenv('HOME'),
                                         getpass.getuser()+'_proxy'),
                    help='Path to grid proxy')
parser.add_argument('--config-file', action='store', default='ps_maps.yaml',
                    help='The name of the yaml config file')
parser.add_argument('--files-per-job', action='store', type=int, default=20,
                    help='Parallelize by making a map per this many files to'
                    ' later coadd. If 0, do all files in one job.')
parser.add_argument('-o', '--output-path', action='store',
                    default=os.path.join('/spt', 'user', getpass.getuser(),
                                         'ps_maps'),
                    help='Top level subdirectory name for map run')
args = parser.parse_args()

# mapmaking yaml file
config = os.path.abspath(args.config_file)

# Location of script to be run
script = os.path.join(os.getenv('SPT3G_SOFTWARE_PATH'), 'std_processing',
                      'mapmakers', 'master_field_mapmaker.py')

# g3files for each point source
ps_list = os.path.join('/spt', 'user', 'agambrel', 'configs',
                       'point_source_obs_{}.npz')
# list of source locs
pts_mask = os.path.join(os.getenv('SPT3G_SOFTWARE_PATH'), 'beams',
                        'mapmaking_scripts', 
                        '1500d_3band_10sigma_ptsrc.txt')
plocs = np.loadtxt(pts_mask)

# which points to make maps for
points = range(10)

aux_input_files = [config, pts_mask]

# Hardware requirements
request_disk = 4*core.G3Units.GB
request_memory = 2*core.G3Units.GB

##################################################

# These requirements were adapted from a similar string in doc/osg/osg_guide.md
# It prevents the same machine from being used on re-submissions,
# requires remote machines to use RHEL 7, and forbids NPX and IIT_CE1,
# both of which have given trouble in the past.
requirements = (
    '((HAS_CVMFS_spt_opensciencegrid_org) && ' +
     '(GLIDEIN_ResourceName =!= "NPX") && ' +
     '(OSGVO_OS_STRING == "RHEL 7") && ' +
     '(GLIDEIN_ResourceName =!= "IIT_CE1")) && (')

# for retries
if args.retries > 0:
    requirements += '&& ('
for r in range(1, args.retries+1):
    requirements += '((TARGET.GLIDEIN_ResourceName =!= '
    requirements += 'MY.MachineAttrGLIDEIN_ResourceName{}) || '.format(r)
    requirements += '(RCC_Factory == "ciconnect"))'
    requirements += ' && 'if r != args.retries else ')'

for point in points:
    print('point {}'.format(point))
    ra, dec = plocs[point][1], plocs[point][2]
    all_obs = []
    # Other inputs to the script
    for field in ['ra0hdec-44.75', 'ra0hdec-52.25', 'ra0hdec-59.75', 
                  'ra0hdec-67.25']:
        obs_file = np.load(ps_list.format(field))
        all_obs += obs_file[str(point)].tolist()
    n_obs = len(all_obs)
    print('nobs: ',n_obs)
    if args.files_per_job == 0:
        n_split = n_obs
    else:
        n_split = args.files_per_job
    for m, start in enumerate(np.arange(0, n_obs, n_split)):
        # Description of job that will be appended to outputs
        job_root = 'ps{}_map{}'.format(point, m)
        # now, intersperse absolute path to cal frames
        cal0 = None
        abs_obs_with_cal = []
        all_cals = []
        obs_files = all_obs[start: start + n_split]
        for obs in obs_files:
            d = os.path.dirname(obs)
            cal = os.readlink(os.path.join(d, 'offline_calibration.g3'))
            if os.path.exists(cal):
                if cal == cal0:
                    abs_obs_with_cal.append(obs)
                else:
                    abs_obs_with_cal.append(cal)
                    abs_obs_with_cal.append(obs)
                    all_cals.append(cal)
                    cal0 = cal

        if len(all_cals) == 0:
            continue
        # on the grid, want relative paths
        root_dir_obs = os.path.commonpath(obs_files)
        root_dir_cals = os.path.commonpath(all_cals)
        rel_obs_with_cal = []
        for f in abs_obs_with_cal:
            if 'cal' in f:
                p = os.path.relpath(f, root_dir_cals)
                if p == '.':
                    p = os.path.basename(f)
                rel_obs_with_cal.append(p)
            else:
                p = os.path.relpath(f, root_dir_obs)
                if p == '.':
                    p = os.path.basename(f)
                rel_obs_with_cal.append(p)
        out_root = os.path.join(args.output_path, 'ps{}'.format(point))
        out_file = 'ps{}_map{:03}.g3'.format(point, m)
        extra_args = '--map-center-ra {} --map-center-dec {}'.format(ra, dec)
        extra_args += ' --config-file {}'.format(os.path.basename(config))
        extra_args += ' -o {}'.format(out_file)
        jobname = job_root
        # Need to use relative paths here, otherwise the remote machine will try
        # to find the input files at their local absolute path
        script_args = [' '.join(rel_obs_with_cal), extra_args]
        # don't submit jobs for maps that already finished
        if os.path.exists(os.path.join(out_root, out_file)):
            continue
        condor_submit(script, create_only=not args.submit, args=script_args,
                      log_root=os.path.join(args.log_dir, jobname),
                      output_root=out_root, retry=args.retries,
                      jobname=jobname, 
                      input_files=abs_obs_with_cal,
                      input_files_dest=rel_obs_with_cal,
                      aux_input_files=aux_input_files,
                      output_files=[out_file],
                      requirements=requirements,
                      request_disk=request_disk,
                      request_memory=request_memory,
                      grid_proxy=args.grid_proxy)
