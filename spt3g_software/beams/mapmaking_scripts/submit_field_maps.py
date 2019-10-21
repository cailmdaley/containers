'''
submit_field_maps.py

This script submits jobs using master_field_mapmaker.py for either SPT data
or Planck mock observations. Any modifications to mapmaking parameters must
be done in the yaml config file, since those arguments overwrite command line
args. If hm arg is specified, will make a mock observed Planck map of this 
half-mission instead of data. 
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
parser.add_argument('--config-file', default='field_maps.yaml',
                    help='The name of the yaml config file')
parser.add_argument('-s','--source', action ='store', nargs='+',
                    default=['ra0hdec-44.75', 'ra0hdec-52.25', 
                             'ra0hdec-59.75', 'ra0hdec-67.25'],
                    help='The name of the source field')
parser.add_argument('--min-obsid', action='store', type=int, default=47000000,
                    help='Minimimum obsid to make a map for')
parser.add_argument('--max-obsid', action='store', type=int, default=57999999,
                    help='Maximum obsid to make a map for')
parser.add_argument('--fullrate', action='store_true', default=False,
                    help='Use full rate bolometer timestreams')
parser.add_argument('--band', action ='store', nargs='+', 
                    default=['90', '150', '220'], 
                    help='The band to observe')
parser.add_argument('--hm', action ='store', nargs='+', 
                    default=None, choices=[1, 2], type=int,
                    help='The Planck half-mission to observe')
parser.add_argument('--ptsrc-mask', action='store',
                    default=os.path.join(os.getenv('SPT3G_SOFTWARE_PATH'),
                                         'sources', 
                                         '1500d_ptsrc_3band_50mJy.txt'),
                    help='List of point sources to mask')
parser.add_argument('-o', '--output-path', action='store', 
                    default=os.path.join('/spt', 'user', getpass.getuser()),
                    help='Top level subdirectory name for map run')
parser.add_argument('--map-run-name', action='store', default='field_maps',
                    help='Name of this map run. Field maps will be in '
                    'args.output_path/args.map_run_name, mock maps in, eg, '
                    'args.output_path/args.map_run_name_planck_hm1')
args = parser.parse_args()

# Planck map to use per band
pmap = {'90': '100', '150': '143', '220': '217'}
config = os.path.abspath(args.config_file)

# Location of script to be run
script = os.path.join(os.getenv('SPT3G_SOFTWARE_PATH'), 'std_processing',
                      'mapmakers', 'master_field_mapmaker.py')

# Hardware requirements
request_disk = 5*core.G3Units.GB
request_memory = 4*core.G3Units.GB

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

if args.fullrate:
    sample_rate = 'fullrate'
else:
    sample_rate = 'downsampled'

if args.hm is not None:
    # Submit a mock run for each simstub
    # make sure directory with simstubs exists
    if not os.path.exists(os.path.join(args.output_path, args.map_run_name)):
        raise IOError('Must have simstubs on disk to run mocks. '
            'None found in {}'.format(os.path.join(args.output_path), 
                                      args.map_run_name))
    mname = 'HFI_SkyMap_{}_2048_R3.01_halfmission-{}_cut_C_G3Units_hpl300.fits'

    for band in args.band:
        # need to split jobs by band since need different Planck map to mock
        # per band
        for hm in args.hm:
            job_root = '{}GHz_HM{}_{}'.format(band, hm, args.map_run_name)
            # Other inputs to the script
            sim_map = os.path.join('/spt', 'user', 'agambrel', 'planck_maps',
                                   mname.format(pmap[band], hm))
            extra_args = ['-m {}'.format(os.path.basename(sim_map)),
                          '--bands-to-use {}'.format(band),
                          '--sim', 
                          '--config-file {}'.format(os.path.basename(config))]
            for source in args.source:
                # Search for all data directories of interest
                all_obs = sorted(
                    glob(os.path.join(args.output_path, args.map_run_name, 
                                      source, '*', 'simstub*.g3')))
                if len(all_obs) == 0:
                    print('No sim stubs to produce maps for in {}'.format(
                        os.path.join(args.output_path, args.map_run_name, 
                                     source)))
                for obs in all_obs:
                    obsid = obs.split('/')[-2]
                    jobname = job_root+'_'+obsid

                    script_args = [os.path.basename(obs), 
                                   '-o {}'.format(jobname+'.g3')]
                    script_args += extra_args
                    # don't submit jobs for maps that already finished
                    out_root = os.path.join(
                        args.output_path, 
                        args.map_run_name+'_planck_hm{}'.format(hm),
                        source, obsid)
                    out_file = os.path.join(out_root, jobname + '.g3')
                    if os.path.exists(out_file):
                        continue
                    condor_submit(script, create_only=not args.submit, 
                                  args=script_args,
                                  log_root=os.path.join(args.log_dir, jobname),
                                  output_root=out_root, retry=args.retries,
                                  jobname=jobname, 
                                  input_files=[obs, sim_map],
                                  aux_input_files=[config, args.ptsrc_mask],
                                  output_files=['{}.g3'.format(jobname)],
                                  requirements=requirements,
                                  request_disk=request_disk,
                                  request_memory=request_memory,
                                  grid_proxy=args.grid_proxy)
else:
    # Submit a data map for each obsid in range given
    job_root = args.map_run_name
    # Other inputs to the script-- always produce a simstub
    extra_args = ['--bands-to-use {}'.format(' '.join(args.band)),
                  '--config-file {}'.format(os.path.basename(config)),
                  '--produce-simstub']
    for source in args.source:
        # Search for all data directories of interest
        all_obs = sorted(
                    glob(os.path.join('/spt', 'data', 'bolodata', sample_rate,
                                      source, '*')))
        for obs in all_obs:
            # check if obs is in range
            obsid = obs.split('/')[-1]
            if np.logical_or(int(obsid) <= args.min_obsid, 
                             int(obsid) >= args.max_obsid):
                continue
            jobname = job_root+'_'+obsid
            cal = os.readlink(os.path.join(obs, 'offline_calibration.g3'))
            # sometimes, the symlink is broken- there is no cal file. don't use
            # these observations
            if not os.path.exists(cal):
                continue
            dat = sorted(glob(os.path.join(obs, '[0-9]*.g3')))
            infiles = [cal] + dat
            infiles_rel = [os.path.basename(f) for f in infiles]
            script_args = infiles_rel + ['-o {}'.format(jobname+'.g3')]
            script_args += extra_args
            # don't submit jobs for maps that already finished
            out_root = os.path.join(args.output_path, args.map_run_name,
                                    source, obsid)
            out_file = os.path.join(out_root, jobname + '.g3')
            if os.path.exists(out_file):
                continue
            condor_submit(script, create_only=not args.submit, 
                          args=script_args,
                          log_root=os.path.join(args.log_dir, jobname),
                          output_root=out_root, retry=args.retries,
                          jobname=jobname, 
                          input_files=infiles,
                          aux_input_files=[config, args.ptsrc_mask],
                          output_files=['{}.g3'.format(jobname),
                                        'simstub_{}.g3'.format(jobname)],
                          requirements=requirements,
                          request_disk=request_disk,
                          request_memory=request_memory,
                          grid_proxy=args.grid_proxy)
                
