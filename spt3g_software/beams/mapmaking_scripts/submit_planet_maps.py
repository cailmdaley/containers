"""
This script looks for observations within the range specified for the 
planet specified and submits jobs to the grid running the planet_maps.py 
script. 

Additional options for planet_maps.py can be passed to this script,
and will be passed through to that script.

To not just write submittable files to disk, make sure to
condor_submit_planet_maps.py --submit

Example: Make a first-order poly-filtered mars map.
python submit_planet_maps.py --submit --poly-order 1 
--outdir /spt/user/agambrel/beams/ --out-tag mars_poly1 --min-obsid 54868035
--max-obsid 54868036

Then make maps third-order poly-filtered mars map per 2018 Mars obs using the
poly1 map to make the point source mask, and make simstubs
python submit_planet_maps.py --submit --poly-order 3 
--outdir /spt/user/agambrel/beams/ --out-tag mars_poly3 
--mask-map /spt/user/agambrel/beams/54868035_mars_poly1.g3 --produce-simstub

This script will not work for sims and shouldn't be used for that purpose since
planet observation sims are very fast. Just call planet_maps.py with the same
poly and mask options as previously but point the script to the sim-stub instead
of timestreams.
"""

import os
import argparse
import glob
import getpass
from spt3g.cluster.condor_tools import condor_submit
from spt3g import core

parser = argparse.ArgumentParser()
parser.add_argument('--submit', action='store_true', default=False,
                    help='Submit jobs. Otherwise, just create files.')
parser.add_argument('--log-dir', action='store', 
                    default=os.path.join('/scratch', getpass.getuser(),
                                         'condor_logs'),
                    help='Absolute path to condor log directory')
parser.add_argument('--out-dir', action='store', 
                    default=os.path.join('/spt', 'user', getpass.getuser()),
                    help='Absolute path to output maps')
parser.add_argument('--out-tag', action='store', 
                    default='mars',
                    help='Tag to add to each map file name written to disk')
parser.add_argument('--source', action='store', default='mars',
                    help='Which planet pixel-raster')
parser.add_argument('--mask-map', action='store', default=None,
                    help='Absolute path to map to use for point source mask')
parser.add_argument('--retries', action='store', type=int, default=0,
                    help='Number of times to retry a job before giving up')
parser.add_argument('--min-obsid', action='store', type=int, default=54868035,
                    help='Minimimum obsid to make a map for')
parser.add_argument('--max-obsid', action='store', type=int, default=56489370,
                    help='Maximum obsid to make a map for')
parser.add_argument('--grid-proxy', action='store', 
                    default=os.path.join(os.getenv('HOME'), 
                                         getpass.getuser()+'_proxy'),
                    help='Path to grid proxy')      
                    
args, extra_args = parser.parse_known_args()
if len(extra_args):
    print('Passing arguments {} to planet_maps.py'.format(extra_args))

script = os.path.abspath('planet_maps.py')
data_dir = os.path.join('/spt', 'data', 'bolodata', 'fullrate',
                        '{}-pixelraster'.format(args.source))

cal_dir = os.path.join('/spt', 'user', 'production', 'calibration', 
                       'calframe', '{}-pixelraster'.format(args.source))

# standard requirements
requirements = (
    '((HAS_CVMFS_spt_opensciencegrid_org) && ' +
     '(GLIDEIN_ResourceName =!= "NPX") && ' +
     '(OSGVO_OS_STRING == "RHEL 7") || (RCC_Factory == "ciconnect"))')

# for retries
if args.retries > 0:
    requirements += '&& ('
for r in range(1, args.retries+1):
    requirements += '((TARGET.GLIDEIN_ResourceName =!= '
    requirements += 'MY.MachineAttrGLIDEIN_ResourceName{}) || '.format(r)
    requirements += '(RCC_Factory == "ciconnect"))'
    requirements += ' && 'if r != args.retries else ')'

# submit a job for each obsid between min_obsid and max_obsid in data_dir
for obsid in os.listdir(data_dir):
    # only use if it's in range specified
    if int(obsid) < args.min_obsid or int(obsid) > args.max_obsid:
        continue

    tod_files = sorted(glob.glob(os.path.join(data_dir, obsid, '[0-9]*.g3')))
    cal_file = os.path.join(cal_dir, '{}.g3'.format(obsid))

    # infiles are absolute paths to every file transfered to the grid
    infiles = [cal_file] + tod_files

    # arguments to script cannot be absolute paths, just file names
    args_in = [os.path.basename(dat) for dat in infiles]

    output_file = '{}_{}.g3'.format(args.out_tag, obsid)
    args_in.append('--output {}'.format(output_file))
    if args.mask_map is not None:
        infiles += [args.mask_map]
        args_in.append('--mask-map {}'.format(os.path.basename(args.mask_map)))

    # assume any unrecognized arguments to this script are arguments for 
    # planet_maps.py
    args_in += extra_args
    jobname = '{}_{}'.format(args.out_tag, obsid)

    condor_submit(script, create_only=not args.submit, args=args_in, 
                  log_root=os.path.join(args.log_dir, jobname),
                  retry=args.retries, output_root=args.out_dir, 
                  jobname=jobname, input_files=infiles, 
                  output_files=[output_file], requirements=requirements, 
                  request_disk=16*core.G3Units.GB,
                  request_memory=6*core.G3Units.GB,
                  grid_proxy=args.grid_proxy)
