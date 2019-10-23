import numpy as np
import os
import argparse
import glob
from spt3g.cluster.condor_tools import condor_submit
from spt3g import core

"""
Submit mars maps to condor. Call with --submit option to actually
submit jobs to the grid. Otherwise, will run in test mode and write
submit scripts to log directories.
"""

# Use x/y detector centroids calculated from per-bolo Mars maps instead
# of those in offline calibration file
mars_offs = False

# Flag saturated detector scans
flag_sat = False

poly_order = 3

parser = argparse.ArgumentParser()
parser.add_argument('--submit', action='store_true')

args = parser.parse_args()

condor_dir = '/scratch/agambrel/condor_logs'
out_dir = os.path.join('/spt', 'user', 'agambrel', 'beams', '201903',
                       'mars_maps_poly{}'.format(poly_order))
script = os.path.join(os.getenv('SPT3G_SOFTWARE_PATH'), 'scratch', 'agambrel',
                      'beams', 'make_mars_maps.py')
data_dir = os.path.join('/spt', 'data', 'bolodata', 'fullrate',
                       'mars-pixelraster')

# If not using a map to get the point source mask, set poly1_dir = None
poly1_dir = os.path.join('/spt', 'user', 'agambrel', 'beams', '201903',
                         'mars_maps_oldgaussfits_poly1')

test=True
if args.submit:
    test = False

more_args = ''
if flag_sat:
    tagf = "_flagsat"
    more_args += " --flag-saturated"
else:
    tagf = ""

if mars_offs:
    cal_dir = os.path.join('/spt', 'user', 'agambrel', 'beams', 'fits',
                           'gauss_fits_condor')

    tagm = "_marsoffs"
    more_args += " --mars-offs"
else:
    cal_dir = os.path.join('/spt', 'user', 'production', 'calibration', 
                           'calframe', 'mars-pixelraster')
    tagm = ""

# standard requirements
requirements = (
    '((HAS_CVMFS_spt_opensciencegrid_org) && ' +
     '(GLIDEIN_ResourceName =!= "NPX") && ' +
     '(OSGVO_OS_STRING == "RHEL 7") || (RCC_Factory == "ciconnect")) && (')

# non-standard... I don't know what these are. Stolen from DDutcher
# repeat 5 times if not successful?
for r in range(1, 5):
    requirements += '((TARGET.GLIDEIN_ResourceName =!= '
    requirements += 'MY.MachineAttrGLIDEIN_ResourceName{}) || '.format(r)
    requirements += '(RCC_Factory == "ciconnect"))'
    requirements += ' && 'if r != 4 else ')'

obsids = ['54868035', '55383550', '55211360',
          '55643104', '56489370']

# submit a job for each mars observation after sept 28 in 2018
for mars_obs in obsids:
    map_files = sorted(glob.glob(os.path.join(data_dir, mars_obs, '[0-9]*.g3')))
    if mars_offs:
        cal_file = os.path.join(cal_dir, 
                                '{}_gauss_fit_nooutliers.g3'.format(mars_obs))
    else:
        cal_file = os.path.join(cal_dir, '{}.g3'.format(mars_obs))
    if poly1_dir is not None:
        poly1_map = os.path.join(poly1_dir, '{}_poly1.g3'.format(mars_obs))
    output = '{}_poly{}{}{}.g3'.format(mars_obs, poly_order, tagm, tagf)
    jobname = 'mars_map_{}{}{}'.format(mars_obs, tagm, tagf)
    print(jobname)
    # infiles are absolute paths to files that will all be transferred
    # to a single directory on the grid along with the script..
    infiles = [cal_file] + map_files

    # args_files are the relative paths those files will have on the grid,
    # so just the file names.
    args_files = [os.path.basename(dat) for dat in infiles]
    args_in = [' '.join(args_files), 
               '--output ' + output,
               '--poly-order {}'.format(poly_order), more_args]
    if poly1_dir is not None:
        infiles += [poly1_map]
        args_in.append('--poly1-map {}'.format(os.path.basename(poly1_map)))

    condor_submit(
        script, create_only=test, args=args_in, 
        log_root=os.path.join(condor_dir, 
                              mars_obs+'poly{}_201903'.format(poly_order)), 
        retry=5, output_root=out_dir, jobname=jobname, input_files=infiles, 
        output_files=[output], requirements=requirements, 
        request_disk=16*core.G3Units.GB, request_memory=6*core.G3Units.GB,
        grid_proxy='/home/agambrel/agambrel_proxy')
