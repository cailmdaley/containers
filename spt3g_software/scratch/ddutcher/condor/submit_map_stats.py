'''
submit_map_stats.py

This is ddutcher's script to submit 1500d map stats jobs.
It is essentially a wrapper to condor_tools.condor_submit(),
and submits jobs individually via a for loop.

Currently configured to run on only one field at a time.
'''

import os
import numpy as np
from glob import glob
import argparse
from spt3g import core
from spt3g.cluster.condor_tools import condor_submit

parser = argparse.ArgumentParser()
parser.add_argument('-s','--source', action ='store')
parser.add_argument('--submit', action = 'store_true')
parser.add_argument('--overwrite', action ='store_true')
pargs = parser.parse_args()

assert(pargs.source is not None)

# Description of job that will be appended to outputs
job_root = 'w1-4Hz_20190310'
    
# Search for all (downsampled, fullrate) data directories of interest
maps = sorted(glob(os.path.join('/spt/user/ddutcher',pargs.source,job_root,
                                '*',job_root+'*.g3')))
#print(maps)

condor_dir = os.path.join('/scratch/ddutcher/condor_logs', 'map_stats_'+job_root)
out_root = os.path.join('/spt/user/ddutcher', pargs.source, 'map_stats')
script = '/home/ddutcher/code/field_maps/map_stats.py'

test = True
if pargs.submit:
    test = False

# These requirements were adapted from a similar string in doc/osg/osg_guide.md
# It prevents the same machine from being used on re-submissions,
# requires remote machines to use RHEL 7, and forbids NPX and IIT_CE1,
# both of which have given trouble in the past.
requirements = '''( ( HAS_CVMFS_spt_opensciencegrid_org ) && ( ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName1 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName2 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName3 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName4 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName5 ) || ( RCC_Factory == "ciconnect" ) ) ) && ( OSGVO_OS_STRING == "RHEL 7" ) && (GLIDEIN_ResourceName =!= "NPX"))'''
    
for obs in maps:
    obsid = obs.split('_')[-1].replace('.g3','')
    jobname = pargs.source+'_'+obsid+'_mapstats'
    if not pargs.overwrite:
        if os.path.isfile(os.path.join('/spt/user/ddutcher', pargs.source,
                                       'map_stats', jobname+'.pkl')):
            continue
        if os.path.isfile(os.path.join('/home/ddutcher/data/2018/map_stats/weight1-4Hz_20190310',
                                       jobname+'.pkl')):
            continue
    infiles = [obs]
    # Need to use basename here, otherwise the remote machine will try
    # to find the input files at their local absolute path
    args_in = [os.path.basename(dat) for dat in infiles]
    extra_args=''
    args = '{infiles} -o {outfile} {extra}'.format(infiles = ' '.join(args_in),
                                                   outfile = jobname+'.pkl',
                                                   extra=extra_args)

    condor_submit(script, create_only=test, args = [args],
                  log_root = condor_dir,
                  output_root = out_root,
                  verbose = False,
                  retry = True,
                  jobname = jobname,
                  input_files= infiles,
                  output_files=[jobname+'.pkl'],
                  requirements = requirements,
                  request_disk=8*core.G3Units.GB,
                  request_memory=4*core.G3Units.GB,
                  grid_proxy='/home/ddutcher/ddutcher_proxy')
    