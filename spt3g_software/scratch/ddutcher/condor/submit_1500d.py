'''
submit_1500d.py

This is ddutcher's script to submit 1500d subfield map jobs.
It is essentially a wrapper to condor_tools.condor_submit(),
and submits jobs individually via a for loop.

Currently configured to run on only one field at a time.
'''
import os
import numpy as np
from glob import glob
import argparse
from spt3g import core, std_processing
from spt3g.cluster.condor_tools import condor_submit

parser = argparse.ArgumentParser()
parser.add_argument('-s','--source', action ='store',
                   help = 'The name of the source field')
parser.add_argument('--submit', action = 'store_true',
                   help = 'Create jobs and also submit them')
parser.add_argument('--overwrite', action ='store_true',
                   help = 'Re-analyze data and overwrite outputs')
pargs = parser.parse_args()

assert(pargs.source is not None)
db = std_processing.status_db.ScanifyDatabase(read_only=True)

try:
    spt3g_software = os.environ['SPT3G_SOFTWARE_PATH']
except KeyError as e:
    print('%s environment variable not set'%e)
    raise

##################################################
########## JOB AND USER SPECIFIC SETTINGS ########

# Description of job that will be appended to outputs
job_root = 'default_maps_2018'

# Location of script to be run
script = os.path.join(spt3g_software,'std_processing/mapmakers/master_field_mapmaker.py')

# List of auxilliary input file locations,
# e.g. small text files like config files or point source lists
aux_input_files = [os.path.join(spt3g_software, 'std_processing/mapmakers/default_config.yaml')]

# Command line args to be passed to the script on remote machine.
# Note that any file paths need to be those on the remote machine.
extra_args = '--config-file default_config.yaml --produce-simstub'

# List of output file name formats. These will be filled in with
# the specific jobname later on.
# This is what tells condor what files to transfer back.
output_files = ['{jobname}.g3', 'simstub_{jobname}.g3']

# Locate input data.
# Search for all (downsampled, fullrate) data directories of interest
all_obs = sorted(glob(
    os.path.join('/spt/data/bolodata/downsampled',pargs.source,'5*')))

# Location for logging and data outputs
condor_dir = os.path.join('/scratch/ddutcher/condor_logs', pargs.source, job_root)
out_root = os.path.join('/spt/user/ddutcher', pargs.source, job_root)

# Grid proxy location
grid_proxy = '/home/ddutcher/ddutcher_proxy'

# Hardware requirements
request_disk = 15*core.G3Units.GB
request_memory = 5*core.G3Units.GB

##################################################

# Filter list of all observations to remove 1) Known bad observations
# 2) Observations without a calframe. 3) Observations which we've already processed
# 4) Observations that aren't ready. Yes theres's some redundancy here.
badObs = np.loadtxt(os.path.join(spt3g_software,'scratch/lowitz/badObs/badObs.csv'),
                    delimiter = ',', skiprows = 1, dtype = 'str')
good_obs = []
for ind, path in enumerate(all_obs):
    obsid = path.split('/')[-1].replace('.g3','')
    if obsid in badObs[:,0]:
        continue
    if not db.check_ready(pargs.source, obsid, rate='downsampled'):
        continue
    if not pargs.overwrite:
        test = glob(os.path.join(out_root, obsid,'*.g3'))
        if len(test)!=0:
            continue
    cal = glob(os.path.join('/spt/user/production/calibration/calframe',
                            pargs.source, obsid+'.g3'))
    if len(cal)!=0:
        good_obs.append(obsid)
print(good_obs)

# These requirements were adapted from a similar string in doc/osg/osg_guide.md
# It prevents the same machine from being used on re-submissions,
# requires remote machines to use RHEL 7, and forbids NPX.
requirements = '''( ( HAS_CVMFS_spt_opensciencegrid_org ) && ( ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName1 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName2 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName3 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName4 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName5 ) || ( RCC_Factory == "ciconnect" ) ) ) && ( OSGVO_OS_STRING == "RHEL 7" ) && (GLIDEIN_ResourceName =!= "NPX"))'''

# Create and submit jobs
create_only = True
if pargs.submit:
    create_only = False
for obs in good_obs:
    jobname = job_root+'_'+obs

    # Glob for all the data files for a given observation
    if int(obs)<40000000:
        continue
    if int(obs)>48731000 and int(obs) < 49231800:
        # For a brief time, the downsampler had a bug.
        # Use the re-processed data in Sasha's directory for these:
        data = glob(os.path.join('/spt/user/arahlin/scanify/downsampled',pargs.source,
                                 obs,'0*.g3'))
    else:           
        data = glob(os.path.join('/spt/data/bolodata/downsampled',pargs.source, 
                                 obs,'0*.g3'))

    # Find the calibration file for the observation
    # Note: use the actual file path and not the symlink
    cal = os.path.join('/spt/user/production/calibration/calframe/',pargs.source,
                       obs+'.g3')
    infiles = [cal] + sorted(data)

    # Need to use basename for args_in, otherwise the remote machine will try
    # to find the input files at their local absolute path
    args_in = [os.path.basename(dat) for dat in infiles]
    args = '{infiles} -o {outfile} {extra}'.format(infiles = ' '.join(args_in),
                                                   outfile = jobname+'.g3',
                                                   extra=extra_args)

    condor_submit(script, create_only=create_only, args = [args],
                  log_root = os.path.join(condor_dir, obs),
                  output_root = os.path.join(out_root, obs),
                  retry = True,
                  jobname = jobname,
                  input_files= infiles,
                  aux_input_files=aux_input_files,
                  output_files=[f.format(jobname=jobname) for f in output_files],
                  requirements = requirements,
                  request_disk=request_disk,
                  request_memory=request_memory,
                  grid_proxy=grid_proxy)