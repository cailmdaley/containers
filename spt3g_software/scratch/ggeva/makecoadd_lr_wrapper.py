# To run this wrapper script
# python makecoadd_lr_wrapper.py -s [source] [-p] --submit [--overwrite]

import sys
sys.path.append('/home/ggeva/spt3g/spt3g_software/build')

'''
Taken from spt3g_software/scratch/ddutcher/condor/submit_1500d.py

It is essentially a wrapper to condor_tools.condor_submit(),
and submits jobs individually via a for loop.
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
parser.add_argument('-p', '--polarized', action='store_true',
           default=False, help='make polarized maps')
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
job_root = 'coadd_lr_0.25res'

# Location of script to be run
script = os.path.join(spt3g_software, 'scratch/ggeva/makecoadd_lr.py')
extra_args = '-x 3 -y 3 -r 0.25 --lr'
if pargs.polarized:
    extra_args += ' -p'
    job_root += '_pol'
aux_input_files = []
output_files = ['{jobname}.g3']

# Search for all (downsampled, fullrate) data directories of interest
all_obs = sorted(glob(
    os.path.join('/spt/data/bolodata/downsampled',pargs.source,'*')))
print('all_obs', all_obs)

# Location for logging and data outputs.These are directories on scott
condor_dir = os.path.join('/scratch/ggeva/condor_logs', pargs.source, job_root)
out_root = os.path.join('/spt/user/ggeva', pargs.source, job_root)
print(condor_dir)
print(out_root)

# Grid proxy location
grid_proxy = '/home/ggeva/ggeva_proxy'

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
        print('badobs: ' + obsid)
        continue
    if not db.check_ready(pargs.source, obsid, rate='downsampled'):
        print('db.check_ready failed: ' + obisd)
        continue
    if not pargs.overwrite:
        test = glob(os.path.join(out_root, obsid,'*.g3'))
        if len(test)!=0:
            print('len(test) = 0: ' + obisd)
            continue
    cal = glob(os.path.join('/spt/user/production/calibration/calframe',
                            pargs.source, obsid+'.g3'))
    if len(cal)!=0:
        good_obs.append(obsid)
    else:
        print('len(cal) = 0: ' + obsid)
print(good_obs)

# These requirements were adapted from a similar string in doc/osg/osg_guide.md
# It prevents the same machine from being used on re-submissions,
# requires remote machines to use RHEL 7, and forbids NPX and IIT_CE1,
# both of which have given trouble in the past.
requirements = '''( ( HAS_CVMFS_spt_opensciencegrid_org ) && ( ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName1 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName2 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName3 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName4 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName5 ) || ( RCC_Factory == "ciconnect" ) ) ) && ( OSGVO_OS_STRING == "RHEL 7" ) && (GLIDEIN_ResourceName =!= "NPX") && (GLIDEIN_ResourceName =!= "IIT_CE1"))'''
 
# Create and submit jobs
create_only = True
if pargs.submit:
    create_only = False
for obs in good_obs:
    jobname = job_root+'_'+obs
    # For a brief time, the downsampler had a bug.
    # Use the re-processed data in Sasha's directory for these:
    if int(obs)<40000000:
        continue
    if int(obs)>48731000 and int(obs) < 49231800:
        data = glob(os.path.join('/spt/user/arahlin/scanify/downsampled',pargs.source,
                                 obs,'0*.g3'))
    else:           
        data = glob(os.path.join('/spt/data/bolodata/downsampled',pargs.source, 
                                 obs,'0*.g3'))
    cal = os.path.join('/spt/user/production/calibration/calframe/',pargs.source,
                       obs+'.g3')
    infiles = [cal] + sorted(data)
    # Need to use basename here, otherwise the remote machine will try
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
                  request_disk=15*core.G3Units.GB,
                  request_memory=2*core.G3Units.GB,
                  grid_proxy=grid_proxy)