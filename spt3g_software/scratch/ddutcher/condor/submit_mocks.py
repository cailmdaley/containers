'''
submit_mocks.py

This is ddutcher's script to submit 1500d mock-observation jobs.
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
parser.add_argument('-s','--source', action ='store',
                   help = 'The name of the source field')
parser.add_argument('--submit', action = 'store_true',
                   help = 'Create jobs and also submit them')
parser.add_argument('--overwrite', action ='store_true',
                   help = 'Re-analyze data and overwrite outputs')
parser.add_argument('-n', '--num', action='store',
                    help = 'The number of the sim to use')
pargs = parser.parse_args()

assert(pargs.source is not None)

try:
    spt3g_software = os.environ['SPT3G_SOFTWARE_PATH']
except KeyError as e:
    print('%s environment variable not set'%e)
    raise

##################################################    
########## JOB AND USER SPECIFIC SETTINGS ########

# Description of job that will be appended to outputs
job_root = 'sim%s_noW201wafCMpoly19mhpf300'%pargs.num

# Map run to be simulated
map_run = 'noW201_wafCM_poly19_mhpf300_lr'

# Location of script to be run
script = '/home/ddutcher/code/field_maps/mock_1500d.py'

# Search for all data directories of interest
all_obs = sorted(glob(os.path.join('/spt/user/ddutcher', pargs.source, map_run,'51*','simstub*.g3')))

# Other inputs to the script
sim_map = '/spt/user/arahlin/lenspix_maps/lensed_cmb_lmax7000_nside8192_interp0.3_method1_pol_1_sim_%s_lensed_map.fits'%pargs.num
extra_args = '-m %s -x 75 -y 50 --lr'%(os.path.basename(sim_map))
aux_input_files = []

# Format of output file names
output_files = ['{jobname}.g3']

# Location for logging and data outputs
condor_dir = os.path.join('/scratch/ddutcher/condor_logs', pargs.source, job_root)
out_root = os.path.join('/spt/user/ddutcher', pargs.source, job_root)

# Grid proxy location
grid_proxy = '/home/ddutcher/ddutcher_proxy'

# Hardware requirements
request_disk = 5*core.G3Units.GB
request_memory = 4*core.G3Units.GB

##################################################

# These requirements were adapted from a similar string in doc/osg/osg_guide.md
# It prevents the same machine from being used on re-submissions,
# requires remote machines to use RHEL 7, and forbids NPX and IIT_CE1,
# both of which have given trouble in the past.
requirements = '''( ( HAS_CVMFS_spt_opensciencegrid_org ) && ( ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName1 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName2 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName3 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName4 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName5 ) || ( RCC_Factory == "ciconnect" ) ) ) && ( OSGVO_OS_STRING == "RHEL 7" ) && (GLIDEIN_ResourceName =!= "NPX") && (GLIDEIN_ResourceName =!= "IIT_CE1"))'''
 
# Create and submit jobs
create_only = True
if pargs.submit:
    create_only = False
for obs in all_obs:
    obsid = obs.split('/')[-2]
    jobname = job_root+'_'+obsid
    # Need to use basename here, otherwise the remote machine will try
    # to find the input files at their local absolute path
    args = '{infiles} -o {outfile} {extra}'.format(infiles = os.path.basename(obs),
                                                   outfile = jobname+'.g3',
                                                   extra=extra_args)

    condor_submit(script, create_only=create_only, args = [args],
                  log_root = os.path.join(condor_dir, obsid),
                  output_root = out_root,
                  retry = True,
                  jobname = jobname,
                  input_files= [obs, sim_map],
                  aux_input_files=aux_input_files,
                  output_files=[f.format(jobname=jobname) for f in output_files],
                  requirements = requirements,
                  request_disk=request_disk,
                  request_memory=request_memory,
                  grid_proxy=grid_proxy)
