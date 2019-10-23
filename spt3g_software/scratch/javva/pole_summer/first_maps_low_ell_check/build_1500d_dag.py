'''
build_1500d_dag.py

This is ddutcher's script to submit 1500d subfield map jobs.
It is essentially a wrapper to condor_tools.condor_submit(),
and submits jobs via a dag file.
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

if pargs.source is None:
    pargs.source = 'ra0hdec-*'
db = std_processing.status_db.ScanifyDatabase(read_only=True)

try:
    spt3g_software = os.environ['SPT3G_SOFTWARE_PATH']
except KeyError as e:
    print('%s environment variable not set'%e)
    raise

##################################################
########## JOB AND USER SPECIFIC SETTINGS ########

# Description of job that will be appended to outputs
job_root = 'adam_gain_match'

# Location of script to be run
script = os.path.join(spt3g_software,'scratch/javva/free_time/gain_match/adam_gains.py')
# script = '/home/ddutcher/code/field_maps/field_autoproc.py'

extra_args = '--config-file test19.yaml'# --produce-simstub'

aux_input_files = ['/home/javva/spt3g_software/std_processing/mapmakers/test19.yaml',
                   '/home/ddutcher/code/spt3g_software/scratch/ddutcher/1500d_3band_10sigma_ptsrc.txt']

output_files = ['{jobname}.g3']#, 'simstub_{jobname}.g3']

# Search for all (downsampled, fullrate) data directories of interest
datarate = 'downsampled'
all_obs1 = sorted(glob(
    os.path.join('/spt/data/bolodata/',datarate, pargs.source,'6*')))
all_obs2 = sorted(glob(
    os.path.join('/spt/data/bolodata/',datarate, pargs.source,'7*')))
all_obs = all_obs2
# Location for logging and data outputs
condor_dir = os.path.join('/home/javva/grid/condor_logs', job_root)
out_root = os.path.join('/spt/user/javva/maps_2019', '{source}', job_root)
dag_script = os.path.join(condor_dir, job_root+ '.dag')

# Grid proxy location
grid_proxy = '/home/javva/javva_proxy'

# Hardware requirements
request_disk = 15*core.G3Units.GB
request_memory = 4*core.G3Units.GB

##################################################

# Filter list of all observations to remove 1) Known bad observations
# 2) Observations without a calframe. 3) Observations which we've already processed
# 4) Observations that aren't ready. Yes theres's some redundancy here.
badObs = np.loadtxt(os.path.join(spt3g_software,'scratch/lowitz/badObs/badObs.csv'),
                    delimiter = ',', skiprows = 1, dtype = 'str')
good_obs = []
print('checking over %s obs'%len(all_obs))
for ind, path in enumerate(all_obs):
    print(ind)
    obsid = path.split('/')[-1].replace('.g3','')
    source = path.split('/')[-2]
    if obsid in badObs[:,0]:
        continue
    if not db.check_ready(source, obsid, rate = datarate):
        continue
    if not pargs.overwrite:
        test = glob(os.path.join(out_root.format(source = source), obsid,'*.g3'))
        if len(test)!=0:
            continue
    cal = glob(os.path.join('/spt/user/production/calibration/calframe',
                            source, obsid+'.g3'))
    if len(cal)!=0:
        good_obs.append((obsid, source))
print([obsid for obsid, source in good_obs])
print('submitting maps for %s maps'%len(good_obs))
# These requirements were adapted from a similar string in doc/osg/osg_guide.md
# It prevents the same machine from being used on re-submissions,
# requires remote machines to use RHEL 7, and forbids NPX.
requirements = '''( ( HAS_CVMFS_spt_opensciencegrid_org ) && ( ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName1 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName2 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName3 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName4 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName5 ) || ( RCC_Factory == "ciconnect" ) ) ) && ( OSGVO_OS_STRING == "RHEL 7" ) && (GLIDEIN_ResourceName =!= "NPX"))'''

# Create and submit jobs
if not os.path.exists(condor_dir):
    os.mkdir(condor_dir)

f = open(dag_script,'w')
idx = -1
for obs, source in (good_obs):
    idx +=1
    print('writing line %s of %s'%(idx, len(good_obs)))
    jobname = job_root+'_'+obs
    # For a brief time, the downsampler had a bug.
    # Use the re-processed data in Sasha's directory for these:
    if int(obs)<40000000:
        continue
    if int(obs)>48731000 and int(obs)<49231800 and datarate=='downsampled':
        data = glob(os.path.join('/spt/user/arahlin/scanify/downsampled',source,
                                 obs,'0*.g3'))
    else:           
        data = glob(os.path.join('/spt/data/bolodata/',datarate,source, 
                                 obs,'0*.g3'))
    cal = os.path.join('/spt/user/production/calibration/calframe/',source,
                       obs+'.g3')
    infiles = [cal] + sorted(data)
    print(infiles)
    # Need to use basename here, otherwise the remote machine will try
    # to find the input files at their local absolute path
    args_in = [os.path.basename(dat) for dat in infiles]
    args = '{infiles} -o {outfile} {extra}'.format(infiles = ' '.join(args_in),
                                                   outfile = jobname+'.g3',
                                                   extra=extra_args)

    condor_submit(script, create_only=True, args = [args],
                  log_root = os.path.join(condor_dir, obs),
                  output_root = os.path.join(
                      out_root.format(source=source), obs),
                  verbose = False,
                  retry = False,
                  jobname = jobname,
                  input_files= infiles,
                  aux_input_files=aux_input_files,
                  output_files=[f.format(jobname=jobname) for f in output_files],
                  requirements = requirements,
                  request_disk=request_disk,
                  request_memory=request_memory,
                  grid_proxy=grid_proxy,
                  globus_uri='gsiftp://ceph-gridftp1.grid.uchicago.edu:2811/cephfs')
    print([f.format(jobname=jobname) for f in output_files])
    f.write('JOB %s %s\n'%(jobname, os.path.join(condor_dir, obs, jobname+'.submit')))
    f.write('RETRY %s 5\n'%jobname)
f.close()

if pargs.submit:
    os.system('condor_submit_dag -update_submit -maxidle 25 %s'%dag_script)

