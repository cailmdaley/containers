'''
build_y1ee_dag.py

This is ddutcher's script to submit 1500d subfield map jobs
for the year 1 EE/TE analysis.
It is essentially a wrapper to condor_tools.condor_submit(),
and submits jobs via a dag file.
'''
import os
import numpy as np
from glob import glob
import argparse
from spt3g import core
from spt3g.std_processing import status_db
from spt3g.cluster.condor_tools import condor_submit

parser = argparse.ArgumentParser()
parser.add_argument('-s','--source', action ='store', default='ra0hdec-*',
                   help = 'The name of the source field')
parser.add_argument('--submit', action = 'store_true', default=False,
                   help = 'Create jobs and also submit them')
parser.add_argument('--overwrite', action ='store_true', default=False,
                   help = 'Re-analyze data and overwrite outputs')
parser.add_argument('--test', action='store_true', default=False,
                    help = 'Only submit one job as a test')
pargs = parser.parse_args()

db = status_db.ScanifyDatabase(read_only=True)

try:
    spt3g_software = os.environ['SPT3G_SOFTWARE_PATH']
except KeyError as e:
    print('%s environment variable not set'%e)
    raise

##################################################
########## JOB AND USER SPECIFIC SETTINGS ########

bands = ['90GHz','150GHz','220GHz']
resps = ['high','low']
directions = ['left','right']

# Location of script to be run
script = os.path.join(spt3g_software,'std_processing/mapmakers/master_field_mapmaker.py')

aux_input_files = ['/home/ddutcher/code/spt3g_software/std_processing/mapmakers/y1_ee.yaml']

output_files = ['{jobname}.g3.gz', 'simstub_{jobname}.g3.gz']

# Search for all (downsampled, fullrate) data directories of interest
datarate = 'downsampled'
all_obs = sorted(glob(os.path.join('/spt/data/bolodata/',datarate, pargs.source,'4*')) + 
                 glob(os.path.join('/spt/data/bolodata/',datarate, pargs.source,'5*')))

# Location for logging
dag_script = '/scratch/ddutcher/condor_logs/y1_ee/y1_ee_20190811.dag'

# Grid proxy location
grid_proxy = '/home/ddutcher/ddutcher_proxy'

# Hardware requirements
request_disk = 15*core.G3Units.GB
request_memory = 2*core.G3Units.GB

##################################################

# Filter list of all observations to remove 1) Known bad observations
# 2) Observations without a calframe. 3) Observations which we've already processed
# 4) Observations that aren't ready. Yes theres's some redundancy here.
badObs = np.loadtxt(os.path.join(spt3g_software,'scratch/lowitz/badObs/badObs.csv'),
                    delimiter = ',', skiprows = 1, dtype = 'str')
good_obs = []
for ind, path in enumerate(all_obs):
    obsid = path.split('/')[-1].replace('.g3','')
    source = path.split('/')[-2]
    if obsid in badObs[:,0]:
        continue
    if not db.check_ready(source, obsid, rate = datarate):
        continue
    cal = glob(os.path.join('/spt/user/production/calibration/calframe',
                            source, obsid+'.g3'))
    if len(cal)!=0:
        good_obs.append((obsid, source))
print([obsid for obsid, source in good_obs])

# These requirements were adapted from a similar string in doc/osg/osg_guide.md
# It prevents the same machine from being used on re-submissions,
# requires remote machines to use RHEL 7, and forbids NPX.
requirements = '''( ( HAS_CVMFS_spt_opensciencegrid_org ) && ( ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName1 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName2 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName3 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName4 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName5 ) || ( RCC_Factory == "ciconnect" ) ) ) && ( OSGVO_OS_STRING == "RHEL 7" ) && (GLIDEIN_ResourceName =!= "NPX"))'''

# If running in test mode, just do one job
if pargs.test:
    good_obs = good_obs[:1]
    bands = bands[:1]
    resps = resps[:1]
    directions = directions[:1]

# Create and submit jobs
if not os.path.exists(os.path.dirname(dag_script)):
    os.mkdir(os.path.dirname(dag_script))

f = open(dag_script,'w')
for obs, source in good_obs:
    for band in bands:
        for resp in resps:
            for direction in directions:
                extra_args = '--config-file y1_ee.yaml --produce-simstub'
                extra_args += ' --split-left-right %s --bands-to-use %s'%(direction, band)
                if resp == 'low':
                    extra_args += ' --wafers-to-include w177 w180 w174 w172'
                elif resp=='high':
                    extra_args += ' --wafers-to-include w181 w188 w203 w176'
                job_root =  '%s_%s_%s_maps'%(resp, band, direction)
                condor_dir = os.path.join(
                    '/scratch/ddutcher/condor_logs', 'y1_ee_20190811', job_root)
                out_root = os.path.join(
                    '/spt/user/ddutcher', '{source}', 'y1_ee_20190811', job_root)
                if not pargs.overwrite:
                    test = glob(os.path.join(out_root.format(source = source), obs,'*.g3*'))
                    if len(test)!=0:
                        continue
                jobname = job_root+'_'+obs
                # Find data to use:
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
                # Need to use basename here, otherwise the remote machine will try
                # to find the input files at their local absolute path
                args_in = [os.path.basename(dat) for dat in infiles]
                args = '{infiles} -o {outfile} {extra}'.format(infiles = ' '.join(args_in),
                                                               outfile = jobname+'.g3.gz',
                                                               extra = extra_args)

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

                f.write('JOB %s %s\n'%(jobname, os.path.join(condor_dir, obs, jobname+'.submit')))
                f.write('RETRY %s 5\n'%jobname)
f.close()

if pargs.submit:
    os.system('condor_submit_dag -update_submit -maxidle 25 %s'%dag_script)