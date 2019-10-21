'''
build_y1ee_dag.py

This is ddutcher's script to submit 1500d subfield simmap jobs
for the year 1 EE/TE analysis.
It is essentially a wrapper to condor_tools.condor_submit(),
and submits jobs via a dag file.
'''
import os
import numpy as np
from glob import glob
import argparse
from spt3g import core
from spt3g.cluster.condor_tools import condor_submit

parser = argparse.ArgumentParser()
parser.add_argument('sim_map_index', action='store', type=str, choices=['1','2','full'],
                   help='The Planck mission to observe.')
parser.add_argument('--sim-map-dir', action='store',
                    default='/spt/user/agambrel/planck_maps',
                    help='Directory containing the sim maps')
parser.add_argument('--bands', nargs='+', action='store', type=str,
                    default=['90', '150', '220'])
parser.add_argument('--dag-script', type=str, default='out.dag',
                    help='Location to put dagfile')
parser.add_argument('--tag', type=str,
                    help='Some label for the sim and output')
parser.add_argument('--submit', action = 'store_true', default=False,
                   help = 'Create jobs and also submit them')
parser.add_argument('--overwrite', action ='store_true', default=False,
                   help = 'Re-analyze data and overwrite outputs')
parser.add_argument('--test', action='store_true', default=False,
                    help = 'Only submit one job as a test')
pargs = parser.parse_args()

try:
    spt3g_software = os.environ['SPT3G_SOFTWARE_PATH']
except KeyError as e:
    print('%s environment variable not set' % e)
    raise
##################################################
########## JOB AND USER SPECIFIC SETTINGS ########
resps = ['high', 'low']
directions = ['left', 'right']
# Location of script to be run
script = os.path.join(spt3g_software, 'std_processing/mapmakers/master_field_mapmaker.py')

aux_input_files = ['/home/ddutcher/code/spt3g_software/std_processing/mapmakers/y1_ee.yaml']

output_files = ['{jobname}.g3.gz']

# Location for logging
dag_script = pargs.dag_script
parent = os.path.dirname(os.path.abspath(dag_script))
if not os.path.exists(parent):
    os.makedirs(parent)

# Grid proxy location
grid_proxy = '/home/ddutcher/ddutcher_proxy'

# Hardware requirements
# Simstub is 21MB, simsky is 352 MB, map is 40MB
request_disk = 1*core.G3Units.GB
request_memory = 2*core.G3Units.GB

##################################################

# These requirements were adapted from a similar string in doc/osg/osg_guide.md
# It prevents the same machine from being used on re-submissions,
# requires remote machines to use RHEL 7, and forbids NPX.
requirements = '''( ( HAS_CVMFS_spt_opensciencegrid_org ) && ( ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName1 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName2 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName3 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName4 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName5 ) || ( RCC_Factory == "ciconnect" ) ) ) && ( OSGVO_OS_STRING == "RHEL 7" ) && (GLIDEIN_ResourceName =!= "NPX"))'''

if pargs.test:
    pargs.bands = pargs.bands[:1]
    resps = resps[:1]
    directions = directions[:1]


if pargs.sim_map_index in ['1', '2']:
    mp_pth = 'HFI_SkyMap_{0}_2048_R3.01_halfmission-{1}_cut_C_G3Units_hpl300.fits'
    label = 'hm'+pargs.sim_map_index
elif pargs.sim_map_index == 'full':
    mp_pth = 'HFI_SkyMap_{0}_2048_R3.01_{1}_cut_C_G3Units_hpl300.fits'
    label = pargs.sim_map_index
else:
    raise ValueError(pargs.sim_map_index)

f = open(dag_script,'w')
for band in pargs.bands:
    if band == '90': freq='100'
    if band == '150': freq='143'
    if band == '220': freq='217'

    sim_map = os.path.join(
        pargs.sim_map_dir, mp_pth.format(freq, pargs.sim_map_index))

    if not os.path.isfile(sim_map):
        raise FileNotFoundError(sim_map)

    for resp in resps:
        for direction in directions:
            pths = sorted(glob(os.path.join(
                '/spt/user/ddutcher/ra0hdec-*/y1_ee_20190811',
                '_'.join([resp, band+'GHz', direction, 'maps']),
                '*',
                'simstub_*.g3*')))
            if len(pths) == 0:
                raise FileNotFoundError(os.path.join(
                    '/spt/user/ddutcher/ra0hdec-*/y1_ee_20190811',
                    '_'.join([resp, band+'GHz', direction, 'maps']),
                    '*',
                    'simstub_*.g3*'))
            if pargs.test:
                pths = pths[:1]
            for simstub in pths:
                obs = simstub.split('/')[-2]
                source = simstub.split('/')[4]
                extra_args = '--config-file y1_ee.yaml'
                extra_args += ' -z --sim --sim-map %s' % os.path.basename(sim_map)
                extra_args += ' --split-left-right %s --bands-to-use %s' % (direction, band)
                if resp == 'low':
                    extra_args += ' --wafers-to-include w177 w180 w174 w172'
                elif resp=='high':
                    extra_args += ' --wafers-to-include w181 w188 w203 w176'

                job_root =  '%s_%s_%s_maps'%(resp,  band+'GHz', direction)
                condor_dir = os.path.join(
                    '/scratch/ddutcher/condor_logs/sims',
                    'planck_' + label,
                    job_root,
                )
                out_root = os.path.join(
                    '/spt/user/ddutcher', '{source}', 'y1_ee_20190811', job_root)
                
                if pargs.tag is not None:
                    jobname = '_'.join([
                        'planck_' + label, pargs.tag, job_root, obs])
                else:
                    jobname = '_'.join([
                        'planck_' + label, job_root, obs])
                if not pargs.overwrite:
                    test = glob(os.path.join(
                        out_root.format(source=source),
                        obs,
                        jobname + '*.g3*'))
                    if len(test)!=0:
                        continue

                # Need to use basename here, otherwise the remote machine will try
                # to find the input files at their local absolute path
                args_in = [os.path.basename(simstub)]
                args = '{infiles} -o {outfile} {extra}'.format(
                    infiles = ' '.join(args_in),
                    outfile = jobname+'.g3.gz',
                    extra = extra_args)

                condor_submit(script, create_only=True, args = [args],
                              log_root=os.path.join(condor_dir, obs),
                              output_root=os.path.join(
                                  out_root.format(source=source), obs),
                              verbose=False,
                              retry=False,
                              jobname=jobname,
                              input_files=[simstub, sim_map],
                              user_code='',
                              aux_input_files=aux_input_files,
                              output_files=[f.format(jobname=jobname) for f in output_files],
                              requirements=requirements,
                              request_disk=request_disk,
                              request_memory=request_memory,
                              grid_proxy=grid_proxy,
                              globus_uri='gsiftp://ceph-gridftp1.grid.uchicago.edu:2811/cephfs')

                f.write('JOB %s %s\n'%(jobname, os.path.join(condor_dir, obs, jobname+'.submit')))
                f.write('RETRY %s 5\n'%jobname)
f.close()

if pargs.submit:
    os.system('condor_submit_dag -update_submit -maxidle 500 %s'%dag_script)