import os
from glob import glob
from spt3g.cluster.condor_tools import condor_submit
import argparse
from spt3g import core

parser = argparse.ArgumentParser()
parser.add_argument('-s','--source', action ='store')
parser.add_argument('--submit', action = 'store_true')
parser.add_argument('--overwrite', action ='store_true')
parser.add_argument('--fullrate', action ='store_true')
pargs = parser.parse_args()

assert(pargs.source is not None)

job_root = 'Poly9_noW201'

all_obs = sorted(glob('/spt/data/bolodata/downsampled/'+pargs.source+'/48*'))
#all_obs = sorted(glob('/spt/user/arahlin/scanify/downsampled/'+pargs.source+'/4*'))
#all_obs = sorted(glob('/spt/user/production/calibration/calframe/ra0hdec-57.5/3*.g3'))
good_obs = []
for ind, path in enumerate(all_obs):
    obsid = path.split('/')[-1].replace('.g3','')
    if int(obsid)>48731000 and int(obsid) < 49231800:
        continue
    if not pargs.overwrite:
        test = glob(os.path.join('/spt/user/ddutcher',pargs.source,'psds',job_root+'_'+obsid+'.g3'))
        if len(test)!=0:
            continue
    cal = glob('/spt/user/production/calibration/calframe/'+pargs.source+'/'+obsid+'.g3')
    if len(cal)!=0:
        good_obs.append(obsid)
print good_obs

condor_dir = '/scratch/ddutcher/condor_logs/'+pargs.source
out_root = '/spt/user/ddutcher/'+pargs.source+'/psds/'
script = '/home/ddutcher/code/spt3g_software/scratch/ddutcher/get_obs_psd.py'

test = True
if pargs.submit:
    test = False
    
if pargs.fullrate:
    rate = 'fullrate'
else:
    rate = 'downsampled'

requirements = '''( ( HAS_CVMFS_spt_opensciencegrid_org ) && ( ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName1 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName2 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName3 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName4 ) || ( RCC_Factory == "ciconnect" ) ) ) && ( OSGVO_OS_STRING == "RHEL 7" ) && (GLIDEIN_ResourceName =!= "NPX") && (GLIDEIN_ResourceName =!= "IIT_CE1"))'''
    
for obs in good_obs:
    jobname = job_root+'_'+str(obs)
    data = glob(os.path.join('/spt/data/bolodata/',rate, pargs.source , 
                                       str(obs),'0*.g3'))
    cal = os.path.join('/spt/user/production/calibration/calframe/',pargs.source,
                        str(obs)+'.g3')
    infiles = [cal] + sorted(data)
    args_in = [os.path.basename(dat) for dat in infiles]
    extra_args='--polyfilt'
    args = '{infiles} -o {outfile} {extra}'.format(infiles = ' '.join(args_in), 
                                           outfile = jobname+'.g3',extra=extra_args)

    condor_submit(script, create_only=test, args = [args],
                   log_root = condor_dir, 
                  output_root = out_root,
                  jobname = jobname,
                  input_files= infiles,
                  clustertools_version='py2-v1',
                  aux_input_files=[],
                  output_files=[jobname+'.g3'],
                  requirements = requirements,
                  request_disk=12*core.G3Units.GB,
                  request_memory=2*core.G3Units.GB,
                  grid_proxy='/home/ddutcher/ddutcher_proxy')
