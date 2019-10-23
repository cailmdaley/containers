import os
import glob
from spt3g.cluster.condor_tools import condor_submit
import argparse
from spt3g import core
import pickle 

parser = argparse.ArgumentParser()
#parser.add_argument('obsids', nargs = '+', type = str)
parser.add_argument('--submit', action = 'store_true')
parser.add_argument('-s','--source', action ='store', 
                    default='ra0hdec-57.5')

pargs = parser.parse_args()
#obsids = pargs.obsids

#obsids = pickle.load(open('/home/javva/spt3g_software/scratch/javva/condor/elnod_matches.pkl','rb'))

sim_files = glob.glob('/spt/user/arahlin/synfast/*')


job_root = 'sim_input_maps_synfast'
condor_dir = os.path.join('/home/javva/grid/condor_logs', pargs.source, job_root)
out_root = os.path.join('/spt/user/javva/lowell/maps/', pargs.source, job_root)
script = '/home/javva/spt3g_software/scratch/javva/pf_k_space_maps/makemaps_sim_donothing_lpf.py'

test = False
'''
if pargs.submit:
    test = False
'''
requirements = '''( ( HAS_CVMFS_spt_opensciencegrid_org ) && ( ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName1 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName2 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName3 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName4 ) || ( RCC_Factory == "ciconnect" ) ) ) && ( OSGVO_OS_STRING == "RHEL 7" ) && (GLIDEIN_ResourceName =!= "NPX") )'''
baddies = 0
for sim in sim_files:
    try:
        jobname = job_root+'_'+sim.split('sim')[1].split('_')[0]
#        print(jobname)
        obs = '-15997020'
        data = glob.glob('/spt/user/nwhitehorn/sptpol/fullrate/'+pargs.source+'/' +str(obs)+'/0*.g3')
#        print(data)
        cal = '/spt/user/nwhitehorn/sptpol/autoproc/calibration/calframe/'+pargs.source+'/'+str(obs)+'.g3'
        a = sim
        infiles = [cal] + sorted(data)
        infiles_all = [cal]+sorted(data)+[a]
#        print(infiles)
        args_in1 = [os.path.basename(dat) for dat in infiles]
        args_in = ['./'+d for d in args_in1]
        extra_args= '-s'
        args = '{infiles} -v -s {extra_extra} -o {outfile}'.format(infiles = ' '.join(args_in), 
                                                                   outfile = jobname+'.g3', extra_extra = './'+os.path.basename(sim))
        print(args)

        condor_submit(script, create_only=test, args = [args],
                  log_root = os.path.join(condor_dir, str(obs)),
                  output_root = os.path.join(out_root, str(obs)),
                  jobname = jobname,
                  input_files= infiles_all,
                  clustertools_version='py3-v2',
                  aux_input_files=['/home/javva/spt3g_software/scratch/javva/better_gain_matching/ptsrc_config_ra0hdec-57p5_both_50mJy.txt'],
                  output_files=[jobname+'.g3'],
                  requirements = requirements,
                  request_disk=15*core.G3Units.GB,
                  request_memory=2*core.G3Units.GB,
                  grid_proxy='/home/javva/javva_proxy')
    except:
        print('THIS OBSID DIDNT WORK ==',obs)
        baddies = baddies+1

print("end of processing. didn't complete ",baddies)
