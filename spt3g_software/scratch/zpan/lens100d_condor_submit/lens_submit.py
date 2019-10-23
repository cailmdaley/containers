import sys, os, getopt, imp, argparse, copy, subprocess
import hashlib#, joblib
import numpy as np
import pylab as pl
import pickle as pk
import datetime
from spt3g.cluster.condor_tools import condor_submit
from spt3g import core
import time
from spt3g import lensing as sl

## ---- parameters
parser = argparse.ArgumentParser(description='Submit condor jobs for lensing')
parser.add_argument('parfile', type=str, nargs=1)
parser.add_argument('--condor_folder',   dest='condor_folder', type=str, action='store', default='/spt/user/panz/lens_condor_new_pipeline/')
parser.add_argument('--do_ivf',    dest='do_ivf',  action='store_true')  # Do inverse-variance filtering for data and sims.
parser.add_argument('--do_mf',    dest='do_mf',  action='store_true')  # Do mean field estimation, and save individual mean field outputs. 
parser.add_argument('--do_avg_mf',    dest='do_avg_mf',  action='store_true')  # Average the individual mean field outputs. 
parser.add_argument('--do_qcl_ncl',       dest='do_qcl_ncl',  action='store_true')  # Do the qcl and ncl estimations.


args    = parser.parse_args()
par     = imp.load_source('par', args.parfile[0])
condor_folder= args.condor_folder
do_ivf= args.do_ivf
do_mf=args.do_mf
do_avg_mf= args.do_avg_mf
do_qcl_ncl=args.do_qcl_ncl



bdir= par.bdir
condor_folder='/spt/user/panz/lens_condor_new_pipeline/'
n_idx= 1
requirements = '''( ( HAS_CVMFS_spt_opensciencegrid_org ) && ( ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName1 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName2 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName3 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName4 ) || ( RCC_Factory == "ciconnect" ) ) ) && ( OSGVO_OS_STRING == "RHEL 7" || ( RCC_Factory == "ciconnect" ) ) && (GLIDEIN_ResourceName =!= "NPX") )'''
script= '/home/panz/lensing100d_new_pipeline/make_qecl_condor.py'
test=False
clear_job=False
log_root= '/scratch/panz/lensing_logs'
qcl_to_run= '--all_qcl'
user_code=''


def mkdir_data_ivf(bdir, condor_folder, folder_name):
    data_ivf_folder= os.path.join(condor_folder, folder_name)
    if not os.path.exists(data_ivf_folder):
       subprocess.check_call('mkdir -p %s'%data_ivf_folder, shell=True)
    subprocess.check_call('find %s'%bdir+' -mindepth 1 -depth -type d -printf "%P\n"| while read dir; do mkdir -p'+' "%s/$dir"; done; '%data_ivf_folder, shell=True)
    filenames=['inputs', 'scripts/*py', 'scripts/*sh', 'data/map/20140418_ra23h30dec55_lens100d/halfs/150/half_A_0.g3',
               'data/map/20140418_ra23h30dec55_lens100d/halfs/150/half_B_0.g3', 'libs/sims_cmb/t1p1_scan/sim_0000_buf.npy',
               'libs/sims_cmb/t1p1_scan/weight.npy', 'libs/sims_cmb/t2p1_scan/sim_0000_buf.npy', 'libs/sims_cmb/t2p1_scan/weight.npy'   ]
    for filename in filenames:
     try:
       subprocess.check_call('cp -rf '+os.path.join(bdir, filename)+' '+os.path.dirname(os.path.join(data_ivf_folder, filename)), shell=True)
     except:
       continue
    subprocess.call('cd %s; tar --warning=no-file-changed -zcf %s.tar.gz %s;'%(condor_folder, folder_name, folder_name), shell=True)
    #print "command"+'tar -zcf %s.tar.gz %s'%(os.path.join(condor_folder,folder_name), data_ivf_folder)
    #subprocess.check_call('tar --warning=no-file-changed   -zcf %s.tar.gz --directory %s %s'%( data_ivf_folder, condor_folder, folder_name), shell=True)
    return '%s.tar.gz'%(os.path.join(condor_folder,folder_name))


def mkdir_ivf(bdir, condor_folder, folder_name, idxs):
    data_ivf_folder= os.path.join(condor_folder, folder_name)
    if not os.path.exists(data_ivf_folder):
       subprocess.check_call('mkdir -p %s'%data_ivf_folder, shell=True)
    subprocess.check_call('find %s'%bdir+' -mindepth 1 -depth -type d -printf "%P\n"| while read dir; do mkdir -p'+' "%s/$dir"; done; '%data_ivf_folder, shell=True)
    filenames=['inputs', 'scripts/*py', 'scripts/*sh', 'data/map/20140418_ra23h30dec55_lens100d/halfs/150/half_A_0.g3',
               'data/map/20140418_ra23h30dec55_lens100d/halfs/150/half_B_0.g3', 'libs/sims_cmb/t1p1_scan/weight.npy',
               'libs/sims_cmb/t2p1_scan/weight.npy', 'libs/sims_cmb/t1_unl_scan/weight.npy']
    for idx in idxs:
        filenames.append('libs/sims_cmb/t1p1_scan/sim_%04d_buf.npy'%idx)
        filenames.append('libs/sims_cmb/t2p1_scan/sim_%04d_buf.npy'%idx)
        filenames.append('libs/sims_cmb/t1_unl_scan/sim_%04d_buf.npy'%idx)
        filenames.append('data/map/20140418_ra23h30dec55_lens100d/halfs/150/half_A_%d.g3'%idx)
        filenames.append('data/map/20140418_ra23h30dec55_lens100d/halfs/150/half_B_%d.g3'%idx)
    for filename in filenames:
     try:
       subprocess.check_call('cp -rf '+os.path.join(bdir, filename)+' '+os.path.dirname(os.path.join(data_ivf_folder, filename)), shell=True)
     except:
       continue
    subprocess.call('cd %s; tar --warning=no-file-changed -zcf %s.tar.gz %s;'%(condor_folder, folder_name, folder_name), shell=True)
    #print "command"+'tar -zcf %s.tar.gz %s'%(os.path.join(condor_folder,folder_name), data_ivf_folder)
    #subprocess.check_call('tar --warning=no-file-changed   -zcf %s.tar.gz --directory %s %s'%( data_ivf_folder, condor_folder, folder_name), shell=True)
    return '%s.tar.gz'%(os.path.join(condor_folder,folder_name))


def mkdir_mf(bdir, condor_folder, folder_name, idxs):
    data_ivf_folder= os.path.join(condor_folder, folder_name)
    if not os.path.exists(data_ivf_folder):
       subprocess.check_call('mkdir -p %s'%data_ivf_folder, shell=True)
    subprocess.check_call('find %s'%bdir+' -mindepth 1 -depth -type d -printf "%P\n"| while read dir; do mkdir -p'+' "%s/$dir"; done; '%data_ivf_folder, shell=True)
    filenames=['inputs', 'scripts/*py', 'scripts/*sh', 'data/map/20140418_ra23h30dec55_lens100d/halfs/150/half_A_0.g3',
               'data/map/20140418_ra23h30dec55_lens100d/halfs/150/half_B_0.g3']
    for idx in idxs:
        filenames.append('qest/run08/cinv150_len_t1p1/sim_%04d_teb_bar.pk'%idx)
        filenames.append('qest/run08/cinv150_len_t1p1_nofg/sim_%04d_teb_bar.pk'%idx)
        filenames.append('qest/run08/cinv150_len_t2p1_nofg/sim_%04d_teb_bar.pk'%idx)
        filenames.append('qest/run08/cinv150_unl_t1p1/sim_%04d_teb_bar.pk'%idx)
    for filename in filenames:
     try:
       subprocess.check_call('cp -rf '+os.path.join(bdir, filename)+' '+os.path.dirname(os.path.join(data_ivf_folder, filename)), shell=True)
     except:
       continue
    subprocess.call('cd %s; tar --warning=no-file-changed -zcf %s.tar.gz %s;'%(condor_folder, folder_name, folder_name), shell=True)
    #print "command"+'tar -zcf %s.tar.gz %s'%(os.path.join(condor_folder,folder_name), data_ivf_folder)
    #subprocess.check_call('tar --warning=no-file-changed   -zcf %s.tar.gz --directory %s %s'%( data_ivf_folder, condor_folder, folder_name), shell=True)
    return '%s.tar.gz'%(os.path.join(condor_folder,folder_name))

               
def mkdir_qcl_ncl(bdir, condor_folder, folder_name, idxs):
    data_ivf_folder= os.path.join(condor_folder, folder_name)
    if not os.path.exists(data_ivf_folder):
       subprocess.check_call('mkdir -p %s'%data_ivf_folder, shell=True)
    subprocess.check_call('find %s'%bdir+' -mindepth 1 -depth -type d -printf "%P\n"| while read dir; do mkdir -p'+' "%s/$dir"; done; '%data_ivf_folder, shell=True)
    filenames=['inputs', 'scripts/*py', 'scripts/*sh', 'data/map/20140418_ra23h30dec55_lens100d/halfs/150/half_A_0.g3',
               'data/map/20140418_ra23h30dec55_lens100d/halfs/150/half_B_0.g3', 'qest/run08/cinv150_len_t1p1/dat_teb_bar.pk',
               'qest/run08/cinv150_len_t1p1_nofg/dat_teb_bar.pk', 'qest/run08/cinv150_len_t2p1_nofg/dat_teb_bar.pk',
               'qest/run08/cinv150_unl_t1p1/dat_teb_bar.pk', 'qest/run08/par_mf100_lx450_lmax3000/qest_len_dd/*mf*',
               'qest/run08/par_mf100_lx450_lmax3000/qest_unl_uu/*mf*'   ]
    def roll_i( i):
        return i - np.mod(i, 2) + np.mod(i + 1,2)
    for idx in idxs:
        idx2= roll_i(idx)
        filenames.append('qest/run08/cinv150_len_t1p1/sim_%04d_teb_bar.pk'%idx)
        filenames.append('qest/run08/cinv150_len_t1p1_nofg/sim_%04d_teb_bar.pk'%idx)
        filenames.append('qest/run08/cinv150_len_t2p1_nofg/sim_%04d_teb_bar.pk'%idx)
        filenames.append('qest/run08/cinv150_unl_t1p1/sim_%04d_teb_bar.pk'%idx)
        filenames.append('qest/run08/cinv150_len_t1p1/sim_%04d_teb_bar.pk'%idx2)
        filenames.append('qest/run08/cinv150_len_t1p1_nofg/sim_%04d_teb_bar.pk'%idx2)
        filenames.append('qest/run08/cinv150_len_t2p1_nofg/sim_%04d_teb_bar.pk'%idx2)
        filenames.append('qest/run08/cinv150_unl_t1p1/sim_%04d_teb_bar.pk'%idx2)
        filenames.append('libs/sims_cmb/t1p1_proj/sim_%04d_*g3'%idx)

    for filename in filenames:
     try:
       subprocess.check_call('cp -rf '+os.path.join(bdir, filename)+' '+os.path.dirname(os.path.join(data_ivf_folder, filename)), shell=True)
     except:
       continue
    subprocess.call('cd %s; tar --warning=no-file-changed -zcf %s.tar.gz %s;'%(condor_folder, folder_name, folder_name), shell=True)
    #print "command"+'tar -zcf %s.tar.gz %s'%(os.path.join(condor_folder,folder_name), data_ivf_folder)
    #subprocess.check_call('tar --warning=no-file-changed   -zcf %s.tar.gz --directory %s %s'%( data_ivf_folder, condor_folder, folder_name), shell=True)
    return '%s.tar.gz'%(os.path.join(condor_folder,folder_name))


def mkdir_data_qcl(bdir, condor_folder, folder_name):
    data_ivf_folder= os.path.join(condor_folder, folder_name)
    if not os.path.exists(data_ivf_folder):
       subprocess.check_call('mkdir -p %s'%data_ivf_folder, shell=True)
    subprocess.check_call('find %s'%bdir+' -mindepth 1 -depth -type d -printf "%P\n"| while read dir; do mkdir -p'+' "%s/$dir"; done; '%data_ivf_folder, shell=True)
    filenames=['inputs', 'scripts/*py', 'scripts/*sh', 'data/map/20140418_ra23h30dec55_lens100d/halfs/150/half_A_0.g3',
               'data/map/20140418_ra23h30dec55_lens100d/halfs/150/half_B_0.g3', 'qest/run08/cinv150_len_t1p1/dat_teb_bar.pk',
               'qest/run08/cinv150_len_t1p1_nofg/dat_teb_bar.pk', 'qest/run08/cinv150_len_t2p1_nofg/dat_teb_bar.pk',
               'qest/run08/cinv150_unl_t1p1/dat_teb_bar.pk', 'qest/run08/par_mf100_lx450_lmax3000/qest_len_dd/*mf*',
               'qest/run08/par_mf100_lx450_lmax3000/qest_unl_uu/*mf*'   ]

    for filename in filenames:
     try:
       subprocess.check_call('cp -rf '+os.path.join(bdir, filename)+' '+os.path.dirname(os.path.join(data_ivf_folder, filename)), shell=True)
     except:
       continue
    subprocess.call('cd %s; tar --warning=no-file-changed -zcf %s.tar.gz %s;'%(condor_folder, folder_name, folder_name), shell=True)
    #print "command"+'tar -zcf %s.tar.gz %s'%(os.path.join(condor_folder,folder_name), data_ivf_folder)
    #subprocess.check_call('tar --warning=no-file-changed   -zcf %s.tar.gz --directory %s %s'%( data_ivf_folder, condor_folder, folder_name), shell=True)
    return '%s.tar.gz'%(os.path.join(condor_folder,folder_name))




if do_ivf:
    #For data
    dir_name= 'lensing_data_ivf'
    if not os.path.exists(os.path.join(condor_folder, dir_name+'.tar.gz')):
       dir_tar=mkdir_data_ivf(bdir, condor_folder, dir_name)
       file_exists=0
    else:
       dir_tar='%s.tar.gz'%(os.path.join(condor_folder,dir_name))
       file_exists=1
    infiles= [dir_tar]
    arguments= '$PWD/%s/scripts/par_run.py --dir_name %s --output_dir %s --do_data_ivf'%(dir_name, dir_name, bdir)
    if (not clear_job) or (clear_job and file_exists==0):
      condor_submit(script, create_only=test, args = [arguments],
                    log_root = log_root,
                    output_root = bdir,
                    jobname = dir_name,
                    input_files= infiles,
                    aux_input_files=[],
                    output_files=[],
                    user_code=user_code,
                    requirements = requirements,
                    request_disk=15*core.G3Units.GB,
                    request_memory=2*core.G3Units.GB,
                    grid_proxy='/home/panz/panz_proxy')
    #For sims
    idx_all=list(set(np.concatenate((par.mc_sims, par.mc_sims_n1, par.mc_sims_unl))))
    for i in range(len(idx_all)//n_idx+1):
      idxs= idx_all[min(i*n_idx, len(idx_all)):min(i*n_idx+n_idx, len(idx_all))]
      if len(idxs)>0:      
       dir_name= 'lensing_ivf_%d_%d'%(idxs[0], idxs[-1])
       if not os.path.exists(os.path.join(condor_folder, dir_name+'.tar.gz')):
         dir_tar=mkdir_ivf(bdir, condor_folder, dir_name, idxs)
         file_exists=0
       else: 
         dir_tar='%s.tar.gz'%(os.path.join(condor_folder,dir_name))
         file_exists=1
       infiles= [dir_tar]
       arguments= '$PWD/%s/scripts/par_run.py --dir_name %s --output_dir %s --do_ivf --start_idx %d --end_idx %d'%(dir_name, dir_name, bdir, idxs[0], idxs[-1])
       # create lensing folder for sim
       # data and files needed, lensing folder, script for running
       # condor_submit module, the arguments of the file need to have output folder, or par.bdir
       if (not clear_job) or (clear_job and file_exists==0):
         condor_submit(script, create_only=test, args = [arguments],
                       log_root = log_root,
                       output_root = bdir,
                       jobname = dir_name,
                       input_files= infiles,
                       aux_input_files=[],
                       output_files=[],
                       user_code=user_code,
                       requirements = requirements,
                       request_disk=15*core.G3Units.GB,
                       request_memory=2*core.G3Units.GB,
                       grid_proxy='/home/panz/panz_proxy')



if do_mf:
    idx_all=[]
    for qcllib in par.qcllibs:
        for qestlib, idxs in [ (qcllib.qeA, qcllib.mc_sims_mfA), (qcllib.qeB, qcllib.mc_sims_mfB) ]:
            if not (idxs is None):
              idx_all=np.concatenate((idx_all, idxs))
    idx_all= list(set(idx_all))
    for i in range(len(idx_all)//n_idx+1):
      idxs= idx_all[min(i*n_idx, len(idx_all)):min(i*n_idx+n_idx, len(idx_all))]
      if len(idxs)>0:      
       dir_name= 'lensing_mf_%d_%d'%(idxs[0], idxs[-1])
       if not os.path.exists(os.path.join(condor_folder, dir_name+'.tar.gz')):
         dir_tar=mkdir_mf(bdir, condor_folder, dir_name, idxs)
         file_exists=0
       else: 
         dir_tar='%s.tar.gz'%(os.path.join(condor_folder,dir_name))
         file_exists=1
       infiles= [dir_tar]
       arguments= '$PWD/%s/scripts/par_run.py --dir_name %s --output_dir %s --do_mf --start_idx %d --end_idx %d %s'%(dir_name, dir_name, bdir, idxs[0], idxs[-1], qcl_to_run)
       # create lensing folder for sim
       # data and files needed, lensing folder, script for running
       # condor_submit module, the arguments of the file need to have output folder, or par.bdir
       if (not clear_job) or (clear_job and file_exists==0):
         condor_submit(script, create_only=test, args = [arguments],
                       log_root = log_root,
                       output_root = bdir,
                       jobname = dir_name,
                       input_files= infiles,
                       aux_input_files=[],
                       output_files=[],
                       user_code=user_code,
                       requirements = requirements,
                       request_disk=15*core.G3Units.GB,
                       request_memory=2*core.G3Units.GB,
                       grid_proxy='/home/panz/panz_proxy')


if do_avg_mf:
    qftkeys = [ 'Phi_TT', 'Phi_EE', 'Phi_TE_set', 'Phi_TB_set', 'Phi_EB_set', 'Phi_set', 'Phi_pol_set' ]
    for qcllib in par.qcllibs:
       for k in qftkeys:
          for qestlib, idxs in [ (qcllib.qeA, qcllib.mc_sims_mfA), (qcllib.qeB, qcllib.mc_sims_mfB) ]:
             if idxs is None:
                continue
             tfname = qestlib.lib_dir.format(prefix="temp") + "/sim_qft_mf_%s_%s.pk" % (k, hashlib.sha1(np.ascontiguousarray(idxs)).hexdigest())
             if 1:#not os.path.exists(tfname):
                print(("caching mf:", tfname))
                qft_mf_avg = sl.utils.AverageObjects()
                def get_qft(idx):
                    fname= qestlib.lib_dir.format(prefix="temp") + "/sim_qft_mf_%s_%s_%04d.pk" % (k, hashlib.sha1(np.ascontiguousarray(idxs)).hexdigest(), idx)
                    qft= pk.load(open(fname, 'rb'))
                    subprocess.check_call('rm -rf %s'%fname, shell=True)
                    return qft
                for idx in idxs:
                      qft_mf_avg.add(get_qft(idx) )
                pk.dump( qft_mf_avg.get(), open(tfname, 'wb') )


if do_qcl_ncl:
    #For data
    dir_name= 'lensing_data_qcl'
    if not os.path.exists(os.path.join(condor_folder, dir_name+'.tar.gz')):
       dir_tar=mkdir_data_qcl(bdir, condor_folder, dir_name)
       file_exists=0
    else:
       dir_tar='%s.tar.gz'%(os.path.join(condor_folder,dir_name))
       file_exists=1
    infiles= [dir_tar]
    arguments= '$PWD/%s/scripts/par_run.py --dir_name %s --output_dir %s --do_data_qcl %s'%(dir_name, dir_name, bdir, qcl_to_run)
    if (not clear_job) or (clear_job and file_exists==0):
      condor_submit(script, create_only=test, args = [arguments],
                    log_root = log_root,
                    output_root = bdir,
                    jobname = dir_name,
                    input_files= infiles,
                    aux_input_files=[],
                    output_files=[],
                    user_code=user_code,
                    requirements = requirements,
                    request_disk=15*core.G3Units.GB,
                    request_memory=2*core.G3Units.GB,
                    grid_proxy='/home/panz/panz_proxy')
    #For sims
    idx_all=list(set(np.concatenate((par.mc_sims_var, par.mc_sims_n0, par.mc_sims_n1, par.mc_sims_qcr_mc, par.mc_sims_unl))))
    for i in range(len(idx_all)//n_idx+1):
      idxs= idx_all[min(i*n_idx, len(idx_all)):min(i*n_idx+n_idx, len(idx_all))]
      if len(idxs)>0:      
       dir_name= 'lensing_qcl_ncl_%d_%d'%(idxs[0], idxs[-1])
       if not os.path.exists(os.path.join(condor_folder, dir_name+'.tar.gz')):
         dir_tar=mkdir_qcl_ncl(bdir, condor_folder, dir_name, idxs)
         file_exists=0
       else:
         dir_tar='%s.tar.gz'%(os.path.join(condor_folder,dir_name))
         file_exists=1
       infiles= [dir_tar]
       arguments= '$PWD/%s/scripts/par_run.py --dir_name %s --output_dir %s --do_qcl_ncl --start_idx %d --end_idx %d %s'%(dir_name, dir_name, bdir, idxs[0], idxs[-1], qcl_to_run)
       # create lensing folder for sim
       # data and files needed, lensing folder, script for running
       # condor_submit module, the arguments of the file need to have output folder, or par.bdir
       if (not clear_job) or (clear_job and file_exists==0):
         condor_submit(script, create_only=test, args = [arguments],
                       log_root = log_root,
                       output_root = bdir,
                       jobname = dir_name,
                       input_files= infiles,
                       aux_input_files=[],
                       output_files=[],
                       user_code=user_code,
                       requirements = requirements,
                       request_disk=15*core.G3Units.GB,
                       request_memory=2*core.G3Units.GB,
                       grid_proxy='/home/panz/panz_proxy')
