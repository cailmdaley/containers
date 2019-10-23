'''
Script for making qest and qecl spectra from a parameter file.
Create data products in the following order:
 - ivfs
 - mean-field estimates from the ivfs
 - qcl spectra
 - ncl noise spectra
'''
import sys, os, getopt, imp, argparse, copy, glob, subprocess
import hashlib
import numpy as np
import pylab as pl
import pickle as pk
import datetime
t0 = datetime.datetime.now()

## ---- parameters
parser = argparse.ArgumentParser(description='Cache C_l^{pp} objects')
parser.add_argument('parfile', type=str, nargs=1)
parser.add_argument('--dir_name', type=str, action='store') # the lensing directory's name
parser.add_argument('--output_dir', type=str, action='store', default='/spt/user/panz/lens_test1')# the local output copy output to
parser.add_argument('--do_data_ivf', dest='do_data_ivf', action='store_true') # the option to do data inverse variance filtering
parser.add_argument('--do_ivf', dest='do_ivf', action='store_true') # the option to do inverse variance filtering
parser.add_argument('--do_mf', dest='do_mf', action='store_true') # the option to do mean field calculation
parser.add_argument('--do_data_qcl', dest='do_data_qcl', action='store_true') # the option to do qcl for data
parser.add_argument('--do_qcl_ncl', dest='do_qcl_ncl', action='store_true') # the option to do qcl and ncl
parser.add_argument('--start_idx', dest='start_idx', type=int, action='store', default=0)
parser.add_argument('--end_idx', dest='end_idx', type=int, action='store', default=1)
parser.add_argument('--all_qcl',    dest='all_qcl',  action='store_true')  # Calculate all qcl estimators.
parser.add_argument('--mv_only',    dest='mv_only',  action='store_true')  # Calculate qcl estimators for MV only.  e.g., For systematics test
parser.add_argument('--pp_only',    dest='pp_only',  action='store_true')  # Calculate qcl estimators for POL only.
parser.add_argument('--jack',       dest='do_jack',  action='store_true')  # for jackknife tests
#parser.add_argument('--mv_c_only',  dest='mv_c_only',action='store_true')  # Calculate qcl curl estimators for MV only.  e.g., For systematics test
#parser.add_argument('--pp_c_only',  dest='pp_c_only',action='store_true')  # Calculate qcl curl estimators for POL only.
args    = parser.parse_args()
all_qcl = args.all_qcl
mv_only = args.mv_only
pp_only = args.pp_only
do_jack = args.do_jack
#mv_c_only = args.mv_c_only # curl estimator
#pp_c_only = args.mv_c_only # curl estimator

print('unzip the dir')
#subprocess.check_call('tar -zxvf %s.tar.gz'%args.dir_name, shell=True)


from spt3g import lensing as sl
par     = imp.load_source('par', args.parfile[0])


# Specify which sims to use for qcl libraries
mc_sims     = par.mc_sims
mc_sims_mf  = par.mc_sims_mf
mc_sims_var = par.mc_sims_var

# unlensed variance (default to mc_sims_var)
if hasattr(par, 'mc_sims_unl'):
    mc_sims_unl = par.mc_sims_unl
else:
    mc_sims_unl = mc_sims_var

# phi amplitude correction (default to mc_sims_var)
if hasattr(par, 'mc_sims_qcr_mc'):
    mc_sims_qcr_mc = par.mc_sims_qcr_mc
elif hasattr(par, 'mc_sims_tf'):
    print("found OBSOLETE var mc_sims_tf in par file.  change this to mc_sims_qcr_mc.")
    mc_sims_qcr_mc = par.mc_sims_tf
else:
    mc_sims_qcr_mc  = par.mc_sims_var

# N0 bias (default to mc_sims_var)
if hasattr(par, 'mc_sims_n0'):
    mc_sims_n0  = par.mc_sims_n0
else:
    mc_sims_n0  = par.mc_sims_var

# N1 bias (default to mc_sims_var)
if hasattr(par, 'mc_sims_n1'):
    mc_sims_n1  = par.mc_sims_n1
else:
    mc_sims_n1  = par.mc_sims_var

# Evaluate these estimators
qftkeys = par.__dict__.get('qftkeys')
qclkeys = par.__dict__.get('qclkeys')

if all_qcl:
    qftkeys = [ 'Phi_TT', 'Phi_EE', 'Phi_TE_set', 'Phi_TB_set', 'Phi_EB_set', 'Phi_set', 'Phi_pol_set' ]
    qclkeys = list(zip( qftkeys, qftkeys)) + [ # all 15 estimators
        ('Phi_TT', 'Phi_EE'), ('Phi_TT', 'Phi_EB_set'), ('Phi_TT', 'Phi_TE_set'), ('Phi_EE', 'Phi_TE_set'), ('Phi_EE', 'Phi_EB_set'),
        ('Phi_TE_set', 'Phi_EB_set'), ('Phi_TT', 'Phi_TB_set'), ('Phi_EE', 'Phi_TB_set'), ('Phi_EB_set', 'Phi_TB_set'), ('Phi_TE_set', 'Phi_TB_set') ]

if mv_only:
    print("*** calculating mv (key=p) only")
    qftkeys = [ 'Phi_set']
    qclkeys = list(zip( qftkeys, qftkeys))

if pp_only:
    print("*** calculating POL (key=pp) and MV (key=p) only")
    qftkeys = [ 'Phi_pol_set', 'Phi_set']
    qclkeys = list(zip( qftkeys, qftkeys))


# Specify which ivflibs to make
if not hasattr(par, 'ivflibs_mc_sims_n1'):
    par.ivflibs_mc_sims_n1  = [par.cinv150_len_nofg, par.cinv150_len_t2_nofg] # evaluate for idxs in mc_sims_mf


# Helper functions

def get_ivf(xxx_todo_changeme1):
    (idx, ivflibs) = xxx_todo_changeme1
    for ivflib in eval(ivflibs):
        print("get sim ivf for libs "+ivflibs+", idx = %u"%idx)
        ivflib.get_sim_teb(idx)

def get_qcl(xxx_todo_changeme2):
    (idx, qcllib) = xxx_todo_changeme2
    for k1, k2 in qclkeys:
        print("loading qcl = ", k1, " k2 = ", k2, " idx = ", idx)
        eval(qcllib).get_sim_qcl_lm(k1, idx, k2=k2)

def get_ncl(xxx_todo_changeme3):
    (idx, qcllib) = xxx_todo_changeme3
    for k1, k2 in qclkeys:
        print(("loading ncl = ", k1, " k2 = ", k2, " idx = ", idx))
        eval(qcllib).get_sim_ncl_lm(k1, idx, k2=k2)

def get_dat_qcl(xxx_todo_changeme4):
    (qcllib, k1, k2) = xxx_todo_changeme4
    print(("loading dat qcl = ", k1, " k2 = ", k2))
    eval(qcllib).get_dat_qcl_lm(k1, k2=k2)
    print(("loading dat qcr = ", k1, " k2 = ", k2))
    qcr = eval(qcllib).get_qcr_lm(k1, k2=k2)

# -----------------
# ivf
# -----------------



# Data
if args.do_data_ivf:
    for ivflib in par.ivflibs: #dat 
        ivflib.get_dat_teb()
    outfiles=[]
    outfiles= np.concatenate((outfiles, glob.glob(par.bdir+'qest/run08/cinv150_len_t1p1/'+'dat_teb_bar.pk')))
    outfiles= np.concatenate((outfiles, glob.glob(par.bdir+'qest/run08/cinv150_unl_t1p1/'+'dat_teb_bar.pk')))
    outfiles= np.concatenate((outfiles, glob.glob(par.bdir+'qest/run08/cinv150_len_t2p1_nofg/'+'dat_teb_bar.pk')))
    outfiles= np.concatenate((outfiles, glob.glob(par.bdir+'qest/run08/cinv150_len_t1p1_nofg/'+'dat_teb_bar.pk')))
    for outfile in outfiles:
       print(outfile)
       output_dir= args.output_dir
       copy_path= 'gsiftp://ceph-gridftp1.grid.uchicago.edu:2811/cephfs'+ output_dir+outfile.split(args.dir_name)[-1]
       print(copy_path)
       subprocess.check_call('globus-url-copy '+outfile+' '+copy_path, shell=True)
    

# Sims
if args.do_ivf:
   opts = [] # sim
   idxs= np.arange(args.start_idx, args.end_idx+1, 1)
   print(('idxs are', idxs))
   mc_sims= list(set(idxs).intersection(mc_sims))
   mc_sims_n1= list(set(idxs).intersection(mc_sims_n1))
   mc_sims_unl= list(set(idxs).intersection(mc_sims))
   # mc_sims_unl= list(set(idxs).intersection(mc_sims_unl))
   for idx in mc_sims:
       opts.append((idx, "par.ivflibs_mc_sims"))
   for idx in mc_sims_n1:
       opts.append((idx, "par.ivflibs_mc_sims_n1"))
   for idx in mc_sims_unl:
       opts.append((idx, "par.ivflibs_mc_sims_unl"))


   for opt in opts:
       get_ivf(opt)
   t1=t0
   t2 = datetime.datetime.now()
   print("ivfs took "+str(t2-t1))
   #copy data back
   #local folder= par.bdir
   #remote folder= output_dir
   #find the local files and put them in a list
   outfiles=[]
   outfiles= np.concatenate((outfiles, glob.glob(par.bdir+'qest/run08/cinv150_len_t1p1/'+'sim_*_teb_bar.pk')))
   outfiles= np.concatenate((outfiles, glob.glob(par.bdir+'qest/run08/cinv150_unl_t1p1/'+'sim_*_teb_bar.pk')))
   outfiles= np.concatenate((outfiles, glob.glob(par.bdir+'qest/run08/cinv150_len_t2p1_nofg/'+'sim_*_teb_bar.pk')))
   outfiles= np.concatenate((outfiles, glob.glob(par.bdir+'qest/run08/cinv150_len_t1p1_nofg/'+'sim_*_teb_bar.pk')))
   for outfile in outfiles:
      print(outfile)
      output_dir= args.output_dir
      copy_path= 'gsiftp://ceph-gridftp1.grid.uchicago.edu:2811/cephfs'+ output_dir+outfile.split(args.dir_name)[-1]
      print(copy_path)
      subprocess.check_call('globus-url-copy '+outfile+' '+copy_path, shell=True)

# -----------------
# mfs
# -----------------
if args.do_mf:
   for qcllib in par.qcllibs:
      for k in qftkeys:
        for qestlib, idxs in [ (qcllib.qeA, qcllib.mc_sims_mfA), (qcllib.qeB, qcllib.mc_sims_mfB) ]:
            def get_qft(xxx_todo_changeme):
                (idx) = xxx_todo_changeme
                return qestlib.get_sim_qft(k, idx)
            idxs_input= np.arange(args.start_idx, args.end_idx+1, 1)
            if idxs is not None:
                idxs_input= list(set(list(idxs_input)).intersection(list(idxs)))
            else:
                continue
            if len(idxs_input)==0:
                continue    
            for idx in idxs_input:
                tfname = qestlib.lib_dir.format(prefix="temp") + "/sim_qft_mf_%s_%s_%04d.pk" % (k, hashlib.sha1(np.ascontiguousarray(idxs)).hexdigest(), idx)
                qft_idx= get_qft(idx)
                pk.dump(qft_idx, open(tfname, 'wb'))    
   outfiles=[]
   outfiles= np.concatenate((outfiles, glob.glob(par.bdir+'qest/run08/par_mf100_lx450_lmax3000/qest_len_dd/*mf*')))
   outfiles= np.concatenate((outfiles, glob.glob(par.bdir+'qest/run08/par_mf100_lx450_lmax3000/qest_unl_uu/*mf*')))
   for outfile in outfiles:
      print(outfile)
      output_dir= args.output_dir
      copy_path= 'gsiftp://ceph-gridftp1.grid.uchicago.edu:2811/cephfs'+ output_dir+outfile.split(args.dir_name)[-1]
      print(copy_path)
      subprocess.check_call('globus-url-copy '+outfile+' '+copy_path, shell=True)


# ----------------
# qcl and ncl
# ----------------

# data
if args.do_data_qcl:
   opts = []
   for k1, k2 in qclkeys:
       opts.append(("par.qecl_len_dd", k1, k2))

   for opt in opts:
       get_dat_qcl(opt)

   outfiles=[]
   outfiles= np.concatenate((outfiles, glob.glob(par.bdir+'qest/run08/par_mf100_lx450_lmax3000/qecl_len_dd/cache_lm*pk')))
   for outfile in outfiles:
      print(outfile)
      output_dir= args.output_dir
      copy_path= 'gsiftp://ceph-gridftp1.grid.uchicago.edu:2811/cephfs'+ output_dir+outfile.split(args.dir_name)[-1]
      print(copy_path)
      subprocess.check_call('globus-url-copy '+outfile+' '+copy_path, shell=True)


# sim
if args.do_qcl_ncl:
   idxs = [109, 117, 124, 142, 148, 156, 167, 205, 211, 248, 261, 262, 269, 316, 317, 321, 325, 329, 335, 339, 341, 349, 351, 352, 361, 363, 366, 370, 377, 378, 384, 412, 416, 419, 420, 429, 438, 440, 448, 449, 454, 458, 461, 462, 463, 465, 478, 479, 482, 489, 490]
   idxs = idxs[45:]
   #idxs= np.arange(args.start_idx, args.end_idx+1, 1)
   mc_sims_var=list(set(list(mc_sims_var)).intersection(list(idxs)))
   mc_sims_n0=list(set(list(mc_sims_n0)).intersection(list(idxs)))
   mc_sims_n1=list(set(list(mc_sims_n1)).intersection(list(idxs)))
   mc_sims_qcr_mc=list(set(list(mc_sims_qcr_mc)).intersection(list(idxs)))
   mc_sims_unl=list(set(list(mc_sims_unl)).intersection(list(idxs)))

   # qcl
   opts = []
   for idx in mc_sims_var:               # for error bars
       opts.append((idx, "par.qecl_len_dd"))

   for idx in mc_sims_n0:                # for n0 bias term
       opts.append((idx, "par.qecl_len_ds"))
       opts.append((idx, "par.qecl_len_ss"))

   for idx in mc_sims_n1:                    # for n1 bias term
       opts.append((idx, "par.qecl_len_ss_nofg"))
       opts.append((idx, "par.qecl_len_ss2_nofg"))

   for idx in mc_sims_qcr_mc:                # for error bars
       opts.append((idx, "par.qecl_len_dk"))

   for idx in mc_sims_unl:               # for unlensed sims
       opts.append((idx, "par.qecl_len_uu"))

   opts.append((-1, "par.qecl_len_dd"))      # for data!

   # for jackknife tests:

   if do_jack:
       opts = []
       for idx in mc_sims_var:               # for error bars
           opts.append((idx, "par.qecl_len_dd"))
       for idx in mc_sims_n0:                # for n0 bias term
           opts.append((idx, "par.qecl_len_ds"))
           opts.append((idx, "par.qecl_len_ss"))
       opts.append((-1, "par.qecl_len_dd"))      # for data!

   for opt in opts:
       get_qcl(opt)

   # ncl
   opts = []
   for idx in mc_sims_var:
       opts.append((idx, "par.qecl_len_dd"))
   opts.append((-1, "par.qecl_len_dd")) # This analyzes the data

   for opt in opts:
        get_ncl(opt)

   outfiles=[]
   outfiles= np.concatenate((outfiles, glob.glob(par.bdir+'qest/run08/par_mf100_lx450_lmax3000/qecl_len_dd/cache_lm*pk')))
   outfiles= np.concatenate((outfiles, glob.glob(par.bdir+'qest/run08/par_mf100_lx450_lmax3000/qecl_len_ds/cache_lm*pk')))
   outfiles= np.concatenate((outfiles, glob.glob(par.bdir+'qest/run08/par_mf100_lx450_lmax3000/qecl_len_ss/cache_lm*pk')))
   outfiles= np.concatenate((outfiles, glob.glob(par.bdir+'qest/run08/par_mf100_lx450_lmax3000/qecl_len_ss2_nofg/cache_lm*pk')))
   outfiles= np.concatenate((outfiles, glob.glob(par.bdir+'qest/run08/par_mf100_lx450_lmax3000/qecl_len_ss_nofg/cache_lm*pk')))
   outfiles= np.concatenate((outfiles, glob.glob(par.bdir+'qest/run08/par_mf100_lx450_lmax3000/qecl_len_uu/cache_lm*pk')))
   outfiles= np.concatenate((outfiles, glob.glob(par.bdir+'qest/run08/par_mf100_lx450_lmax3000/qxcl_len_dd/cache_lm*pk')))
   for outfile in outfiles:
      print(outfile)
      output_dir= args.output_dir
      copy_path= 'gsiftp://ceph-gridftp1.grid.uchicago.edu:2811/cephfs'+ output_dir+outfile.split(args.dir_name)[-1]
      print(copy_path)
      subprocess.check_call('globus-url-copy '+outfile+' '+copy_path, shell=True)
