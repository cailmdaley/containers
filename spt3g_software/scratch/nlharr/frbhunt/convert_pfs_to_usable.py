from spt3g import core, mapmaker, frbutils
import pickle, copy
import numpy as np

inj_fluences = [ 8, 11.6, 16, 32, 64, 128, 256]
inj_fluence_tags = map(lambda v: str(v).replace('.','p'), inj_fluences)

#info_tag = 'S3p0_2p0W4p0H0p5Nv0p0Xv1p5E-04P1p1E-03Ng0Gr2Co6p5Fs1' #type 0
#info_tag = 'S3p0_2p0W4p0H1p0Nv0p0Xv1p5E-04P1p1E-03Ng0Gr2Co6p5Fs1' #type 1
#pf_fn_lst = map(lambda f: '/home/nlharr/tmp/ldfs/pfs/pf_phinal_%s_ldfs_%s.pkl'%(f, info_tag), inj_fluence_tags) 
#sat_filter_fn = '/home/nlharr/tmp/ldfs/full_sat_filtered_count_map_%s.pkl' %info_tag
#ldf_fn = '/home/nlharr/tmp/ldfs/ldfs_phinal_ldfs_%s.pkl'%info_tag


#info_tag = 'S3p0_2p0W4p0H0p8Nv0p0Xv1p5E-04P1p1E-03Ng0Gr2Co6p5Fs1'
#pf_fn_lst = map(lambda f: '/home/nlharr/tmp/ldfs/pfs/pf_phinal_NEW_ALL_SQ_GLTCH_G99_p80_%s_ldfs_%s.pkl'%(f, info_tag), inj_fluence_tags) 
#sat_filter_fn = '/home/nlharr/tmp/ldfs/full_sat_filtered_count_map_NEW_ALL_SQ_thinner_G99_p80_ldfs_%s.pkl' %info_tag
#ldf_fn = '/home/nlharr/tmp/ldfs/ldfs_phinal_NEW_ALL_SQ_GLTCH_G99_p80_ldfs_%s.pkl'%info_tag




info_tag = 'S3p0_2p0W4p5H0p8Nv0p0Xv1p5E-04P1p1E-03Ng0Gr2Co6p5Fs1'
#pf_fn_lst = map(lambda f: '/home/nlharr/tmp/ldfs/pfs/pf_phinal_phinal_v2_G98_%s_ldfs_%s.pkl'%(f, info_tag), inj_fluence_tags) 


#pf_fn_lst = map(lambda f: '/home/nlharr/tmp/ldfs/pfs/pf_phinal_phinal_v2_5ms_G98_%s_ldfs_%s.pkl'%(f, info_tag), inj_fluence_tags) 

pf_fn_lst = map(lambda f: '/home/nlharr/tmp/ldfs/pfs_new/pf_phinal_phinal_new_1ms_G98_%s_ldfs_%s.pkl'%(f, info_tag), inj_fluence_tags) 




sat_filter_fn = '/home/nlharr/tmp/ldfs/full_sat_filtered_thinner_dfs_phinal_phinal_v2_G98_ldfs_%s.pkl' %info_tag
ldf_fn = '/home/nlharr/tmp/ldfs/ldfs_phinal_phinal_v2_G98_ldfs_%s.pkl'%info_tag


#info_tag = 'S4p0_2p5W4p0H0p8Nv0p0Xv1p5E-04P1p0E-03Ng0Gr2Co6p5Fs1'
#pf_fn_lst = map(lambda f: '/home/nlharr/tmp/ldfs/pfs/pf_phinal_NEW_ALL_SQ_GLTCH_G99_1p0_%s_ldfs_%s.pkl'%(f, info_tag), inj_fluence_tags) 
#sat_filter_fn = '/home/nlharr/tmp/ldfs/full_sat_filtered_count_map_NEW_ALL_SQ_thinner_G99_1p0_ldfs_%s.pkl' %info_tag
#ldf_fn = '/home/nlharr/tmp/ldfs/ldfs_phinal_NEW_ALL_SQ_GLTCH_G99_1p0_ldfs_%s.pkl'%info_tag


#59403.0030908271
#60881.916030847686

sigs_to_use = [7, 9, 13]
#sigs_to_use = [8, 9, 13]

#sigs_to_use = [8]
#sigs_to_use = [13]
#sigs_to_use = [13]
#sigs_to_use = [7,9,13]
exposure_ldf_index = 0

##############analysis code

total_pf = {}

for inj_fluence, pf_fn in zip(inj_fluences, pf_fn_lst):
    tmp = pickle.load(open(pf_fn))
    sigs = tmp[0]
    ldfs = tmp[1]
    for sig, ldf in zip(sigs, ldfs):
        total_pf[(sig,inj_fluence)] = float(ldf.n_found)/ldf.n_scans

pf_usable = {}
for sig in sigs:
    pf_usable[sig] = []
    for flu in inj_fluences:
        pf_usable[sig].append(total_pf[(sig,flu)])

for s in pf_usable.keys():
    pf_usable[s] = np.array(pf_usable[s])
counts = pickle.load(open(sat_filter_fn))

tmp = pickle.load(open(ldf_fn))
significances = tmp[0]
real_ldfs = tmp[1]

exposure_ldf = real_ldfs[exposure_ldf_index]

sigs_to_use = sorted(sigs_to_use)
sig_indices = map(significances.index, sigs_to_use)

#counts, pfs, ldfs
fcounts = {}
fpf_usable = {}
fldfs = {}

for i in range(len(sigs_to_use)):
    fldfs[sigs_to_use[i]] = copy.copy(real_ldfs[sig_indices[i]])

for i in range(len(sigs_to_use)-1):
    fcounts[sigs_to_use[i]] = counts[sigs_to_use[i]] - counts[sigs_to_use[i+1]]
    fpf_usable[sigs_to_use[i]] = pf_usable[sigs_to_use[i]]-pf_usable[sigs_to_use[i+1]]
    fldfs[sigs_to_use[i]] -= fldfs[sigs_to_use[i+1]]
fcounts[sigs_to_use[-1]] = counts[sigs_to_use[-1]]
fpf_usable[sigs_to_use[-1]] = pf_usable[sigs_to_use[-1]]

print(fcounts)
print(fpf_usable)
print(fldfs)

pickle.dump( ( sigs_to_use, inj_fluences, fcounts, fpf_usable, exposure_ldf, fldfs), 
             open("bg_analysis_bundle.pkl", 'w'))
