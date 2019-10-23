from spt3g import core, mapmaker, frbutils
import pickle, copy
import numpy as np
import scipy.interpolate
import scipy.ndimage

inj_fluences = [8,11.6,16,32,64,128,256]
inj_fluence_tags = map(lambda v: str(v).replace('.','p'), inj_fluences)

info_tag = 'S3p0_2p0W4p5H0p8Nv0p0Xv1p5E-04P1p1E-03Ng0Gr2Co6p5Fs1'
#info_tag = 'S3p0_2p0W4p0H1p0Nv0p0Xv1p5E-04P1p1E-03Ng0Gr2Co6p5Fs1' #type 1

pf_fn_lst = map(lambda f: '/home/nlharr/tmp/ldfs/pfs_new/pf_phinal_phinal_new_1ms_G98_%s_ldfs_%s.pkl'%(f, info_tag), inj_fluence_tags) 
sat_filter_fn = '/home/nlharr/tmp/ldfs/full_sat_filtered_count_map_%s.pkl' %info_tag
ldf_fn = '/home/nlharr/tmp/ldfs/ldfs_phinal_ldfs_%s.pkl'%info_tag
sigs_to_use = range(7,14)
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


pf_grid = np.zeros((len(sigs_to_use), len(inj_fluences)))

for i, s in enumerate(sigs_to_use):
    pf_grid[i,:] = pf_usable[s]


intfunc = scipy.interpolate.interp2d( inj_fluences, sigs_to_use, pf_grid,
                                      #kind = 'linear')
                           kind = 'cubic')



calc_fluences = np.arange(4,32)
out_pf_grid = np.zeros((len(sigs_to_use), len(calc_fluences)))

for i, s in enumerate(sigs_to_use):
    for j, f in enumerate(calc_fluences):
        out_pf_grid[i,j] = intfunc(f,s)


ref_val = intfunc(32,7)

max_og = np.max(out_pf_grid)

import matplotlib.pyplot as plt
plt.imshow(out_pf_grid, interpolation = 'nearest', cmap = 'Greys_r', extent=[3,32,16,7])


plt.contour(out_pf_grid, levels=[ref_val/3.0], extent=[3,32,7,16],
           colors='r', linewidths=2, linestyles='dashed')



plt.xlabel('Signal Fluence (Jy ms)')
plt.ylabel('Significance Cutoff')
plt.title('Fraction of Events Found for 1 ms FWHM Signal')
plt.show()
