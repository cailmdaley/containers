from spt3g import core, coordinateutils, mapmaker, mapspectra
import glob as glob
import numpy, pickle, sys, copy
import numpy as np
import pylab
import scipy.ndimage as ndimage
try:
	FileNotFoundError
except NameError:
	FileNotFoundError = IOError

'''
This code evaluates the power spectra of the input and output sims. And option is to include a k-space mask in the process.
'''
inputs_list = glob.glob('/spt/user/javva/lowell/maps/ra0hdec-57.5/sim_input_maps_sasha_noninterp_lpf/-15997020/*.g3')
outputs_dir = '/spt/user/javva/lowell/maps/ra0hdec-57.5/sims_maps_filtered_nonsense_p4_sasha_nointerp_npf/-15997020/*'
psd_dir = '/spt/user/javva/lowell/2dpsd/bun'#directory to the 2dpsds of the filtered sims
k_weights_sim_mask ='simkspacetf_100_correctapod_median_pluscor.pkl'
do_k_weights = False
bundle_number = 1 #eventually I'll put this in a loop

#some stuff for k-weights mask
stf = pickle.load(open(k_weights_sim_mask,'rb'))
q_ft_sn = stf['q']
u_ft_sn = stf['u']
t_ft_sn = stf['t']

i = bundle_number
f = list(core.G3File(map_dir+'bundles_3gpipe_'+str(i)+'.g3'))[0]
ep = pickle.load(open(psd_dir+str(i)+".pkl","rb"))
k_weights_t = ep['t']
k_weights_q = ep['q']
k_weights_u = ep['u']

        

def ps_from_map(fr, apod_mask, qu=False, coadd=None,ell_weights_2d_t = None, ell_weights_2d_e = None, ell_weights_2d_b = None):
        apod_good, msg  = mapspectra.apodmask.validate_apodization_mask(apod_mask, fr['Wpol'])
        if not apod_good:
                core.log_error(msg)

        t,q,u = mapmaker.mapmakerutils.remove_weight(fr['T'], fr['Q'], fr['U'], fr['Wpol'])
        if coadd is not None:
                t -= coadd[0]
                q -= coadd[1]
                u -= coadd[2]

        qf, uf = mapspectra.basicmaputils.flatten_pol(q,u)

        cen_ells = numpy.linspace(10, 2500, 200)
        ell_bins = mapspectra.basicmaputils.get_reg_spaced_ell_bins(cen_ells
        t_cls,e_cls,b_cls = mapspectra.basicmaputils.get_map_cls(t,qf,uf, apod_mask, ell_bins, qu=qu,ell_weights_2d_t = ell_weights_2d_t, ell_weights_2d_e = ell_weights_2d_e, ell_weights_2d_b = ell_weights_2d_b)

        return cen_ells, t_cls, e_cls, b_cls
apod= pickle.load(open("/spt/user/javva/lowell/good_apodization_mask.pkl","rb"))
numpy.asarray(apod)[numpy.asarray(apod) < 0] = 0

tcls1 = {}
qcls1 = {}
ucls1 = {}
for i in inputs_list:
        num = i.split('/')[-1].split('_')[-1]
        op = glob.glob(outputs_dir+num)
        if len(op)<1:
                continue 
                print('nope')
        
        fi = list(core.G3File(i))[0]
        fo = list(core.G3File(op[0]))[0]

        print(i)
        ci,ti,qi,ui = ps_from_map(fi, numpy.asarray(apod), qu=True)
        if do_k_weights:
            co,to,qo,uo = ps_from_map(fo, numpy.asarray(apod), qu=True,ell_weights_2d_t=(1./t_ft_sn)*(np.abs(k_weights_t**2)),ell_weights_2d_e = (1./q_ft_sn)*(np.abs(k_weights_q**2)),ell_weights_2d_b = (1./u_ft_sn)*(np.abs(k_weights_u**2)), coadd = None)
        else:                                                        
            co,to,qo,uo = ps_from_map(fo, numpy.asarray(apod), qu=True)
        tcls1[num] = {}
        qcls1[num] = {}
        ucls1[num] = {}
        tcls1[num]['ci'] = ci
        qcls1[num]['qi'] = qi**0.5 / (core.G3Units.arcmin * core.G3Units.uK)
        ucls1[num]['ui'] = ui**0.5 / (core.G3Units.arcmin * core.G3Units.uK)
        tcls1[num]['ti'] = ti**0.5 / (core.G3Units.arcmin * core.G3Units.uK)
        tcls1[num]['co'] = co
        qcls1[num]['qo'] = qo**0.5 / (core.G3Units.arcmin * core.G3Units.uK)
        ucls1[num]['uo'] = uo**0.5 / (core.G3Units.arcmin * core.G3Units.uK)
        tcls1[num]['to'] = to**0.5 / (core.G3Units.arcmin * core.G3Units.uK)
        dict_nu = {}
        dict_nu['T'] = tcls1
        dict_nu['Q'] = qcls1
        dict_nu['U']= ucls1
        with open('individual_tfs_lpf.pkl', 'wb') as handle:
                pickle.dump(dict_nu, handle, protocol=pickle.HIGHEST_PROTOCOL)

tfs = {}
avg = []
for m in range(len(dict_nu['Q'][num]['qo'])):
        vals = [(dict_nu['Q'][i]['qo']/dict_nu['Q'][i]['qi'])[m] for i in list(dict_nu['Q'].keys())]
        avg = np.append(avg, np.mean(vals))
tfs['Q'] = avg
avg = []
for m in range(len(dict_nu['U'][num]['uo'])):
        vals = [(dict_nu['U'][i]['uo']/dict_nu['U'][i]['ui'])[m] for i in list(dict_nu['U'].keys())]
        avg = np.append(avg, np.mean(vals))
tfs['U'] = avg
avg = []
for m in range(len(dict_nu['T'][num]['to'])):
        vals = [(dict_nu['T'][i]['to']/dict_nu['T'][i]['ti'])[m] for i in list(dict_nu['T'].keys())]
        avg = np.append(avg, np.mean(vals))
tfs['T'] = avg

with open('new_transfer_function_correction.pkl', 'wb') as handle:
        pickle.dump(tfs, handle, protocol=pickle.HIGHEST_PROTOCOL)

