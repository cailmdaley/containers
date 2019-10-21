import numpy as np
import pickle as pk
import os
from glob import glob
import sptpol_software.analysis.pointing.pointing_tools as pt


def gather_offsets(offset_dir = '/data55/bleeml/SPTpol_500d_Nov16/pointing_checks/2015_try_Dec31',
                   snr_cut = 6., hit_cut=5, offset_cut=3.,
                   search_str='smap_set*.txt',
                   file_start=0, file_stop=None):
    #offset_dir = '/data55/bleeml/SPTpol_500d_Nov16/pointing_checks/2015_try2'
    #offset_dir = '/data55/bleeml/SPTpol_500d_Nov16/pointing_checks/2016'

    psfiles = np.sort(glob(os.path.join(offset_dir,search_str)))

    ara = [np.nan]
    adec = [np.nan]
    dra = [np.nan]
    ddec = [np.nan]
    err_ra = [np.nan]
    err_dec = [np.nan]

    if not file_stop:
        file_stop = len(psfiles)

    for i in range(file_start,file_stop):
        this_ara, this_adec, this_dra, this_ddec, this_err_ra, this_err_dec = pt.get_point_source_offsets(psfiles[i], 
                                                                                                          snr_cut=snr_cut, 
                                                                                                          offset_cut=offset_cut)
        ara += this_ara                                        
        adec += this_adec
        dra += this_dra                                   
        ddec += this_ddec                                       
        err_ra += this_err_ra
        err_dec += this_err_dec

        deltas = {}
        errors = {}
    for i in range(1,len(ara)):
        if not (ara[i],adec[i]) in deltas.keys():
            deltas[(ara[i],adec[i])] = [(dra[i],ddec[i])]
            errors[(ara[i],adec[i])] = [(err_ra[i],err_dec[i])]
        else:            
            deltas[(ara[i],adec[i])] += [(dra[i],ddec[i])]
            errors[(ara[i],adec[i])] += [(err_ra[i],err_dec[i])]


    for key in deltas.keys():
        deltas[key] = np.array(deltas[key])
        errors[key] = np.array(errors[key])
        if len(deltas[key]) < int(hit_cut):
            deltas.pop(key, None)
            errors.pop(key, None)

    #Make inverse-variance weights
    weights = {}
    for key in errors.keys():
        weights[key] = 1./errors[key]**2.

    #Calculate weighted mean and std.
    mean_deltas = {}
    std_deltas = {}
    for key in deltas.keys():
        mean_deltas[key] = np.average(deltas[key],axis=0, weights=weights[key])
        diffs = deltas[key] - mean_deltas[key]
        std_deltas[key] = np.sqrt(np.sum(weights[key]*diffs**2, axis=0)/np.sum(weights[key],axis=0))


    #Offsets and jitter in arcseconds
    offsets = []
    jitter = []
    for key in deltas.keys():
        offsets.append(mean_deltas[key]*3600.)
        jitter.append(std_deltas[key]*3600.)

    offsets = np.array(offsets)
    jitter = np.array(jitter)

    mean_offset = np.mean(offsets)
    mean_jitter = np.mean(jitter)

    ara = []
    adec = []
    ra_offsets = []
    dec_offsets = []
    err = []
    for key in deltas.keys():
        ara.append(float(key[0]))
        adec.append(float(key[1]))
        ra_offsets.append(mean_deltas[key][0])
        dec_offsets.append(mean_deltas[key][1])
        err.append(np.mean(std_deltas[key]))

    return deltas, -np.array(ara), -np.array(adec), -np.array(ra_offsets), np.array(dec_offsets), np.array(err), offsets, jitter
                    
