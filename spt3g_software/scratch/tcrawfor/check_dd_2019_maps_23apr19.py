#import numpy as np
#import scipy
#from scipy import ndimage
#import pickle
#from spt3g import core, std_processing, gcp
#import os
#import glob
##from spt3g.mapspectra import basicmaputils as bmu
#from spt3g.mapmaker import mapmakerutils as mm
#from spt3g.mapmaker import summingmaps as sm
#from spt3g.simulations import quick_flatsky_routines as qfr
#from copy import copy
#
## get files
#fields = ['ra0hdec-44.75','ra0hdec-52.25','ra0hdec-59.75','ra0hdec-67.25']
#file_dict = {}
#for field in fields:
#    files = glob.glob('/spt/user/ddutcher/'+field+'/test_2019_maps/*/test_2019_maps*.g3')
#    files.sort()
#    file_dict[field] = files
#
## dictionary for individual-map stats
#stats_dict_onemap = {}
#for field in fields:
#    stats_dict_onemap[field] = {}
#
## dictionary for coadd stats
#stats_dict_coadd = {}
##nweeks = 1. # make a new coadd every this many weeks
#nweeks = 2./7. # make a new coadd every this many weeks
#
## hard-coded pixel ranges for individual-map stats
#xmin = {}
#xmin[fields[0]] = 1050
#xmin[fields[1]] = 1050
#xmin[fields[2]] = 1050
#xmin[fields[3]] = 1050
#xmax = {}
#xmax[fields[0]] = 1200
#xmax[fields[1]] = 1200
#xmax[fields[2]] = 1200
#xmax[fields[3]] = 1200
#ymin = {}
#ymin[fields[0]] = 1050
#ymin[fields[1]] = 825
#ymin[fields[2]] = 600
#ymin[fields[3]] = 375
#ymax = {}
#ymax[fields[0]] = 1200
#ymax[fields[1]] = 975
#ymax[fields[2]] = 750
#ymax[fields[3]] = 525
#
#ymin_coadd = 425
#ymax_coadd = 1175
#xmin_coadd = 750
#xmax_coadd = 1500
#
## loop over individual maps, get stats on each, decide if it can go in
## the coadd, do so if so.
#allfiles = []
#for field in fields:
#    for file1 in file_dict[field]:
#        if '66469280' not in file1 and '67626036' not in file1:
#            allfiles.append(file1)
#allfiles = np.asarray(allfiles)
#allobsids = np.asarray([file1.split('/')[-2] for file1 in allfiles])
#sobs = np.argsort(allobsids)
#allobsids = allobsids[sobs]
#allfiles = allfiles[sobs]
#obsid_min = np.float(allobsids[0])
#last_coadd = obsid_min
#rms_min_qu = 75.
#rms_max_qu = 750.
#maxwt = 1e4
#
### !!!
###allfiles = allfiles[227:240]
###allobsids = allobsids[227:240]
##allfiles = allfiles[227:]
##allobsids = allobsids[227:]
##allfiles = allfiles[112:]
##allobsids = allobsids[112:]
##allfiles = allfiles[338:350]
##allobsids = allobsids[338:350]
### !!!
#
#coadd_started = False
#mfkeys = ['T','Q','U','Wpol']
#for file1,obsid1 in zip(allfiles,allobsids):
#    print(obsid1)
#    if np.int(obsid1) < 70000000:
#        continue
#    for field in fields:
#        if file1 in file_dict[field]:
#            thisfield = field
#    print(thisfield)
#    stats_dict_onemap[thisfield][obsid1] = {}
#    f1 = core.G3File(file1)
#    frame150l = None
#    frame150r = None
#    for frame in f1:
#        if frame.type is core.G3FrameType.Map:
#            if frame['Id'] == 'Left150GHz':
#                frame150l = frame
#            if frame['Id'] == 'Right150GHz':
#                frame150r = frame
#    if frame150l is not None and frame150r is not None:
#        lframe2 = core.G3Frame(core.G3FrameType.Map)
#        for key in mfkeys:
#            lframe2[key] = frame150l[key]
#        rframe2 = core.G3Frame(core.G3FrameType.Map)
#        for key in mfkeys:
#            rframe2[key] = frame150r[key]
#        stats_dict_onemap[thisfield][obsid1]['weight_slice'] = np.asarray(frame150l['Wpol'].TT)[:,1125] + np.asarray(frame150r['Wpol'].TT)[:,1125]
#        mm.RemoveWeightModule(frame150l)
#        mm.RemoveWeightModule(frame150r)
#        ttemp = (np.asarray(frame150l['T'])[ymin[thisfield]:ymax[thisfield],xmin[thisfield]:xmax[thisfield]]-np.asarray(frame150r['T'])[ymin[thisfield]:ymax[thisfield],xmin[thisfield]:xmax[thisfield]])/2.
#        qtemp = (np.asarray(frame150l['Q'])[ymin[thisfield]:ymax[thisfield],xmin[thisfield]:xmax[thisfield]]-np.asarray(frame150r['Q'])[ymin[thisfield]:ymax[thisfield],xmin[thisfield]:xmax[thisfield]])/2.
#        utemp = (np.asarray(frame150l['U'])[ymin[thisfield]:ymax[thisfield],xmin[thisfield]:xmax[thisfield]]-np.asarray(frame150r['U'])[ymin[thisfield]:ymax[thisfield],xmin[thisfield]:xmax[thisfield]])/2.
#        # define "noise" as sqrt(cl) in 3000 < l < 5000 converted to uK-arcmin
#        cltemp1 = qfr.cl_flatsky(ttemp/core.G3Units.microkelvin,reso_arcmin=2.,delta_ell=100.,hann=True)
#        whint = np.where(np.logical_and(cltemp1['ell'] > 3000.,cltemp1['ell'] < 5000.))
#        stats_dict_onemap[thisfield][obsid1]['uk_arcmin_t'] = np.sqrt(np.mean(cltemp1['cl']['TT'][whint]))/(np.pi/180./60.)
#        cltemp2 = qfr.cl_flatsky(qtemp/core.G3Units.microkelvin,reso_arcmin=2.,delta_ell=100.,hann=True)
#        stats_dict_onemap[thisfield][obsid1]['uk_arcmin_q'] = np.sqrt(np.mean(cltemp2['cl']['TT'][whint]))/(np.pi/180./60.)
#        cltemp3 = qfr.cl_flatsky(utemp/core.G3Units.microkelvin,reso_arcmin=2.,delta_ell=100.,hann=True)
#        stats_dict_onemap[thisfield][obsid1]['uk_arcmin_u'] = np.sqrt(np.mean(cltemp3['cl']['TT'][whint]))/(np.pi/180./60.)
#        qutemp = (stats_dict_onemap[thisfield][obsid1]['uk_arcmin_q'] + stats_dict_onemap[thisfield][obsid1]['uk_arcmin_u'])/2.
#        print(qutemp)
#        print(np.max(stats_dict_onemap[thisfield][obsid1]['weight_slice']))
#        print(" ")
#        if qutemp > rms_min_qu and qutemp < rms_max_qu and np.max(stats_dict_onemap[thisfield][obsid1]['weight_slice']) < maxwt:
#            if coadd_started:
#                t_coadd_l += lframe2['T']
#                q_coadd_l += lframe2['Q']
#                u_coadd_l += lframe2['U']
#                wpol_coadd_l += lframe2['Wpol']
#                t_coadd_r += rframe2['T']
#                q_coadd_r += rframe2['Q']
#                u_coadd_r += rframe2['U']
#                wpol_coadd_r += rframe2['Wpol']
#            else:
#                t_coadd_l = lframe2['T']
#                q_coadd_l = lframe2['Q']
#                u_coadd_l = lframe2['U']
#                wpol_coadd_l = lframe2['Wpol']
#                t_coadd_r = rframe2['T']
#                q_coadd_r = rframe2['Q']
#                u_coadd_r = rframe2['U']
#                wpol_coadd_r = rframe2['Wpol']
#                coadd_started = True
#            if np.float(obsid1) - last_coadd > 86400.*7.*nweeks:
#                print("last coadd:"+str(last_coadd))
#                stats_dict_coadd[obsid1] = {}
#                coadd_l_fordiff = core.G3Frame(core.G3FrameType.Map)
#                coadd_l_fordiff['T'] = t_coadd_l
#                coadd_l_fordiff['Q'] = q_coadd_l
#                coadd_l_fordiff['U'] = u_coadd_l
#                coadd_l_fordiff['Wpol'] = wpol_coadd_l
#                mm.RemoveWeightModule(coadd_l_fordiff) 
#                coadd_r_fordiff = core.G3Frame(core.G3FrameType.Map)
#                coadd_r_fordiff['T'] = t_coadd_r
#                coadd_r_fordiff['Q'] = q_coadd_r
#                coadd_r_fordiff['U'] = u_coadd_r
#                coadd_r_fordiff['Wpol'] = wpol_coadd_r
#                mm.RemoveWeightModule(coadd_r_fordiff)
#                stats_dict_coadd[obsid1]['weight_slice'] = np.asarray(coadd_l_fordiff['Wpol'].TT)[:,1125] + np.asarray(coadd_r_fordiff['Wpol'].TT)[:,1125]
#                ttemp = (np.asarray(coadd_l_fordiff['T'])[ymin_coadd:ymax_coadd,xmin_coadd:xmax_coadd]-np.asarray(coadd_r_fordiff['T'])[ymin_coadd:ymax_coadd,xmin_coadd:xmax_coadd])/2.
#                qtemp = (np.asarray(coadd_l_fordiff['Q'])[ymin_coadd:ymax_coadd,xmin_coadd:xmax_coadd]-np.asarray(coadd_r_fordiff['Q'])[ymin_coadd:ymax_coadd,xmin_coadd:xmax_coadd])/2.
#                utemp = (np.asarray(coadd_l_fordiff['U'])[ymin_coadd:ymax_coadd,xmin_coadd:xmax_coadd]-np.asarray(coadd_r_fordiff['U'])[ymin_coadd:ymax_coadd,xmin_coadd:xmax_coadd])/2.
#                cltemp4 = qfr.cl_flatsky(ttemp/core.G3Units.microkelvin,reso_arcmin=2.,delta_ell=50.,hann=True)
#                whint = np.where(np.logical_and(cltemp4['ell'] > 3000.,cltemp4['ell'] < 5000.))
#                cltemp5 = qfr.cl_flatsky(qtemp/core.G3Units.microkelvin,reso_arcmin=2.,delta_ell=50.,hann=True)
#                cltemp6 = qfr.cl_flatsky(utemp/core.G3Units.microkelvin,reso_arcmin=2.,delta_ell=50.,hann=True)
#                stats_dict_coadd[obsid1]['t_cutout'] = ttemp
#                stats_dict_coadd[obsid1]['uk_arcmin_t'] = np.sqrt(np.mean(cltemp4['cl']['TT'][whint]))/(np.pi/180./60.)
#                stats_dict_coadd[obsid1]['q_cutout'] = qtemp
#                stats_dict_coadd[obsid1]['uk_arcmin_q'] = np.sqrt(np.mean(cltemp5['cl']['TT'][whint]))/(np.pi/180./60.)
#                stats_dict_coadd[obsid1]['u_cutout'] = utemp
#                stats_dict_coadd[obsid1]['uk_arcmin_u'] = np.sqrt(np.mean(cltemp6['cl']['TT'][whint]))/(np.pi/180./60.)
#                last_coadd = np.float(obsid1)
#                print("last coadd:"+str(last_coadd))
#
#obsids_field0 = []
#uk_arcmin_t_field0 = []
#uk_arcmin_q_field0 = []
#uk_arcmin_u_field0 = []
#weight_slice_field0 = []
#for obsid01 in stats_dict_onemap[fields[0]].keys():
#    if len(stats_dict_onemap[fields[0]][obsid01]) > 0:
#        obsids_field0.append(np.float(obsid01))
#        uk_arcmin_t_field0.append(stats_dict_onemap[fields[0]][obsid01]['uk_arcmin_t'])
#        uk_arcmin_q_field0.append(stats_dict_onemap[fields[0]][obsid01]['uk_arcmin_q'])
#        uk_arcmin_u_field0.append(stats_dict_onemap[fields[0]][obsid01]['uk_arcmin_u'])
#        weight_slice_field0.append(stats_dict_onemap[fields[0]][obsid01]['weight_slice'])
#obsids_field1 = []
#uk_arcmin_t_field1 = []
#uk_arcmin_q_field1 = []
#uk_arcmin_u_field1 = []
#weight_slice_field1 = []
#for obsid01 in stats_dict_onemap[fields[1]].keys():
#    if len(stats_dict_onemap[fields[1]][obsid01]) > 0:
#        obsids_field1.append(np.float(obsid01))
#        uk_arcmin_t_field1.append(stats_dict_onemap[fields[1]][obsid01]['uk_arcmin_t'])
#        uk_arcmin_q_field1.append(stats_dict_onemap[fields[1]][obsid01]['uk_arcmin_q'])
#        uk_arcmin_u_field1.append(stats_dict_onemap[fields[1]][obsid01]['uk_arcmin_u'])
#        weight_slice_field1.append(stats_dict_onemap[fields[1]][obsid01]['weight_slice'])
#obsids_field2 = []
#uk_arcmin_t_field2 = []
#uk_arcmin_q_field2 = []
#uk_arcmin_u_field2 = []
#weight_slice_field2 = []
#for obsid01 in stats_dict_onemap[fields[2]].keys():
#    if len(stats_dict_onemap[fields[2]][obsid01]) > 0:
#        obsids_field2.append(np.float(obsid01))
#        uk_arcmin_t_field2.append(stats_dict_onemap[fields[2]][obsid01]['uk_arcmin_t'])
#        uk_arcmin_q_field2.append(stats_dict_onemap[fields[2]][obsid01]['uk_arcmin_q'])
#        uk_arcmin_u_field2.append(stats_dict_onemap[fields[2]][obsid01]['uk_arcmin_u'])
#        weight_slice_field2.append(stats_dict_onemap[fields[2]][obsid01]['weight_slice'])
#obsids_field3 = []
#uk_arcmin_t_field3 = []
#uk_arcmin_q_field3 = []
#uk_arcmin_u_field3 = []
#weight_slice_field3 = []
#for obsid01 in stats_dict_onemap[fields[3]].keys():
#    if len(stats_dict_onemap[fields[3]][obsid01]) > 0:
#        obsids_field3.append(np.float(obsid01))
#        uk_arcmin_t_field3.append(stats_dict_onemap[fields[3]][obsid01]['uk_arcmin_t'])
#        uk_arcmin_q_field3.append(stats_dict_onemap[fields[3]][obsid01]['uk_arcmin_q'])
#        uk_arcmin_u_field3.append(stats_dict_onemap[fields[3]][obsid01]['uk_arcmin_u'])
#        weight_slice_field3.append(stats_dict_onemap[fields[3]][obsid01]['weight_slice'])
#
#pickle.dump(stats_dict_onemap,open('/spt/user/tcrawfor/public/check_dd_2019_maps_stats_onemap_30apr19.pkl','w'))
#pickle.dump(stats_dict_coadd,open('/spt/user/tcrawfor/public/check_dd_2019_maps_stats_coadd_30apr19.pkl','w'))



coadd_l_fordiff = core.G3Frame(core.G3FrameType.Map)
coadd_l_fordiff['T'] = t_coadd_l
coadd_l_fordiff['Q'] = q_coadd_l
coadd_l_fordiff['U'] = u_coadd_l
coadd_l_fordiff['Wpol'] = wpol_coadd_l
mm.RemoveWeightModule(coadd_l_fordiff) 
coadd_r_fordiff = core.G3Frame(core.G3FrameType.Map)
coadd_r_fordiff['T'] = t_coadd_r
coadd_r_fordiff['Q'] = q_coadd_r
coadd_r_fordiff['U'] = u_coadd_r
coadd_r_fordiff['Wpol'] = wpol_coadd_r
mm.RemoveWeightModule(coadd_r_fordiff)

stats_dict_final = {}
stats_dict_final['weight_slice'] = np.asarray(coadd_l_fordiff['Wpol'].TT)[:,1125] + np.asarray(coadd_r_fordiff['Wpol'].TT)[:,1125]
for thisfield in fields:
    stats_dict_final[thisfield] = {}
    ttemp = (np.asarray(coadd_l_fordiff['T'])[ymin[thisfield]:ymax[thisfield],xmin[thisfield]:xmax[thisfield]]-np.asarray(coadd_r_fordiff['T'])[ymin[thisfield]:ymax[thisfield],xmin[thisfield]:xmax[thisfield]])/2.
    qtemp = (np.asarray(coadd_l_fordiff['Q'])[ymin[thisfield]:ymax[thisfield],xmin[thisfield]:xmax[thisfield]]-np.asarray(coadd_r_fordiff['Q'])[ymin[thisfield]:ymax[thisfield],xmin[thisfield]:xmax[thisfield]])/2.
    utemp = (np.asarray(coadd_l_fordiff['U'])[ymin[thisfield]:ymax[thisfield],xmin[thisfield]:xmax[thisfield]]-np.asarray(coadd_r_fordiff['U'])[ymin[thisfield]:ymax[thisfield],xmin[thisfield]:xmax[thisfield]])/2.
    # define "noise" as sqrt(cl) in 3000 < l < 5000 converted to uK-arcmin
    cltemp1 = qfr.cl_flatsky(ttemp/core.G3Units.microkelvin,reso_arcmin=2.,delta_ell=100.,hann=True)
    whint = np.where(np.logical_and(cltemp1['ell'] > 3000.,cltemp1['ell'] < 5000.))
    stats_dict_final[thisfield]['uk_arcmin_t'] = np.sqrt(np.mean(cltemp1['cl']['TT'][whint]))/(np.pi/180./60.)
    cltemp2 = qfr.cl_flatsky(qtemp/core.G3Units.microkelvin,reso_arcmin=2.,delta_ell=100.,hann=True)
    stats_dict_final[thisfield]['uk_arcmin_q'] = np.sqrt(np.mean(cltemp2['cl']['TT'][whint]))/(np.pi/180./60.)
    cltemp3 = qfr.cl_flatsky(utemp/core.G3Units.microkelvin,reso_arcmin=2.,delta_ell=100.,hann=True)
    stats_dict_final[thisfield]['uk_arcmin_u'] = np.sqrt(np.mean(cltemp3['cl']['TT'][whint]))/(np.pi/180./60.)

