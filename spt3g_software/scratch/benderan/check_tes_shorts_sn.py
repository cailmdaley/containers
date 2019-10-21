from spt3g import core, dfmux,  std_processing
import os
import numpy as np



caldir='/poleanalysis/sptdaq/calresult/calibration/calibrator'
cal_files=os.listdir(caldir)
date_start=std_processing.time_to_obsid('20181231_040052')
date_stop=std_processing.time_to_obsid('20190102_040052')


gfl=[]
for fl in cal_files:
    if (np.int32(fl.split('.')[0]) > date_start) and (np.int32(fl.split('.')[0]) < date_stop):
        gfl.append(fl)
gfl.sort()
gfl=np.array(gfl)
cal_freq=np.zeros(len(gfl))


for kk, gff in enumerate(gfl):
    frames=list(core.G3File(os.path.join(caldir,gff)))
    cal_frame=frames[0]
    if kk==0:
        bolos=cal_frame['CalibratorResponseSN'].keys()    
        cal_sn=np.zeros((len(bolos),len(gfl)))
        n_alive=np.zeros(len(gfl))
    cal_freq[kk]=cal_frame['CalibratorResponseFrequency']/core.G3Units.Hz
    print gff, cal_frame['CalibratorResponseFrequency']/core.G3Units.Hz
    for ii,bolo in enumerate(bolos):
        cal_sn[ii,kk]=cal_frame['CalibratorResponseSN'][bolo]
    n_alive[kk]=len(np.where(cal_sn[:,kk]>20)[0])

gobs=np.where(cal_freq < 10)[0]


bolos=np.array(bolos)


blist_short_master=np.array([])
autoprobe_exclude_dir='/home/benderan/hardware_maps_southpole/2019/global/exclude/autoprobe_exclude/short'
probe_files=os.listdir(autoprobe_exclude_dir)
for fl in probe_files:
    if 'short' in fl:
        blist_short=create_blist_from_csv(filein=os.path.join(autoprobe_exclude_dir,fl))
        blist_short_master=np.append(blist_short_master,np.array(blist_short))


blist_inds=[]
for kk in blist_short_master:
    if len(np.where(bolos==kk)[0])>0:
        blist_inds.append(np.where(bolos==kk)[0][0])
blist_inds=np.array(blist_inds)

for obs in gobs:
    tmp_slice=cal_sn[blist_inds,obs]
    print '%i NaN,  %i with S/N < 20, %i with  S/N > 20' % (len(np.where(np.isnan(tmp_slice))[0]), len(np.where(tmp_slice <20)[0]), len(np.where(tmp_slice >=20)[0]))



