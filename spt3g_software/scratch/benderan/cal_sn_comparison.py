from spt3g import core, dfmux, std_processing
import numpy as np
import os

calibrator_dir='/spt/user/production/calibration/calibrator/'
cal_files=os.listdir(calibrator_dir)

ngfl=[]
gfl=[]
for fl in cal_files:
    #    if (np.int32(fl.split('.')[0]) >8007450) and (np.int32(fl.split('.')[0]) < 15007453):
    #        gfl.append(fl)
    if (np.int32(fl.split('.')[0]) >75024859) and (np.int32(fl.split('.')[0])<81024859):
        ngfl.append(fl)
gfl.sort()
ngfl.sort()
gfl=np.array(gfl)
ngfl=np.array(ngfl)



ngfl=ngfl[np.arange(0,len(ngfl),15)]

#frames=list(core.G3File(os.path.join(calibrator_dir,'15796205.g3')))

for kk, gff in enumerate(ngfl):
    frames=list(core.G3File(os.path.join(calibrator_dir,gff)))
    cal_frame=frames[0]
    if kk==0:
        bolos=cal_frame['CalibratorResponseSN'].keys()    
        cal_sn_new_noise=np.zeros((len(bolos),len(ngfl)))
        cal_s_new_noise=np.zeros((len(bolos),len(ngfl)))
        n_alive=np.zeros(len(ngfl))
        cal_freq=np.zeros(len(ngfl))
    #cal_freq=4.
    #if np.abs(cal_frame['CalibratorResponseFrequency']/core.G3Units.Hz-cal_freq)/cal_freq > 0.01:
    #    continue
    if 'CalibratorResponseFrequency' not in cal_frame.keys():
        continue
    print(gff, cal_frame['CalibratorResponseFrequency']/core.G3Units.Hz)
    if cal_frame['CalibratorResponseFrequency']/core.G3Units.Hz > 4.1:
        continue
    for ii,bolo in enumerate(bolos):
        if bolo in cal_frame['CalibratorResponseSN'].keys():
            cal_sn_new_noise[ii,kk]=cal_frame['CalibratorResponseSN'][bolo]
            cal_s_new_noise[ii,kk]=cal_frame['CalibratorResponse'][bolo]
    n_alive[kk]=len(np.where(cal_sn_new_noise[:,kk]>20)[0])
    cal_freq[kk]=cal_frame['CalibratorResponseFrequency']/core.G3Units.Hz

cal_sn_old_noise=np.zeros((len(bolos),len(gfl)))-1
cal_s_old_noise=np.zeros((len(bolos),len(gfl)))-1
bias_freq=np.zeros(len(bolos))
for kk,gff in enumerate(gfl):
    print(kk)
    frames=list(core.G3File(os.path.join(calibrator_dir,gff)))
    cal_frame=frames[0]
    if ('CalibratorResponseSN' not in cal_frame) or ('CalibratorResponseFrequency' not in cal_frame):
        continue
    if np.abs(cal_frame['CalibratorResponseFrequency']/core.G3Units.Hz-6.)/6. > 0.01:
        continue
    else:
        print(cal_frame['CalibratorResponseFrequency']/core.G3Units.Hz)
    for ii,bolo in enumerate(bolos):
        if bolo in cal_frame['CalibratorResponseSN'].keys():
            cal_sn_old_noise[ii,kk]=cal_frame['CalibratorResponseSN'][bolo]
            cal_s_old_noise[ii,kk]=cal_frame['CalibratorResponse'][bolo]
        bias_freq[ii]=np.float(bolo.split('.')[-1])/1e3


cal_sn_max=np.zeros(len(bolos))
cal_s_max=np.zeros(len(bolos))
cal_sn_med=np.zeros(len(bolos))
cal_sn_min=np.zeros(len(bolos))
max_id=np.zeros(len(bolos),dtype=np.int32)
for ii in range(len(bolos)):
    ginds=np.intersect1d(np.intersect1d(np.where(np.isnan(cal_sn_old_noise[ii,:])==False)[0],np.where(cal_sn_old_noise[ii,:] !=-1)[0]),np.where(np.isinf(cal_sn_old_noise[ii,:])==False)[0])
    if len(ginds)>0:
        cal_sn_max[ii]=np.max(cal_sn_old_noise[ii,ginds])
        max_id[ii]=np.min(np.where(cal_sn_old_noise[ii,:]==cal_sn_max[ii])[0])
        cal_s_max[ii]=cal_s_old_noise[ii,max_id[ii]]
        cal_sn_med[ii]=np.median(cal_sn_old_noise[ii,ginds])
        cal_sn_min[ii]=np.min(cal_sn_old_noise[ii,ginds])



#frames=list(core.G3File(os.path.join('/spt/data/bolodata/downsampled/calibrator', ngfl[0].split('.')[0],'0000.g3')))


elnod_obsid='16523872'
elnod_dir='/spt/user/production/calibration/elnod'
elnod=list(core.G3File(os.path.join(elnod_dir,elnod_obsid+'.g3')))

elnod_datadir='/spt/data/bolodata/downsampled/elnod'
elnod_data=list(core.G3File(os.path.join(elnod_datadir,elnod_obsid,'0000.g3')))

elnod_sn=np.zeros(len(bolos))
for ii,bb in enumerate(bolos):
    elnod_sn[ii]=elnod[0]['ElnodSNSlopes'][bb]


caldata=list(core.G3File('/spt/user/production/calibration/calframe/RCW38-pixelraster/36985606.g3'))
cal_frame=caldata[0]

kys=cal_frame['CalibratorResponseSN'].keys()
cal_sn=np.zeros(len(kys))

wafer=np.zeros(len(cal_frame['CalibratorResponseSN'].keys()))
band=np.zeros(len(cal_frame['CalibratorResponseSN'].keys()))
pixel=np.zeros(len(cal_frame['CalibratorResponseSN'].keys()))
pol=np.zeros(len(cal_frame['CalibratorResponseSN'].keys()))
xoffset=np.zeros(len(cal_frame['CalibratorResponseSN'].keys()))
yoffset=np.zeros(len(cal_frame['CalibratorResponseSN'].keys()))

for ii,ky in enumerate(kys):
    cal_sn[ii]=cal_frame['CalibratorResponseSN'][ky]
    if ky in cal_frame['BolometerProperties'].keys():
        tmp=cal_frame['BolometerProperties'][ky]
        phys_name=tmp.physical_name
        wtmp=phys_name.split('_')[0]
        wafer[ii]=np.int32(wtmp.split('w')[1])
        band[ii]=tmp.band
        pixel[ii]=np.int32(phys_name.split('_')[1].split('.')[0])
        if phys_name.split('.')[-1] =='y':
            pol[ii]=1
        xoffset[ii]=tmp.x_offset
        yoffset[ii]=tmp.y_offset

alive_inds=np.where(cal_sn > 20.)[0]
for jj in alive_inds:
    xcent=xoffset[jj]/core.G3Units.arcmin
    ycent=yoffset[jj]/core.G3Units.arcmin
    if band[jj]==900.:
        continue
        col='r'
        if pol[jj]==1:
            phi=90.
        else:
            phi=0
    elif band[jj]==1500:
        continue
        col='g'
        phi_offset=0. #30.
        if pol[jj]==1:
            phi=phi_offset+90.
        else:
            phi=phi_offset
    elif band[jj]==2200:

        col='b'
        phi_offset=0. #60.
        if pol[jj]==1:
            phi=phi_offset+90.
        else:
            phi=phi_offset
    plt.plot(np.array([xcent-np.cos(phi/180*np.pi)*0.5,xcent+np.cos(phi/180*np.pi)*0.5]),np.array([ycent- np.sin(phi/180*np.pi)*0.5,ycent+ np.sin(phi/180*np.pi)*0.5]),col)
    plt.xlim(-70,70)
    plt.ylim(-50,60)
    plt.xticks([])
    plt.yticks([])






wafers=np.unique(wafer)
x90=[]
y90=[]
x150=[]
y150=[]
x220=[]
y220=[]
for ww in wafers:
    ginds=np.where(wafer==ww)[0]
    pixels=np.unique(pixel[ginds])
    for pp in pixels:
        gpix=np.intersect1d(ginds,np.where(pixel==pp)[0])
        binds=np.where(band[gpix]==900.)[0]
        if len(binds)==2:
            x90.append(xoffset[gpix][0])
            y90.append(yoffset[gpix][0])
        binds=np.where(band[gpix]==1500.)[0]
        if len(binds)==2:
            x150.append(xoffset[gpix][0])
            y150.append(yoffset[gpix][0])
        binds=np.where(band[gpix]==2200.)[0]
        if len(binds)==2:
            x220.append(xoffset[gpix][0])
            y220.append(yoffset[gpix][0])


plt.figure()
for ww in wafers:
    ginds=np.where(wafer==ww)[0]
    xsets=xoffset[ginds]
    ysets=yoffset[ginds]
    xmean=np.mean(xsets[np.where(np.isnan(xsets)==False)[0]])
    ymean=np.mean(ysets[np.where(np.isnan(ysets)==False)[0]])
    plt.text(xmean,ymean, np.str(np.int(ww)))
    plt.grid('on')

plt.text(-0.012,0.016,'W181',color='r')
plt.text(-0.0,0.014,'W174',color='r')
plt.text(0.01,0.012,'W188',color='r')
plt.text(-0.018,0.01,'W203',color='r')
plt.text(-0.005,0.0065,'W172',color='r')
plt.text(0.002,0.003,'W177',color='r')
plt.text(0.015,0.002,'W176',color='r')
plt.text(-0.014,-0.011,'W201',color='r')
plt.text(-0.007,-0.012,'W180',color='r')
plt.text(0.012,-0.01,'W187',color='r')


