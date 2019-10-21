import numpy as np
import matplotlib.pyplot as plt
from spt3g import core, calibration
import os

# some useful noise functions
def readout_noise(x, readout):
    return np.sqrt(readout)*np.ones(len(x))
def photon_noise(x, photon, tau):
    return np.sqrt(photon / (1 + 2*np.pi*((x*tau)**2)))
def atm_noise(x, A, alpha):
    return np.sqrt(A * (x)**(-1*alpha))
def noise_model(x, readout, A, alpha, photon, tau):
    return np.sqrt(readout + (A * (x)**(-1*alpha)) + photon / (1 + 2*np.pi*((x*tau)**2)))
def horizon_model(x, readout, A, alpha):
    return np.sqrt(readout + (A * (x)**(-1*alpha)))
def knee_func(x, readout, A, alpha, photon, tau):
    return (A * (x)**(-1*alpha)) - photon / (1 + 2*np.pi*((x*tau)**2)) - readout
def horizon_knee_func(x, readout, A, alpha):
    return (A * (x)**(-1*alpha)) - readout

horizondatafile='/home/adama/SPT/spt_analysis/20190329_gainmatching/horizon_noise_77863968_bender_ltd.g3'


fr = list(core.G3File(horizondatafile))[1]

band_numbers = {90.: 1, 150.: 2, 220.: 3}
subplot_numbers = {90.: 1, 150.: 1, 220.: 1}
band_labels = {90:'95', 150:'150', 220:'220'}

fig, ax = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(12,4))
fig.subplots_adjust(wspace=0)
for jband, band in enumerate([90., 150., 220.]):
    
    group = '{:.1f}_w180'.format(band)

    ff_diff = np.array(fr['AverageASDDiff']['frequency']/core.G3Units.Hz)
    ff_sum = np.array(fr['AverageASDSum']['frequency']/core.G3Units.Hz)
    asd_diff = np.array(fr['AverageASDDiff'][group]) / np.sqrt(2.)
    asd_sum = np.array(fr['AverageASDSum'][group]) / np.sqrt(2.)

    par = fr["AverageASDDiffFitParams"][group]
    ax[jband].loglog(ff_sum[ff_sum<75], asd_sum[ff_sum<75],
                     label='pair sum (measured)', color='0.6')
    ax[jband].loglog(ff_diff[ff_diff<75], asd_diff[ff_diff<75],
                     label='pair difference (measured)', color='k')
    ax[jband].loglog(ff_sum, atm_noise(ff_sum, par[1], par[2]) / np.sqrt(2.),
                     'C0:', label='low-frequency noise',linewidth=2)
    ax[jband].loglog(ff_sum, readout_noise(ff_sum, par[0]) / np.sqrt(2.),
                     'C2:', label='white noise',linewidth=2)
    ax[jband].loglog(ff_sum, horizon_model(ff_sum, *list(par)) / np.sqrt(2.),
                     'C1--', label='total noise model',linewidth=2)

    ax[jband].set_title('{} GHz'.format(band_labels[band]),fontsize=14)
    ax[jband].set_xlabel('Frequency (Hz)',fontsize=14)
    ax[jband].tick_params(labelsize=14)
    ax[jband].grid()
ax[0].set_ylabel('NEI (pA/$\sqrt{Hz}$)',fontsize=14)

plt.ylim([5,1000])
plt.legend()
#xtks=plt.xticks()
#plt.xticks(xtks,fontsize=14)
#ytks=plt.yticks
#plt.yticks(ytks,fontsize=14)
plt.tight_layout()
plt.savefig('/home/benderan/spt_analysis/w180_horizon_noise_ltd_v2.pdf')

horizondatafile='/home/benderan/spt_analysis/horizon_noise_77863968_bender_ltd_perbolo_only.g3'
fr = list(core.G3File(horizondatafile))[1]
bolokeys=fr['ASDFitParams'].keys()
param1=np.zeros(len(bolokeys))
param2=np.zeros(len(bolokeys))
param3=np.zeros(len(bolokeys))
onefknee=np.zeros(len(bolokeys))  # defined as where wnl = 1/f noise level in PSD (corner frequency)

for ii,kk in enumerate(bolokeys):
    if len(fr['ASDFitParams'][kk])>0:
        param1[ii]=fr['ASDFitParams'][kk][0]
        param2[ii]=fr['ASDFitParams'][kk][1]
        param3[ii]=fr['ASDFitParams'][kk][2]
        onefknee[ii]= (param2[ii]/param1[ii])**(1./param3[ii])
        # the model is ASD = np.sqrt(param1 + param2*frequency^(-1*param3))

wnl=np.sqrt(param1)

obsband=np.zeros(len(bolokeys))
wafer=np.zeros(len(bolokeys),dtype='S4')
calsn=np.zeros(len(bolokeys))
biasfreq=np.zeros(len(bolokeys))
biasamp=np.zeros(len(bolokeys))
bias_state=np.zeros(len(bolokeys),dtype='S10')
# load in the offline calibration file used in the analysis
calfile='/spt/data/bolodata/downsampled/noise/73798315/offline_calibration.g3'
caldata=list(core.G3File(calfile))[0]

fulldata=list(core.G3File(os.path.join(calfile.strip('offline_calibration.g3'),'0000.g3')))
for ff in fulldata:
    if ff.type==core.G3FrameType.Wiring:
        wmap=ff['WiringMap']
    if ff.type==core.G3FrameType.Scan:
        hkdata=ff['DfMuxHousekeeping']

master_chlist,master_freq,master_amp,master_state=make_hk_map(hkdata)
namb=0
gpsd=np.zeros(len(bolokeys))+1
for ii,kk in enumerate(bolokeys):
    obsband[ii]=caldata['BolometerProperties'][kk].band
    wafer[ii]=caldata['BolometerProperties'][kk].wafer_id
    calsn[ii]=caldata['CalibratorResponseSN'][kk]
    if np.all(np.abs(fr['ASD'][kk])<1e-10):
        gpsd[ii] =0
    bfreq,camp,state=get_dfmux_information(kk,wmap,master_chlist,master_freq,master_amp,master_state)
    biasfreq[ii]=bfreq
    biasamp[ii]=camp
    bias_state[ii]=state
    if bfreq==-1:
        namb+=1

calsn[np.isnan(calsn)] = 0.

# let's isolate the untuned etc
ginds=np.intersect1d(np.intersect1d(np.intersect1d(np.where(biasfreq>0),np.where(biasamp >0)),np.where(bias_state == b'tuned')),np.where(calsn >20.))

ginds=np.intersect1d(np.intersect1d(np.where(calsn>20.)[0],np.where(gpsd ==1)[0]),np.where(wnl > 1e-2)[0])

inds90=np.where(obsband==900.)[0]
inds150=np.where(obsband==1500.)[0]
inds220=np.where(obsband==2200.)[0]




plt.figure()
plt.plot(biasfreq[ginds],wnl[ginds],'k.')
plt.xlabel('Bias Frequency (MHz)')
plt.ylabel('White Noise Level (pA/$\sqrt{Hz}$)')
plt.ylim(0,100)


plt.savefig('/home/benderan/spt_analysis/nei_vs_bfreq_horizon_77863968.png')

plt.figure()
plt.plot(biasfreq[ginds],onefknee[ginds],'k.')
plt.xlabel('Bias Frequency (MHz)')
plt.ylabel('1/f knee (Hz)')
plt.ylim(0,1)
plt.savefig('/home/benderan/spt_analysis/onefknee_vs_bfreq_horizon_77863968.png')

plt.figure()
plt.plot(biasfreq[ginds],np.sqrt(param2[ginds]),'k.')
plt.xlabel('Bias Frequency (MHz)')
plt.ylabel('$\sqrt{1/f Amplitude}$')
plt.ylim(0,40)
plt.savefig('/home/benderan/spt_analysis/onefamp_vs_bfreq_horizon_77863968.png')

wnlbins=np.arange(0,35,1)
plt.figure()
plt.hist(wnl[ginds],wnlbins,label='Measured')#,histtype='step',linewidth=2.)
plt.grid('on')
plt.xticks(np.arange(0,40,5),fontsize=14)
plt.yticks(np.arange(0,2000,500),fontsize=14)
plt.ylim(0,1700)
ylm=plt.ylim()
plt.plot(np.array([1,1])*18.9, np.array([0,ylm[1]]),'--',label='150 GHz\n Photon Noise')
plt.legend(fontsize=12)
#plt.legend((p1,p2),('Measured', '150 GHz Photon Noise'),loc='best')
plt.xlabel('White Noise Level (pA/$\sqrt{Hz}$)',fontsize=14)
plt.ylabel('Number of Readout Channels',fontsize=14)
plt.xlim(0,35)
plt.tight_layout()
plt.savefig('/home/benderan/spt_analysis/nei_hist_ltd_horizon_77863968.png')

onefkneebins=np.arange(0,0.4,0.01)
plt.figure()
plt.hist(onefknee[ginds],onefkneebins)
plt.grid('on',which='both')
plt.xlim(0,0.4)
plt.xticks(np.arange(0,0.5,0.1),fontsize=14)
from matplotlib.ticker import MultipleLocator
plt.axes().xaxis.set_minor_locator(MultipleLocator(0.05))
plt.yticks(np.arange(0,3000,500),fontsize=14)
plt.xlabel('1/f Knee Frequency (Hz)',fontsize=14)
plt.ylabel('Number of Readout Channels',fontsize=14)
plt.tight_layout()
plt.savefig('/home/benderan/spt_analysis/onefknee_hist_ltd_horizon_77863968.png')

wnlbins=np.arange(0,35,1)
plt.figure()
plt.hist(wnl[np.intersect1d(ginds,inds90)],wnlbins,label='95 GHz',alpha=0.5)
plt.hist(wnl[np.intersect1d(ginds,inds150)],wnlbins,label='150 GHz',alpha=0.5)
plt.hist(wnl[np.intersect1d(ginds,inds220)],wnlbins,label='220 GHz',alpha=0.5)
#plt.hist(wnl[np.intersect1d(ginds,inds90)],wnlbins,histtype='step',linewidth=2,label='95 GHz')
#plt.hist(wnl[np.intersect1d(ginds,inds150)],wnlbins,label='150 GHz',histtype='step',linewidth=2)
#plt.hist(wnl[np.intersect1d(ginds,inds220)],wnlbins,label='220 GHz',histtype='step',linewidth=2)

plt.xticks(np.arange(0,40,5),fontsize=14)
plt.yticks(np.arange(0,2000,500),fontsize=14)
plt.ylim(0,1700)
ylm=plt.ylim()
plt.plot(np.array([1,1])*15.7,np.array([0,ylm[1]]),'--',color='C0')
plt.plot(np.array([1,1])*18.8,np.array([0,ylm[1]]),'--',color='C1')
plt.plot(np.array([1,1])*18.3, np.array([0,ylm[1]]),'--',color='C2')
#plt.plot(np.array([1,1])*15.7,np.array([0,ylm[1]]),'--',label='95 GHz\n Photon Noise',color='C0')
#plt.plot(np.array([1,1])*18.8,np.array([0,ylm[1]]),'--',label='150 GHz\n Photon Noise',color='C1')
#plt.plot(np.array([1,1])*18.3, np.array([0,ylm[1]]),'--',label='220 GHz\n Photon Noise',color='C2')

plt.legend(fontsize=12)
plt.xlabel('White Noise Level (pA/$\sqrt{Hz}$)',fontsize=14)
plt.ylabel('Number of Readout Channels',fontsize=14)
plt.xlim(0,35)
plt.tight_layout()
plt.savefig('/home/benderan/spt_analysis/nei_hist_3band_ltd_horizon_77863968.png')

print('median white noise levels')
print(np.median(wnl[np.intersect1d(ginds,inds90)]))
print(np.median(wnl[np.intersect1d(ginds,inds150)]))
print(np.median(wnl[np.intersect1d(ginds,inds220)]))
print('median 1/f knee')
print(np.median(onefknee[np.intersect1d(ginds,inds90)]))
print(np.median(onefknee[np.intersect1d(ginds,inds150)]))
print(np.median(onefknee[np.intersect1d(ginds,inds220)]))



onefkneebins=np.arange(0,0.4,0.01)
plt.figure()
plt.hist(onefknee[np.intersect1d(ginds,inds90)],onefkneebins,label='95 GHz',alpha=0.5)
plt.hist(onefknee[np.intersect1d(ginds,inds150)],onefkneebins,label='150 GHz',alpha=0.5)
plt.hist(onefknee[np.intersect1d(ginds,inds220)],onefkneebins,label='220 GHz',alpha=0.5)
plt.grid('on',which='both')
plt.xlim(0,0.4)
plt.xticks(np.arange(0,0.5,0.1),fontsize=14)
from matplotlib.ticker import MultipleLocator
plt.axes().xaxis.set_minor_locator(MultipleLocator(0.05))
plt.yticks(np.arange(0,2000,500),fontsize=14)
plt.xlabel('1/f Knee Frequency (Hz)',fontsize=14)
plt.ylabel('Number of Readout Channels',fontsize=14)
plt.tight_layout()
plt.legend()
plt.savefig('/home/benderan/spt_analysis/onefknee_hist_3band_ltd_horizon_77863968.png')


param2bins=np.arange(0,np.max(param2),1000)
plt.figure()
plt.hist(param2[np.intersect1d(ginds,inds90)],param2bins,label='95 GHz',alpha=0.5)
plt.hist(param2[np.intersect1d(ginds,inds150)],param2bins,label='150 GHz',alpha=0.5)
plt.hist(param2[np.intersect1d(ginds,inds220)],param2bins,label='220 GHz',alpha=0.5)
plt.tight_layout()
plt.legend()




def get_dfmux_information(bolokey,wmap,master_chlist,master_freq,master_amp,master_state):
    #crate=wmap[bolokey].crate_serial
    #slot=wmap[bolokey].board_slot
    squid=wmap[bolokey].module
    channel=wmap[bolokey].channel
    serial=wmap[bolokey].board_serial
    search_str=str(serial)+'/'+str(squid+1)+'/'+str(channel+1)
    search_result=np.where(master_chlist==search_str)[0]
    if len(search_result)==1:
        bfreq=master_freq[search_result[0]]
        camp=master_amp[search_result[0]]
        state=master_state[search_result[0]]
    else:
        print('ambiguous search: '+bolokey+':'+search_str)
        bfreq=-1
        camp=-1
        state='-1'
    return bfreq,camp,state


def make_hk_map(hkdata):
    mezz=[1,2]
    module=[1,2,3,4]

    boards=hkdata.keys()
    master_chlist=[]
    master_freq=[]
    master_amp=[]
    master_state=[]

    for kk in boards:
        for mz in mezz:
            for md in module:
                tmp= hkdata[kk].mezz[mz].modules[md]
                chlist=tmp.channels.keys()
                for ch in chlist:
                    sq=(mz-1)*4+md
                    master_chlist.append(hkdata[kk].serial.lstrip('0')+'/'+str(sq)+'/'+str(ch))
                    master_amp.append(tmp.channels[ch].carrier_amplitude)  ### note that this is in normalized units
                    master_freq.append(tmp.channels[ch].carrier_frequency/core.G3Units.MHz)   # in Mhz
                    master_state.append(tmp.channels[ch].state)

    master_chlist=np.array(master_chlist)
    master_freq=np.array(master_freq)
    master_amp=np.array(master_amp)
    master_state=np.array(master_state)

    return master_chlist,master_freq,master_amp,master_state
    

