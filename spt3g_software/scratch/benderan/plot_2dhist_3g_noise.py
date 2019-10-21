import numpy as np
import matplotlib.pyplot as plt
from spt3g import core, calibration
import os



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


def plot_2dhist_noise():
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
    binedges_f=np.arange(1.4,5.5,0.2)
    binedges_n=np.arange(0,50,1) 
    plt.hist2d(biasfreq[ginds],wnl[ginds],[binedges_f,binedges_n],cmap='bone')
    plt.ylim(0,30)
    plt.xlabel('Bias Frequency (MHz)')
    plt.ylabel('White Noise Level (pA/$\sqrt{Hz}$)')
    plt.colorbar(label = 'Number of Detetors')
    plt.clim(0,300)
    return
