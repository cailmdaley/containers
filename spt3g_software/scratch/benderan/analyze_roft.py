from spt3g import core, dfmux, std_processing
from spt3g.util.extractdata import extract_keys
import numpy as np
import matplotlib.pyplot as plt
import os,pickle,pdb
import pandas as pd

def analyze_and_plot_roft(obsid_files='',outputfile='',plot_r_t=False, bololist=[]):
    #obsid_files='/poleanalysis/benderan/roft_62754248_c5_s16.g3'
    frames=list(core.G3File(obsid_files))

    for ff in frames:
        if ff.type==core.G3FrameType.Wiring:
            wmap=ff['WiringMap']

    bolonames=np.array(wmap.keys())
    crates=np.array([])
    slots=np.array([])
    squids=np.array([])
    channels=np.array([])
    bias_freq=np.array([])


    for bb in bolonames:
        wtmp=wmap[bb]
        crates=np.append(crates,wtmp.crate_serial)
        slots=np.append(slots,wtmp.board_slot)
        squids=np.append(squids,wtmp.module)
        channels=np.append(channels,wtmp.channel)

    ts=extract_keys(obsid_files,'CalTimestreams_I_subset')

    bkeys=ts.keys()
    if bololist=='all':
        bkeys=np.array(bkeys)
    elif len(bololist)>0:
        bkeys=np.intersect1d(bololist,np.array(bkeys))


    if len(bkeys)== 0:
        raise ValueError('incorrect data file selected for this list')


    rstart=np.zeros(len(bkeys))
    rstart_sig=np.zeros(len(bkeys))
    rstop=np.zeros(len(bkeys))
    rstop_sig=np.zeros(len(bkeys))

    nsamps_start=np.arange(10000)
    nsamps_stop=np.arange(len(ts[bkeys[0]])-10000,len(ts[bkeys[0]]))
    # PathStringForBolo

    if plot_r_t:
        plt.figure()
        for ind,bb_select in enumerate(bkeys):
            plt.plot(ts[bb_select])
    for ind,bb_select in enumerate(bkeys):
        rstart[ind]=np.median(ts[bb_select][nsamps_start])
        rstart_sig[ind]=np.std(ts[bb_select][nsamps_start])
        rstop[ind]=np.median(ts[bb_select][nsamps_stop])
        rstop_sig[ind]=np.std(ts[bb_select][nsamps_stop])
                
        plt.xlabel('Time')
        plt.ylabel('Resistance ($\Omega$)')
        plt.title('Crate '+str(wmap[bb_select].crate_serial)+', Slot '+str(wmap[bb_select].board_slot))            

    rdiff=rstart-rstop
    plt.hist(rdiff,20)


    rdiff_med=np.median(rdiff)
    binds=np.where(np.abs(rdiff-rdiff_med) > 0.5)[0]
    lowr_start_inds=np.intersect1d(np.where(rstart < 1.4)[0], np.where(rstop > 0.3)[0])
    highr_stop_inds=np.where(rstop > 0.8)[0]

    binds=np.unique(np.union1d(np.union1d(binds,lowr_start_inds),highr_stop_inds))
    if bololist=='all':
        binds=range(len(bkeys))
    if len(bololist)>0:
        binds=range(len(bkeys))

    plt.figure()
    blist=[]
    for bind in binds:
        plt.plot(ts[bkeys[bind]])
        plt.xlabel('Time')
        plt.ylabel('Resistance ($\Omega$)')
        plt.title(bkeys[bind])
        plt.show()
        tmpin=raw_input('Bad (y), else press enter to continue: ')
        if tmpin=='y':
            blist.append(bkeys[bind])
        plt.clf()
        
    print np.array(blist)
    raw_data={'name':blist}
    df=pd.DataFrame(raw_data,columns=['name'])
    df.to_csv(outputfile,index=False,sep='\t',encoding='utf-8')
    return



def create_blist_from_csv(filein='',physical_name=True):
    df=pd.read_csv(filein,sep='\t')
    frames=list(core.G3File('/spt_data/bolodata/downsampled/userdata/62886453/nominal_online_cal.g3'))
    calframe=frames[0]
    search_keys=[]
    blist=[]
    if physical_name:
        for kk in range(len(df['wafer'])):
            search_keys.append(df['wafer'][kk]+'_'+df['physical_name'][kk])
            
        
        bolos_master=calframe['NominalBolometerProperties'].keys()
        for bb in bolos_master:
            if calframe['NominalBolometerProperties'][bb].physical_name in search_keys:
                blist.append(bb)

    else:
        for kk in range(len(df['name'])):
            blist.append(df['name'][kk])
        
    blist=np.array(blist)
    return blist
