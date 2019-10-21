import pickle
import os, pdb
import matplotlib.pyplot as plt
import numpy as np

def plot_noise_vs_zsquid(datafile):
    f=open(datafile)
    data=pickle.load(f)
    f.close()
    kys=data.keys()
    zsquid=np.zeros(len(kys))
    med_noise=np.zeros(len(kys))
    for ind,ky in enumerate(kys):
        if np.mod(ind,500)==0:
            print ind
        zsquid[ind]=data[ky]['squid_z']
        med_noise[ind]=data[ky]['noise']['median_noise']
        #units = pA/rt(Hz)
    plt.figure()
    plt.plot(zsquid,med_noise,'k.')
    plt.xlabel('SQUID Transimpedance ($\Omega$)')
    plt.ylabel('Median Noise (pA /$\sqrt{Hz}$')
    return zsquid,med_noise

def plot_noise_vs_biasfreq(datafile,color_by_comb=False,mezzs=''):
    f=open(datafile)
    data=pickle.load(f)
    f.close()
    kys=data.keys()
    freqs=[]
    noise=[]
    mezz=[]
    module=[]
    camp=[]
    for ky in data.keys():
        freqs.append(data[ky]['frequency'])
        noise.append(data[ky]['noise']['median_noise'])
        mezz.append(ky.split('/')[2])
        module.append(ky.split('/')[3])
        camp.append(data[ky]['carrier_amp'])
#    plt.figure()
    mezz=np.array(mezz)
    module=np.array(module)
    mezz_master=['1','2']
    mod_master=['1','2','3','4']
    color=['ko','ro','bo','go','k+']

    if color_by_comb:
        for mm in mezz_master:
            cind=0
            plt.subplot(2,1,np.int32(mm))
            for modmod in mod_master:
                ginds=np.intersect1d(np.where(mezz == mm)[0],np.where(module==modmod)[0])

                if len(ginds)>0:
                    plt.plot(np.array(freqs)[ginds]/1e6,np.array(noise)[ginds],color[cind],markeredgecolor=color[cind][0])
                cind+=1
            plt.xlabel('Bias Frequency (MHz)')
            plt.ylabel('Median Noise Level(pA/$\sqrt{Hz}$)')
            plt.grid('on')
            plt.ylim(0,400)
                
    else:
        ginds=np.intersect1d(np.where(np.array(freqs) >100)[0],np.where(np.array(camp)>0)[0])
        if len(mezzs)>0:
            ginds=np.intersect1d(ginds,np.where(np.array(mezz) == mezzs)[0])
        plt.plot(np.array(freqs)[ginds]/1e6,np.array(noise)[ginds],'ko')
    plt.xlabel('Bias Frequency (MHz)')
    plt.ylabel('Median Noise Level(pA/$\sqrt{Hz}$)')
    plt.grid('on')
    return freqs, noise

def compare_2_noise_runs(run1file,run2file,mezzs=''):
    
    f=open(run1file)
    data1=pickle.load(f)
    f.close()
    kys1=data1.keys()
    f=open(run2file)
    data2=pickle.load(f)
    f.close()
    kys2=data2.keys()
    kys=np.intersect1d(np.array(kys1),np.array(kys2))
    freqs=[]
    noise1=[]
    noise2=[]
    camp=[]
    mezz=[]
    for ky in kys:
        noise1.append(data1[ky]['noise']['median_noise'])
        noise2.append(data2[ky]['noise']['median_noise'])
        freqs.append(data1[ky]['frequency'])
        camp.append(data1[ky]['carrier_amp'])
        mezz.append(ky.split('/')[2])
    ginds=np.intersect1d(np.where(np.array(freqs) >100)[0],np.where(np.array(camp)>0)[0])

    if len(mezzs)>0:
        ginds=np.intersect1d(ginds,np.where(np.array(mezz) == mezzs)[0])
#    pdb.set_trace()
    plt.figure()
    plt.plot(np.array(noise2)[ginds],np.array(noise1)[ginds],'ko')
    plt.ylabel('Median Noise from '+run1file+ '(pA/rt(Hz)')
    plt.xlabel('Median noise from '+run2file+ '(pA/rt(Hz)')
    xlm=plt.xlim()

    plt.plot(np.array([xlm[0],xlm[1]]),np.array([xlm[0],xlm[1]]),'r--')
    plt.grid('on')
    return

def compare_2_noise_runs_hist(run1file,run2file):
        
    f=open(run1file)
    data1=pickle.load(f)
    f.close()
    kys1=data1.keys()
    f=open(run2file)
    data2=pickle.load(f)
    f.close()
    kys2=data2.keys()
    kys=np.intersect1d(np.array(kys1),np.array(kys2))
    freqs=[]
    noise1=[]
    noise2=[]
    camp=[]
    for ky in kys:
        noise1.append(data1[ky]['noise']['median_noise'])
        noise2.append(data2[ky]['noise']['median_noise'])
        freqs.append(data1[ky]['frequency'])
        camp.append(data1[ky]['carrier_amp'])

    ginds=np.intersect1d(np.where(np.array(freqs) >100)[0],np.where(np.array(camp)>0)[0])

    bins=np.arange(0,800,10)
    plt.hist(noise1,bins,color='k',histtype='step',linewidth=2)
    plt.hist(noise2,bins,color='r',histtype='step',linewidth=2)
    plt.xlabel('Noise Level (pA/$\sqrt{Hz}$)')
    plt.ylabel('Number of Bolometers')
    return
