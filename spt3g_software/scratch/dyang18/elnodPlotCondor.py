import matplotlib
matplotlib.use('Agg')
import argparse as ap
from spt3g import core, std_processing, mapmaker, dfmux, calibration, xtalk, todfilter, coordinateutils
import scipy.stats
import os, numpy
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib import cm
import pandas as pd
from spt3g.std_processing import obsid_to_g3time
from glob import glob
import tarfile
from astropy.time import Time

# Pass in the needed parameters
P = ap.ArgumentParser(description='Various Elnod Plots',
              formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('-o','--obsid', action='store', default=None,
           help='observation id')
P.add_argument('-d','--daterange', action='store', default=None,
           help='whether this is a timeseries plot')
args = P.parse_args()

#extract all the observation files
tar = tarfile.open('elnod.tar.gz', "r:gz")
tar.extractall()
tar.close()
tar = tarfile.open('elnodBolometerProperties.tar.gz', "r:gz")
tar.extractall()
tar.close()

#this need to be changed depending on the name of the file on condor
obsp=['ElnodDrifts', 'ElnodRSquared', 'ElnodSNSlopes', 'ElnodSigmaSlopes','ElnodSlopes','ElnodVariances']
wafer_lst=['W136', 'W139', 'W142', 'W147', 'W148', 'W152', 'W153', 'W157', 'W158', 'W162']

# Helper functions for stats things
def autoRange(a, nSigma=3.5,remove_percent=3.0):
    '''
    comments leftover from copy pasting
    Throw out top and bottom remove_percent% of outliers
    Return (vmin, vmax, ticks) were
    (vmin, vmax) ~= (mean-nSigma*sigma, mean+nSigma*sigma)
    unless either side of this exceeds the actual min/max,
    in which case return that.
    The approximatly-equal uses matplotlib's ticker.AutoLocator to round to nice values, 
    since the histogram will extend the limits to these anyway
    NOTE: Unfortunately, vmin and vmax are sometimes slightly too big: 1.00000000000000008882
          Which is why clip() and the Histograming part of makePlots() must shrink this range by epsilon
    for bool or and int types with a small range, lock the min and max exactly to take full advantage of the colors.
    '''
    (mina, maxa, mean, sigma, median) = finiteStats(a,remove_percent)  # returns all nan if there are no finite

    low  = max(mina, mean-nSigma*sigma)
    high = min(maxa, mean+nSigma*sigma)
    return [low,high]
def finiteStats(a,remove_percent):
    '''
    Throw out top and bottom remove_percent% of outliers (unless boolean)
    return (min, max, mean, sigma) of the finite
    '''
    # Throw out nan
    
    if len(a)==0:
        return (np.nan, np.nan, np.nan, np.nan, np.nan)
    amid=a[ (a >= np.percentile(a,remove_percent)) * 
                (a <= np.percentile(a,100-remove_percent)) ]
    #amid  =[x for x in a if x >= np.percentile(a,remove_percent) and x <= np.percentile(a,100-remove_percent)]
    mina   = np.min(a)
    maxa   = np.max(a)
    mean   = np.mean(amid)  # Note empty arrays return nan
    sigma  = np.std(amid)
    median = np.median(amid)
    return (mina, maxa, mean, sigma, median)



def getcolor(n):
    colors = np.linspace(0,1,1000)
    color_list = cm.rainbow(colors).tolist()
    color=[]
    for i in color_list:
        color.append(tuple(i))
    if n<0:
        return (1.0,0,0)
    elif n>=999:
        return color[-1]
    else:
        return color[n]

def StackedHistogram(request):
 
    obsid=request['observation']
    varname=request['parameter'] #'CalibratorResponse'
    #this grabs the data and bolometer properties frame
    for fr in core.G3File(request['input']):
        data_frame=fr
    #need to be changed according to how I decide to treat bolometer files
    for fr in core.G3File(request['boloprop']):
        bolo_props=fr['NominalBolometerProperties']

    

    
    #temporarily set ith 3 frequencies because some bolos are giving negative
    #frequency which is causing issues
    #bands = np.unique([bolo_props[bname].band for bname in bolo_props.keys()])
    #bands.sort()
    bands=[900.0,1500.0,2200.0]
    wafers = np.unique([bname.split('/')[0] for bname in bolo_props.keys()])
    wafers.sort()
    
    
    data = {k:v for k,v in data_frame[varname].iteritems() if k in bolo_props.keys()}
    
    
    #establishing color conventions
    colors = np.linspace(0,1,len(bands)*len(wafers))
    color_list = cm.rainbow(colors).tolist()
    color=[]
    for i in color_list:
        color.append(tuple(i))
    color.reverse()
    #print(len(color))
    #using autoRange, assuming  data_frame[varname] is a dictionary keyed by bolometer names and contain
    # values
    plotrange=autoRange(np.array((data_frame[varname]).values()))
    # plot data by band
    f = plt.figure()
    bindata=[] 
    bindatakey=[]
    histbins = np.linspace(plotrange[0], plotrange[1], 101)
    for band in bands:
        for wafer in wafers:
            bindatakey.append(str(band/core.G3Units.GHz)+'/'+str(wafer))
            bindata.append([data[bname] for bname in data.keys() if bname.split('/')[0] == wafer and bolo_props[bname].band==band])
    plt.hist([np.clip(bindata[i],plotrange[0],plotrange[1]) for i in range(len(bindata))], bins=histbins, color=color, stacked =True, label=bindatakey)
    plt.legend(bbox_to_anchor=(0.9,1.05),loc=2,prop={'size':7})
    plt.xlim(plotrange[0],plotrange[1])
    f.savefig(obsid+'elnod'+'StackedHistogram'+'.png')

def FocalPlane(request):
    obsid=request['observation']
    varname=request['parameter'] #CalibratorResponse'
    #this grabs the right calframe
    for fr in core.G3File(request['input']):
        data_frame=fr
    # need to be changed according to how I treat bolometer property file.
    for fr in core.G3File(request['boloprop']):
        bolo_props=fr['NominalBolometerProperties']

                     
    bands = np.unique([bolo_props[bname].band for bname in bolo_props.keys()])
    #using autoRange, assuming  data_frame[varname] is a dictionary keyed by bolometer names and contain
    # values
    bands.sort()
    data = {k:v for k,v in (data_frame[varname]).iteritems() if k in bolo_props.keys()}
    for i in range(len(bands)):
     
        bindata=np.array([data[bname]for bname in data.keys() if bolo_props[bname].band==bands[i]])
        xdata=np.array([bolo_props[bname].x_offset for bname in data.keys() if bolo_props[bname].band==bands[i]])
        ydata=np.array([bolo_props[bname].y_offset for bname in data.keys() if bolo_props[bname].band==bands[i]])
        if len(xdata)>0 and len(ydata)>0:
            
            f, ax1 = plt.subplots(2,gridspec_kw = {'height_ratios':[10, 1]},figsize=(15,15))

            pltrange=autoRange(bindata)
            margin=[0.002,0.002]
            stepsize=(pltrange[1]-pltrange[0])/1000
            for j in range(len(xdata)):
                ax1[0].add_artist(matplotlib.patches.Rectangle(xy=(xdata[j],ydata[j]),width=0.0001,height=0.0005,color=getcolor(int((bindata[j]-pltrange[0])/stepsize))))
            ax1[0].set_xlim([np.nanmin(xdata)-margin[0],np.nanmax(xdata)+margin[0]])
            ax1[0].set_ylim([np.nanmin(ydata)-margin[1],np.nanmax(ydata)+margin[1]])
            cmap = matplotlib.cm.rainbow
            norm = matplotlib.colors.Normalize(vmin=pltrange[0], vmax=pltrange[1])
            cb = matplotlib.colorbar.ColorbarBase(ax1[1], cmap=cmap,
                                    norm=norm,
                                    orientation='horizontal')
            f.savefig(str(obsid)+'elnod'+str(bands[i]/core.G3Units.GHz)+'GHz'+'FocalPlane'+'.png')

def HistogramOverTime(request):
    bins=100
    obstype='elnod'
    obsparameter=request['parameter'] #CalibratorResponse'

    #end of required parameters
   
        # grab the right files and observation numbers
    filelst=request['input']
    obslist=request['observation']
    dater=request['daterange']
    if len(filelst)!=0:
        obslist.sort()
        filelst.sort()



        # obstain a list of all bolometers in a specific wafer            
        bolos=[]

        for fr in core.G3File(filelst):
            #this part is still hard coded in.
            bp=fr[obsparameter]
            bolos=bp.keys()
            #print([fr['CalibratorResponse'][bolo] for bolo in bolos if bolo in fr['CalibratorResponse']])
            break 

        frs=[]
        newobslist=[]
        #make an updated list of frames
        for i in range(len(filelst)):

            for fr in core.G3File(filelst[i]):
                if obsparameter in fr:
                    frs.append(fr[obsparameter].values())
                    newobslist.append(obslist[i])
                else: print('Thing not in frame')
        frs=pd.DataFrame(frs)
        frs=frs[frs!=0]
        frs=frs.as_matrix()
        frs=np.absolute(frs)

        flatfrs=frs.flatten()
        flatfrs=flatfrs[np.logical_not(np.isnan(flatfrs))]

        vmin,vmax=autoRange(flatfrs)

        histo_edges=[np.histogram(np.clip(i,vmin,vmax),bins=bins,range=(vmin,vmax)) for i in frs]
        histo_lines = [he[0] for he in histo_edges]
        histomax = np.percentile([max(h) for h in histo_lines], 50)
        step=(vmax-vmin)/bins
        f=plt.figure(figsize=(10,20))
        plt.imshow(histo_lines,interpolation='none',vmin=vmin, vmax=histomax)
        plt.yticks([int(i) for i in np.linspace(0,len(frs)-1,num=30)],[str(obsid_to_g3time(newobslist[j]))[:17] for j in [int(i) for i in np.linspace(0,len(frs)-1,num=30)]],fontsize=5)
        plt.xticks(np.linspace(0,100,5),[int(step*i) for i in np.linspace(0,100,5)],rotation='vertical',fontsize=5)
        plt.colorbar(fraction=0.02, pad=0.02, aspect=40)
        plt.xlabel(obsparameter)
        f.savefig(str(dater)+'elnod'+str(obsparameter)+'HistogramOverTime.png')



def SingleWaferArray(request):
    obstype='elnod'

    obsparameter=request['parameter'] 
    
    filelst=request['input']
    obslist=request['observation']
    dater=request['daterange']
    if len(filelst)!=0:
        obslist.sort()
        filelst.sort()

        frs=[]
        newobslist=[]
        #make an updated list of frames
        for i in range(len(filelst)):
            d=core.G3File(filelst[i])
            for fr in d:
                if obsparameter in fr:
                    frs.append(fr[obsparameter])
                    newobslist.append(obslist[i])
                
        for wafer in wafer_lst:
        # obstain a list of all bolometers in a specific wafer            
            bolos=[]
            for fr in core.G3File(filelst):

                bp=fr[obsparameter]
                bolos=[key for key in bp.keys() if (key.split('/'))[0]==wafer]
                
                break 



        # make and fill in data
            arr=numpy.zeros(shape = (len(frs), len(bolos)))
            for i in range(len(arr)):
                for j in range(len(arr[0])):
                    if bolos[j] in frs[i]:
                        arr[i][j]=(frs[i])[bolos[j]]
                    if not np.isfinite(arr[i][j]):
                        arr[i][j]=0
            if len(arr)!=0:
                f=plt.figure(figsize=(15,10))
                plt.imshow(arr,interpolation='none')
                plt.yticks([int(i) for i in np.linspace(0,len(frs)-1,num=30)],[str(obsid_to_g3time(newobslist[j])) for j in [int(i) for i in np.linspace(0,len(frs)-1,num=30)]],fontsize=5)

                plt.xticks([int(i) for i in np.linspace(0,len(bolos)-1,num=100)],[bolos[j] for j in [int(i) for i in np.linspace(0,len(bolos)-1,num=100)]],rotation='vertical',fontsize=5)

                plt.colorbar(fraction=0.03)
                plt.savefig(str(dater)+'elnod'+wafer+str(obsparameter)+'Array.png')

if args.daterange is not None:

    input_files=glob('*.g3')
    obslist=[]
    for fil in input_files:
        a=fil.split('B')[0]
        obsid=int(a.split('.')[0])
            
            
        x=obsid_to_g3time(obsid)
        if str(args.daterange) in Time(x.mjd,format='mjd').iso:
            obslist.append(obsid)
            
    obslist.sort()
    filelst=[str(obsid)+'.g3' for obsid in obslist]
    filelst = [fil for fil in filelst if os.path.isfile(fil)]
    for p in obsp:
        
        r={}
        r['observation']=obslist
        r['input']=filelst
        r['daterange']=args.daterange
        r['parameter']=p
    # need to convert boloprop file to a new format
        HistogramOverTime(r)
        SingleWaferArray(r)
else:
        # otherwise we try a bunch of different dates
    input_files=glob('*.g3')
    t=[str(x)+'-'+y for x in [2017,2018,2019] for y in ['01','02','03','04','05','06','07','08','09','10','11','12']]
    for time in t:
        obslist=[]
        for fil in input_files:
            a=fil.split('B')[0]
            obsid=int(a.split('.')[0])
            
            x=obsid_to_g3time(obsid)
            if time in Time(x.mjd,format='mjd').iso:
                obslist.append(obsid)

        obslist.sort()
        filelst=[str(obsid)+'.g3' for obsid in obslist]
        filelst = [fil for fil in filelst if os.path.isfile(fil)]
        for p in obsp:
            
            r={}
            r['observation']=obslist
            r['input']=filelst
            r['daterange']=time
            r['parameter']=p
        # need to convert boloprop file to a new format
            HistogramOverTime(r)
            SingleWaferArray(r)

input_files=glob('*.g3')
obslist=[int((fil.split('B')[0]).split('.')[0]) for fil in input_files]
if args.obsid is None:
    for obs in obslist:
        for p in obsp:
            r={}
            r['observation']=str(obs)
            r['input']=str(obs)+'.g3'
            r['parameter']=p
            r['boloprop']=str(obs)+'BolometerProperties.g3'
            StackedHistogram(r)
            FocalPlane(r)
else:
    for p in obsp:
        
        r={}
        r['observation']=str(args.obsid)
        r['input']=str(args.obsid)+'.g3'
        r['boloprop']=str(args.obsid)+'BolometerProperties.g3'
        r['parameter']=p
        StackedHistogram(r)
        FocalPlane(r)
    
#Now tar the outputfiles
outputfilelst=glob('*.png')
tar = tarfile.open(args.daterange+'elnod'+args.obsid +".tar.gz", "w:gz")
for fil in outputfilelst:
    tar.add(fil)
tar.close()
