from movingAvg import movingAvg
from spt3g import core, gcp
from spt3g import std_processing
from spt3g.util.extractdata import extract_keys, Accumulator
import os
import datetime
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import pydfmux

#stolen from ndhutils, now with 5000% more comments
def arcfiles_for_time(start, stop, arcdir = '/spt_data/arc'):
    ''' 
    returns sorted list of strings containing paths to arcfiles between start and stop.  

    start=time to start looking for arcfiles
    stop = time to stop looking for arfciles
    arcdir = (string) base directory where aout=rcfiles are stored.  Defaults to /spt_data/arc

    time formatting: (G3Time object) example: start=core.G3Time( '20180101_125601')
    to get seconds: (start.time)/core.G3Units.s

    '''

    #minifunction takes path and returns the last/bottom filename or directory, with the .<whatever> stripped off.  e.g. /home/lowitz/foo.py ==> foo.  In this case, returns the timestamp of the file, then turns it into a g3 time object
    arc2time = lambda path: core.G3Time(os.path.basename(path).split('.')[0])

    #returns full path to first arcfile before specified time
    arcfiles = [std_processing.FindARCFiles.ARCForTime(start, arcdir)]

    #build list of files
    count=0
    err=False
    while arc2time(arcfiles[-1]).time < stop.time and err==False:
        try:
            nextarc=std_processing.FindARCFiles.NextARCFile(arcfiles[-1])
            count+=1
            arcfiles.append(nextarc)
            if np.mod(count,20)==0:
                print(nextarc)
        except:
            err=True
            continue
    return arcfiles[:-1]


def fridgeData(start, stop, ds=1, save=True, arcdir='/spt_data/arc', savepath='/home/lowitz/fridgeData'):
    #get sorted list of arcfiles
    arcs=arcfiles_for_time(start, stop, arcdir)
    print('Done making arcfile list.  Processing...')

    keys={'time': ['array', 'frame', 'utc'], 
          'UChead':['array', 'cryo', 'temperature', 0, 0], 
          'UCstage':['array', 'cryo', 'temperature', 0, 10], 
          'ICstage':['array', 'cryo', 'temperature', 0, 12],
          'IChead':['array', 'cryo', 'temperature', 0, 1],
          'He4head':['array', 'cryo', 'temperature', 0, 2],
          'UCpump':['array', 'cryo', 'temperature', 0, 6],
          'UCppwr':['array', 'cryo', 'heater_dac', 0, 5],
          'UCswitch':['array', 'cryo', 'temperature', 0, 9],
          'ICswitch':['array', 'cryo', 'temperture', 0, 8],
          'head50K':['array', 'cryo', 'temperature', 0, 15],
          'LCtower':['array', 'cryo', 'temperature', 0, 11]}

    reader=gcp.ARCFileReader(arcs)
    
    if isinstance(keys, dict):
        accums = {}
        for k in keys:
            accums[k] = Accumulator(keys[k])
    else:
        accums = {'value': Accumulator(keys)}

    pipe = core.G3Pipeline()
    pipe.Add(reader)
    pipe.Add(gcp.ARCExtract)
    for a in accums.values():
        pipe.Add(a)
    pipe.Run()
    out={}
    for k, acc in accums.items():
        out[k] = acc.extract_values()[::ds]
    if 'value' in out.keys() and len(out.keys()) == 1:
        out = out['value']
    if save:
        try:
            name=start.GetFileFormatString()+'_'+stop.GetFileFormatString()
            path=savepath+'/'+name+'.pkl'
            pydfmux.save_pickle(out,path)
        except:
            print('WARNING: saving failed')

    return out

#adapted from spt3g_software/scratch/ndhuang/find_el_glitches.py
def ptdata(start, stop, ds=100, save=True, arcdir='/spt_data/arc', savepath='/home/lowitz/pt415_data'):
    '''
    collects pulse tube compressor pressure data from GCP between start and stop into a dict.

    start: (G3Time object) time to start looking for arcfiles.  eg: spt3g.core.G3Time('20190101_000000')
    stop: (G3Time object) time to stop looking for arcfiles (will fail if too recent).  Same format as start
    ds= downsample rate. e.g. ds=100 downsamples by a factor of 100.  default 100
    save: (boolean) whether or not to save output
    arcdir: (string) where to look for arcfiles
    savepath: (string) where to save the output (not used if save==False)
    '''

    #get sorted list of arcfiles
    arcs=arcfiles_for_time(start, stop, arcdir)
    print('Done making arcfile list.  Processing...')

    keys={'time': ['array', 'frame', 'utc'], 
          'optH':['array', 'pt415', 'pressure_high', 0], 
          'optL':['array', 'pt415', 'pressure_low', 0], 
          'rxH':['array', 'pt415', 'pressure_high', 1], 
          'rxL':['array', 'pt415', 'pressure_low',1]}
    reader=gcp.ARCFileReader(arcs)

    if isinstance(keys, dict):
        accums = {}
        for k in keys:
            accums[k] = Accumulator(keys[k])
    else:
        accums = {'value': Accumulator(keys)}

    pipe = core.G3Pipeline()
    pipe.Add(reader)
    pipe.Add(gcp.ARCExtract)
    for a in accums.values():
        pipe.Add(a)
    pipe.Run()
    out={}
    for k, acc in accums.items():
        out[k] = acc.extract_values()[::ds]
    if 'value' in out.keys() and len(out.keys()) == 1:
        out = out['value']
    if save:
        try:
            name=start.GetFileFormatString()+'_'+stop.GetFileFormatString()
            path=savepath+'/'+name+'.pkl'
            pydfmux.save_pickle(out,path)
        except:
            print('WARNING: saving failed')

    return out
    

    #data=extract_keys(keys,arcs)

def ptplotter(out):
    '''
    takes output from ptdata and makes a nice plot
    '''
    
    #make times into a nice format
    times=[datetime.utcfromtimestamp(x.time / core.G3Units.s) for x in out['time']]

    plt.figure('high')
    plt.plot(times,out['optH'], 'r',label='Optical High')
    plt.plot(times, out['rxH'], 'g', label='Rx High')
    plt.legend()
    plt.tight_layout()

    plt.figure('low')
    plt.plot(times,out['optL'], 'r',label='Optical Low')
    plt.plot(times, out['rxL'], 'g',label='Rx Low')
    plt.legend()
    plt.tight_layout()

    plt.show()

def pwrPlot(data, key='UCppwr', align=.3, smooth=0):
    for run in data:
        
        startTime=data[run]['time'][0].time
        times=[]
        for t in data[run]['time']:
            times.append((t.time-startTime)/core.G3Units.s)

        temps=data[run][key]
        idx=[x for x,i in enumerate(temps) if i>align][0]

        temps_later=temps[9000:]
        idx=[x for x,i in enumerate(temps_later) if i>align][0]
        idx=idx-5600

        print(idx)

        alignedStartTime=times[idx]
        alignedTimes=[]
        for t in times[idx:]:
            alignedTimes.append((t-alignedStartTime)/60)
        
        
        temps_to_plot=temps[idx:]
        if smooth>0:
            temps_to_plot=movingAvg(temps_to_plot, smooth)
            
        plt.plot(alignedTimes,temps_to_plot, label=run)


    plt.title(key)
    plt.legend()
    plt.xlabel('Time [min]')
    plt.ylabel('power [mW]')
    plt.tight_layout()
    plt.show()



def FridgePlot(data, therm='UCstage', align=.3, smooth=0):

'''
Plots one thermometer for many fridge cycles
data is a dict of several outputs from fridgeData().  Typially, one dict entry per fridge cycle.
therm is which thermometer you want to plot
align:  all of the cycles will the lined up according to when therm *first* hits this temperature in K
smooth: if you want to smooth out noisy thermometry data, this is the smoothing window.  If smooth=0, no smoothing is applied
'''

    for run in data:
        
        startTime=data[run]['time'][0].time
        times=[]
        for t in data[run]['time']:
            times.append((t.time-startTime)/core.G3Units.s)

        temps=data[run][therm]
        idx=[x for x,i in enumerate(temps) if i>align][0]

        if therm=='UCppwr':
            temps_later=temps[9000:]
            idx=[x for x,i in enumerate(temps_later) if i>align][0]
            idx=idx-5690

        print(idx)

        alignedStartTime=times[idx]
        alignedTimes=[]
        for t in times[idx:]:
            alignedTimes.append((t-alignedStartTime)/60)
        
        
        temps_to_plot=temps[idx:]
        if smooth>0:
            temps_to_plot=movingAvg(temps_to_plot, smooth)
            
        plt.plot(alignedTimes,temps_to_plot, label=run)


    plt.title(therm)
    plt.legend()
    plt.xlabel('Time [min]')
    plt.ylabel('Temperature [K]')
    plt.tight_layout()
    plt.show()

def FridgePlotSingle(data, therms=['UCstage', 'ICstage', 'UChead', 'LCtower', 'UCswitch'], title=''):

'''
plots several thermometers from a single fridge cycle

data is a dict, output straight from fridgeData().
Therms is a list of thermometer names from the dict made in fridgeData().  These are the only ones that will be plotted. 
title= title for the plot
'''
    
    for therm in therms:
        startTime=data['time'][0].time
        print(data['time'][0].Description())
        times=[]
        for t in data['time']:
            times.append(((t.time-startTime)/core.G3Units.s)/60)
        temps=data[therm]


        plt.plot(times,temps, label=therm)

    plt.legend()
    plt.title(title)
    plt.xlabel('Time [min]')
    plt.ylabel('Temperature [K]')
    plt.tight_layout()
    plt.show()
        


    

    
