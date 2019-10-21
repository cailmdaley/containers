#import numpy as np
#import scipy
#from scipy import ndimage
#import pickle
#import argparse as ap
#from spt3g import core, std_processing, gcp
#import os
#
#def grab_data(frame, cmap_dict, data1 = [], data2 = [], data3 = []):
#
#    if frame.type != core.G3FrameType.Timepoint:
#        return
#    try:
#        bdata = []
#        data1.append(frame['CalibratorOn'])
#        data2.append(frame['EventHeader'])
#        for key in cmap_dict:
#            inds = cmap_dict[key]
#            bdata.append(frame['DfMux'][inds[0]][inds[1]][inds[2]])
#        data3.append(bdata)
#    except:
#        pass
#
##files = ['/buffer/bolodata/20170204_161934.g3']
#files = ['/buffer/bolodata/20170205_183856.g3','/buffer/bolodata/20170205_183946.g3']
##'/buffer/bolodata/20170204_200355.g3','/buffer/bolodata/20170204_200445.g3']
#
#dir1 = '/poleanalysis/sptdaq/20170129_rawdata/'
#fname0 = dir1 + '20170129_194237.g3'
#f1 = core.G3File(fname0)
#wframe = f1.next()
#wmap = wframe['WiringMap']
#bnames = np.asarray(wmap.keys())
#cframe = f1.next()
#bp = cframe['NominalBolometerProperties']
#
## !!! this gets all bolos
#bolos2get = bnames
##bolos2get = bnames[67:72]
#
#cmap_dict = {}
#for key in bolos2get:
#    wmk = wmap[key]
#    ilist = [wmk.board_serial, wmk.module, wmk.channel*2]
#    cmap_dict[key] = ilist
#
##files = args.input
#
#data1 = []
#data2 = []
#data3 = []
#
#for fname in files:
#    for frame in core.G3File(fname):
#        grab_data(frame, cmap_dict, data1 = data1, data2 = data2, data3 = data3)
#
#bdata = np.zeros([len(data3[0]),len(data3)])
#for i in np.arange(len(data3)):
#    bdata[:,i] = data3[i]
#
#mjds = [ttime.mjd for ttime in data2]
#start_time = data2[0]
#stop_time = data2[len(data2)-1]
#
#def grab_boards(f, boardlist = []):
#    key = 'antenna0'
#    try:
#        boardlist.append(f[key])
#    except:
#        pass
#
#arcdir='/spt_data/arc/'
#data4 = []
#pipe1 = core.G3Pipeline()
#pipe1.Add(std_processing.ARCTimerangeReader, start_time=start_time, stop_time=stop_time, basedir=arcdir)
#pipe1.Add(gcp.ARCExtract)
#pipe1.Add(grab_boards, boardlist=data4)
#pipe1.Run()
#el = np.zeros(len(data4)*100)
#elmjds = np.zeros(len(data4)*100)
#for i in np.arange(len(data4)):
#    el[100*i:100*(i+1)] = data4[i]['tracker']['actual'][1]
#    for j in np.arange(100):
#        elmjds[100*i+j] = (data4[i]['tracker']['utc'][0][j]).mjd
#
#el_interp = np.interp(mjds-elmjds[0], elmjds-elmjds[0], el)
#template = 1./np.sin(el_interp)
#template -= np.mean(template)
#
## find actual elnod
#elsm = ndimage.gaussian_filter1d(el_interp/core.G3Units.deg,20)
#d_el = np.diff(elsm)
#dtime = (np.max(mjds) - np.min(mjds))*86400./np.float(len(template))
#deldt = d_el/dtime
#whelnod = np.where(np.abs(deldt) > 0.01)
#minwh = np.min(whelnod[0])
#maxwh = np.max(whelnod[0])
template = template[minwh:maxwh]
bdata = bdata[:,minwh:maxwh]

offsets = np.zeros(len(bolos2get))
elnod_response = np.zeros(len(bolos2get))
resid_rms = np.zeros(len(bolos2get))
elnod_sn = np.zeros(len(bolos2get))
for i in np.arange(len(bolos2get)):
#for i in np.arange(26):
    if np.sum(bdata[i,:]) == 0.:
        continue
    slope, intercept, r_value, p_value, std_err = \
        scipy.stats.linregress(template,bdata[i,:])
    offsets[i] = intercept
    elnod_response[i] = slope
    model = template*slope + intercept
    resids = bdata[i,:] - model
    resid_rms[i] = np.std(resids)
    if np.abs(std_err) > 0.:
        elnod_sn[i] = np.abs(slope/std_err)
#    elnod_sn[i] = np.sqrt(np.sum((template*slope/resid_rms[i])**2))

names = cmap_dict.keys()
elnod_dir = {}
for i in np.arange(len(names)):
    dirtemp = {}
    dirtemp['elnod_response'] = elnod_response[i]
    dirtemp['elnod_sn'] = elnod_sn[i]
    dirtemp['physical_name'] = bp[names[i]]
    elnod_dir[names[i]] = dirtemp


