from spt3g import core, dfmux
from spt3g.scratch.tcrawfor import tctools
import numpy as np
import pickle

# kludge to get conversion factors
file2 = '/spt_data/bolodata/fullrate/calibrator/2747056/0000.g3'
f1 = core.G3File(file2)
f1.next()
wframe = f1.next()
frame2 = f1.next()
converter = dfmux.unittransforms.ConvertTimestreamUnits(Input='RawTimestreams_I')
converter(wframe)
converter(frame2)
convdict = (converter.convfactors)[core.G3TimestreamUnits.Counts]

## get elnod data
#elnod_dir = pickle.load(open('/home/tcrawfor/elnod_dir_02feb17.pkl'))
#names = elnod_dir.keys()
#bolos2get = []
#for key in names:
#    if elnod_dir[key]['elnod_sn'] > 20.:
#        bolos2get.append(key)
#
### !!!
##bolos2get_orig = (np.asarray(bolos2get)).copy()
##bolos2get = bolos2get[0:10]
## !!!
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
#dir1 = '/buffer/bolodata/'
#fname1 = dir1 + '20170201_140202.g3'
#f1 = core.G3File(fname1)
#wframe = f1.next()
#wmap = wframe['WiringMap']
#bnames = np.asarray(wmap.keys())
#cframe = f1.next()
#bp = cframe['NominalBolometerProperties']
#
#cmap_dict = {}
#for key in bolos2get:
#    fkey = filter(None,key.split('.'))
#    wmk = wmap[key]
#    ilist = [wmk.board_serial, wmk.module, wmk.channel*2]
#    cmap_dict[key] = ilist
#names = cmap_dict.keys()
#
#data1 = []
#data2 = []
#data3 = []
#
#files = [fname1]
#for fname in files:
#    for frame in core.G3File(fname):
#        grab_data(frame, cmap_dict, data1 = data1, data2 = data2, data3 = data3)
#
#bdata = np.zeros([len(data3[0]),len(data3)])
#for i in np.arange(len(data3)):
#    bdata[:,i] = data3[i]
#
#npts_psd = 1024
#psds = np.zeros([len(names),npts_psd])
#for i in np.arange(len(names)):
#    thisbdata = bdata[i,:] - np.mean(bdata[i,:])
#    psdtemp = tctools.quick_pspec(thisbdata,npts_psd=npts_psd)
#    psds[i,:] = psdtemp['psd']
#
#freqs = psdtemp['freqs']

whwhite = np.where(np.logical_and(freqs > 13., freqs < 17.))
whitenoise_cts = np.zeros(len(names))
for i in np.arange(len(names)):
    psdtemp = psds[i,:]
    whitenoise_cts[i] = np.sqrt(np.mean(psdtemp[whwhite]**2))

psd_dict = {}
for i in np.arange(len(names)):
    name = names[i]
    tempdict = {}
    tempdict['psd_w_rthz'] = psds[i,:]*convdict[name]
    tempdict['whitenoise_w_rthz'] = whitenoise_cts[i]*convdict[name]
    tempdict['freqs'] = freqs
    psd_dict[name] = tempdict


