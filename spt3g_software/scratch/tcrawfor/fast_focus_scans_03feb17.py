#import numpy as np
#import scipy
#import pickle
#from spt3g import core, std_processing, gcp
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
#elnod_dir = pickle.load(open('/home/tcrawfor/elnod_dir_03feb17.pkl'))
#bolos2get = []
#for key in elnod_dir.keys():
#    if np.isfinite(elnod_dir[key]['elnod_sn']) and elnod_dir[key]['elnod_sn'] > 3000.:
#        bolos2get.append(key)
#bolos2get = np.asarray(bolos2get)
#names = bolos2get.copy()
#
#cmap_dict = {}
#for key in bolos2get:
#    key2 = key
#    wmk = wmap[key2]
#    ilist = [wmk.board_serial, wmk.module, wmk.channel*2]
#    cmap_dict[key2] = ilist
#
#dir2 = '/buffer/bolodata/'
##fname1 = dir1 + '20170129_194920.g3'
##fname1 = dir2 + '20170131_190445.g3'
##fname2 = dir2 + '20170131_190535.g3'
##fname1 = dir2 + '20170131_210816.g3'
##fname2 = dir2 + '20170131_210906.g3'
##files = [fname1, fname2]
##files = [fname1]
#files = [dir2 + '20170202_155143.g3',dir2 + '20170202_155234.g3',dir2 + '20170202_155324.g3',dir2 + '20170202_155415.g3']
#
#data1 = []
#data2 = []
#data3 = []
#
#for fname in files:
#    for frame in core.G3File(fname):
#        grab_data(frame, cmap_dict, data1 = data1, data2 = data2, data3 = data3)
#
bdata = np.zeros([len(data3[0]),len(data3)])
for i in np.arange(len(data3)):
    bdata[:,i] = data3[i]

cnames = np.asarray(cmap_dict.keys())
bdata_orig = bdata.copy()
for i in np.arange(len(names)):
    whn = np.where(cnames == names[i])
    j = whn[0][0]
    bdata[i,:] = bdata_orig[j,:] - np.mean(bdata_orig[j,:])

mjds = [ttime.mjd for ttime in data2]
start_time = data2[0]
stop_time = data2[len(data2)-1]

def grab_boards(f, boardlist = []):
    key = 'antenna0'
    try:
        boardlist.append(f[key])
    except:
        pass

arcdir='/spt_data/arc/'
data4 = []
pipe1 = core.G3Pipeline()
pipe1.Add(std_processing.ARCTimerangeReader, start_time=start_time, stop_time=stop_time, basedir=arcdir)
pipe1.Add(gcp.ARCExtract)
pipe1.Add(grab_boards, boardlist=data4)
pipe1.Run()
el = np.zeros(len(data4)*100)
az = np.zeros(len(data4)*100)
lst = np.zeros(len(data4)*100)
elmjds = np.zeros(len(data4)*100)
for i in np.arange(len(data4)):
    el[100*i:100*(i+1)] = data4[i]['tracker']['actual'][1]
    az[100*i:100*(i+1)] = data4[i]['tracker']['actual'][0]
    lst[100*i:100*(i+1)] = data4[i]['tracker']['lst'][0]
    for j in np.arange(100):
        elmjds[100*i+j] = (data4[i]['tracker']['utc'][0][j]).mjd

el_interp = np.interp(mjds-elmjds[0], elmjds-elmjds[0], el)
az_interp = np.interp(mjds-elmjds[0], elmjds-elmjds[0], az)
lst_interp = np.interp(mjds-elmjds[0], elmjds-elmjds[0], lst)
lst_interp_hr = lst_interp/core.G3Units.h*60. # I think lst is wrong cal file
lst_interp_deg = lst_interp_hr*15.
ra_interp = np.mod(az_interp/core.G3Units.deg+lst_interp_deg,360.)
dec_interp = -el_interp/core.G3Units.deg
#
#radec0 = [264.2,-21.6]
#reso_arcmin = 0.5
#xoff = (ra_interp - radec0[0])*np.cos(el_interp)*60./reso_arcmin
#yoff = (dec_interp - radec0[1])*60./reso_arcmin
#mapshape = [400,400]
#xpix = np.round(xoff+mapshape[1]/2.)
#ypix = np.round(yoff+mapshape[0]/2.)
#
#maps = np.zeros([len(names),mapshape[0],mapshape[1]])
#hits = np.zeros(mapshape)
#for i in np.arange(len(xpix)):
#    hits[ypix[i],xpix[i]] += 1.
#    for j in np.arange(len(names)):
#        maps[j,ypix[i],xpix[i]] += bdata[j,i]
