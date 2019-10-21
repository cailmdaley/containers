import numpy as np
import scipy
import pickle
from spt3g import core, std_processing, gcp
from spt3g.scratch.tcrawfor import tctools

def grab_data(frame, cmap_dict, data1 = [], data2 = [], data3 = []):

    if frame.type != core.G3FrameType.Timepoint:
        return
    try:
        bdata = []
        data1.append(frame['CalibratorOn'])
        data2.append(frame['EventHeader'])
        for key in cmap_dict:
            inds = cmap_dict[key]
            bdata.append(frame['DfMux'][inds[0]][inds[1]][inds[2]])
        data3.append(bdata)
    except:
        pass

files = ['/buffer/bolodata/20170205_192145.g3']

dir1 = '/poleanalysis/sptdaq/20170129_rawdata/'
fname0 = dir1 + '20170129_194237.g3'
f1 = core.G3File(fname0)
wframe = f1.next()
wmap = wframe['WiringMap']
cframe = f1.next()
bp = cframe['NominalBolometerProperties']

bolos2get = pickle.load(open('gnames2_06feb17.pkl'))

cmap_dict = {}
for key in bolos2get:
    wmk = wmap[key]
    ilist = [wmk.board_serial, wmk.module, wmk.channel*2]
    cmap_dict[key] = ilist
names = cmap_dict.keys()

data1 = []
data2 = []
data3 = []

for fname in files:
    for frame in core.G3File(fname):
        grab_data(frame, cmap_dict, data1 = data1, data2 = data2, data3 = data3)

bdata = np.zeros([len(data3[0]),len(data3)])
for i in np.arange(len(data3)):
    bdata[:,i] = data3[i]

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
az = np.zeros(len(data4)*100)
el = np.zeros(len(data4)*100)
lst = np.zeros(len(data4)*100)
elmjds = np.zeros(len(data4)*100)
for i in np.arange(len(data4)):
    az[100*i:100*(i+1)] = data4[i]['tracker']['actual'][0]
    el[100*i:100*(i+1)] = data4[i]['tracker']['actual'][1]
    lst[100*i:100*(i+1)] = data4[i]['tracker']['lst'][0]
    for j in np.arange(100):
        elmjds[100*i+j] = (data4[i]['tracker']['utc'][0][j]).mjd

az_interp = np.interp(mjds-elmjds[0], elmjds-elmjds[0], az)
el_interp = np.interp(mjds-elmjds[0], elmjds-elmjds[0], el)
lst_interp = np.interp(mjds-elmjds[0], elmjds-elmjds[0], lst)
lst_hr = lst_interp/core.G3Units.h*60. # I think lst is wrong in cal file
lst_deg = lst_hr*15.
ra = np.mod(az_interp/core.G3Units.deg+lst_deg,360.)
dec = -el_interp/core.G3Units.deg

