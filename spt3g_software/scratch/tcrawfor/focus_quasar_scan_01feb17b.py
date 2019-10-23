import numpy
import scipy
import pickle
import glob
from spt3g import core, std_processing, gcp
from spt3g.scratch.tcrawfor import tctools

np = numpy

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

dir1 = '/poleanalysis/sptdaq/20170129_rawdata/'
fname0 = dir1 + '20170129_194237.g3'
f1 = core.G3File(fname0)
wframe = f1.next()
wmap = wframe['WiringMap']
bnames = np.asarray(wmap.keys())
cframe = f1.next()
bp = cframe['NominalBolometerProperties']

f2 = core.G3File('/home/tcrawfor/spt_code/spt3g_software/calibration/scripts/elnod_output_2668507.g3')
eframe = f2.next()
names = eframe['ElnodSlopes'].keys()
elnod_s = np.array([eframe['ElnodSlopes'][key] for key in names])
elnod_ss = np.array([eframe['ElnodSigmaSlopes'][key] for key in names])
elnod_ss[(np.where(elnod_ss == 0.))[0]] = 1e12
elnod_sn = elnod_s / elnod_ss
n2get = 10
stemp = np.argsort(elnod_sn)
rstemp = stemp[::-1]
bolos2get = (np.asarray(names))[rstemp[0:10]]

cmap_dict = {}
for key in bolos2get:
    key2 = key
    wmk = wmap[key2]
    ilist = [wmk.board_serial, wmk.module, wmk.channel*2]
    cmap_dict[key2] = ilist

dir2 = '/buffer/bolodata/'
files = glob.glob('/buffer/bolodata/20170131_2*')
files.sort()
files = files[92:142]

data1 = []
data2 = []
data3 = []

for fname in files:
    for frame in core.G3File(fname):
        grab_data(frame, cmap_dict, data1 = data1, data2 = data2, data3 = data3)

#bdata = np.zeros([len(data3[0]),len(data3)])
#for i in np.arange(len(data3)):
#    bdata[:,i] = data3[i]

