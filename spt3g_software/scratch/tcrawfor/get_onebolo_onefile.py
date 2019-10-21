import numpy as np
from spt3g import core, std_processing, gcp

def grab_data(frame, inds, data1 = []):

    if frame.type != core.G3FrameType.Timepoint:
        return
    try:
        data1.append(frame['DfMux'][inds[0]][inds[1]][inds[2]])
    except:
        pass

#bolo2get = 'W158/2017.W158.4.34.6.3360'
#file1 = ''

dir1 = '/poleanalysis/sptdaq/20170129_rawdata/'
fname0 = dir1 + '20170129_194237.g3'
f1 = core.G3File(fname0)
wframe = f1.next()
wmap = wframe['WiringMap']
bnames = np.asarray(wmap.keys())
cframe = f1.next()
bp = cframe['NominalBolometerProperties']

wmk = wmap[bolo2get]
ilist = [wmk.board_serial, wmk.module, wmk.channel*2]

data1 = []
for frame in core.G3File(file1):
    grab_data(frame, ilist, data1 = data1)

