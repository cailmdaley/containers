from spt3g import core, dfmux, std_processing, gcp
#from spt3g.scratch.tcrawfor import tctools
from spt3g.util import tctools
import numpy as np
import pickle

converter = dfmux.unittransforms.ConvertTimestreamUnits(Input='RawTimestreams_I')
#dir1 = '/spt_data/bolodata/fullrate/calibrator/3355750/'
dir1 = '/spt/data/bolodata/downsampled/calibrator/3355750/'
file2 = dir1 + 'nominal_online_cal.g3'
f2 = core.G3File(file2)
bp = f2.next()['NominalBolometerProperties']
names = bp.keys()
file3 = dir1 + '0000.g3'
f3 = core.G3File(file3)
oframe = f3.next()
wframe = f3.next()
converter(wframe)
dframe = f3.next()
converter(dframe)

psd_dict = {}
npts_psd = 1024
for name in names:
    try:
        thisbdata = dframe['CalTimestreams'][name]
        thisbdata -= np.mean(thisbdata)
        psdtemp = tctools.quick_pspec(thisbdata,rate=76.3,npts_psd=npts_psd)
        psd_dict[name] = psdtemp['psd']
    except:
        pass

freqs = psdtemp['freqs']
wndict = {}
#frange = [2,8]
frange = [11,17]
whf = np.where(np.logical_and(freqs > frange[0], freqs < frange[1]))
for name in names:
    try:
        thispsd = psd_dict[name]
        wndict[name] = np.sqrt(np.mean(thispsd[whf[0]]**2))
    except: 
        pass
