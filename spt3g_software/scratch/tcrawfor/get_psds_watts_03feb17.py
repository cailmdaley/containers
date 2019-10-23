from spt3g import core, dfmux
from spt3g.scratch.tcrawfor import tctools
import numpy as np
import pickle

file2 = '/spt_data/bolodata/fullrate/RCW38/2653414/nominal_online_cal.g3'
f2 = core.G3File(file2)
bp = f2.next()['NominalBolometerProperties']

names = pickle.load(open('goodnames_03feb17.pkl'))
nbolo = len(names)

converter = dfmux.unittransforms.ConvertTimestreamUnits(Input='RawTimestreams_I')
file3 = '/spt_data/bolodata/fullrate/calibrator/2815395/0000.g3'
f3 = core.G3File(file3)
oframe = f3.next()
wframe = f3.next()
converter(wframe)
dframe = f3.next()
converter(dframe)

psd_dict = {}
npts_psd = 1024
for name in names:
    thisbdata = dframe['CalTimestreams'][name]
    thisbdata -= np.mean(thisbdata)
    psdtemp = tctools.quick_pspec(thisbdata,npts_psd=npts_psd)
    psd_dict[name] = psdtemp['psd']

freqs = psdtemp['freqs']
wndict = {}
frange = [2,8]
whf = np.where(np.logical_and(freqs > frange[0], freqs < frange[1]))
for name in names:
    thispsd = psd_dict[name]
    wndict[name] = np.sqrt(np.mean(thispsd[whf[0]]**2))
