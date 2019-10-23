from spt3g import core, dfmux, std_processing, gcp
from spt3g.scratch.tcrawfor import tctools
import numpy as np
import pickle

file2 = '/spt_data/bolodata/fullrate/RCW38/2653414/nominal_online_cal.g3'
f2 = core.G3File(file2)
bp = f2.next()['NominalBolometerProperties']

#names = pickle.load(open('goodnames_03feb17.pkl'))
#nbolo = len(names)
file1 = '/poleanalysis/cal_results/calibrator/3090145.g3'
f1 = core.G3File(file1)
cframe = f1.next()
file1 = '/poleanalysis/cal_results/elnod/3091165.g3'
f2 = core.G3File(file1)
enframe = f2.next()
names = []
for name in cframe['CalibratorResponseSN'].keys():
    ert = enframe['ElnodSlopes'][name]
    est = enframe['ElnodSigmaSlopes'][name]
    esnt = 0.
    if est > 0.:
        esnt = np.abs(ert/est)
    if cframe['CalibratorResponseSN'][name] > 10. and esnt > 300.:
        names.append(name)

converter = dfmux.unittransforms.ConvertTimestreamUnits(Input='RawTimestreams_I')
file3 = '/spt_data/bolodata/fullrate/calibrator/3090937/0000.g3'
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
#frange = [2,8]
frange = [11,19]
whf = np.where(np.logical_and(freqs > frange[0], freqs < frange[1]))
for name in names:
    thispsd = psd_dict[name]
    wndict[name] = np.sqrt(np.mean(thispsd[whf[0]]**2))
