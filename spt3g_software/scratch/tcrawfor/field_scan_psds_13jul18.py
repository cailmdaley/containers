from spt3g import core, dfmux, std_processing, todfilter
from spt3g.util import tctools
import numpy as np
import pickle
from scipy import ndimage

obsid_str = '47718888'

file2 = '/spt/data/bolodata/downsampled/ra0hdec-67.25/'+obsid_str+'/nominal_online_cal.g3'
f2 = core.G3File(file2)
bp = f2.next()['NominalBolometerProperties']
names = bp.keys()
nbolo = len(names)

# get noise data & convert to watts
file3 = '/spt/data/bolodata/downsampled/ra0hdec-67.25/'+obsid_str+'/0000.g3'
f3 = core.G3File(file3)
noisedict1 = {}
psd_dict2 = {}
freqs1 = np.arange(1024)/1024.*152.6/2.
for i in np.arange(nbolo):
    noisedict1[names[i]] = []
    psd_dict2[names[i]] = np.zeros(1024)
converter = dfmux.unittransforms.ConvertTimestreamUnits(Input='RawTimestreams_I')
i = 0
nfnoise = 20
while i < nfnoise:
    frame = f3.next()
    if frame.type is core.G3FrameType.Wiring:
        converter(frame)
    if frame.type is core.G3FrameType.Scan:
        if 'Turnaround' not in frame:
            converter(frame)
            filtered_data = todfilter.polyutils.poly_filter_g3_timestream_map(frame['CalTimestreams'], 4)
#        psd_dict_temp, freqs_temp =  todfilter.dftutils.get_psd_of_ts_map(frame['CalTimestreams'])
            psd_dict_temp, freqs_temp =  todfilter.dftutils.get_psd_of_ts_map(filtered_data)
            freqs_temp_hz = freqs_temp/core.G3Units.Hz
            for name in names:
                try:
                    #                noisedict1[name].append(frame['CalTimestreams'][name])
                    noisedict1[name].append(filtered_data[name])
                    psd_dict2[name] += np.interp(freqs1,freqs_temp_hz,psd_dict_temp[name])/np.float(nfnoise)
                except:
                    pass
            i += 1

noisedict = {}
for name in names:
    ntemp = noisedict1[name]
    npts_tot = 0
    for i in np.arange(len(ntemp)):
        npts_tot += len(ntemp[i])
    noisedict[name] = np.zeros(npts_tot)
    npts_running = 0
    for i in np.arange(len(ntemp)):
        thisnpts = len(ntemp[i])
        noisedict[name][npts_running:npts_running+thisnpts] = ntemp[i]
        npts_running += thisnpts

psd_dict = {}
for name in noisedict.keys():
    pdtemp = tctools.quick_pspec(noisedict[name],rate=152.6,npts_psd=1024)
    psd_dict[name] = pdtemp['psd']
freqs = pdtemp['freqs']

