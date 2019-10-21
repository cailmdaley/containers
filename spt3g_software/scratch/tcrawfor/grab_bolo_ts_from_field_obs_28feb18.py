from spt3g import core, dfmux, std_processing, todfilter
from spt3g.util import tctools
import numpy as np
import pickle, glob
from scipy import ndimage

obsid_str = '34473303'
#obsid_str = '34596211'


# get calframe and find best 100 bolos in each band
cf1 = core.G3File('/spt/user/production/calibration/calframe/ra0hdec-57.5/'+obsid_str+'.g3').next()
gb = {}
bands = ['90','150','220']
file2 = '/spt/data/bolodata/downsampled/RCW38-pixelraster/'+obsid_str+'/nominal_online_cal.g3'
bp = core.G3File('/spt/data/bolodata/downsampled/ra0hdec-57.5/'+obsid_str+'/nominal_online_cal.g3').next()['NominalBolometerProperties']
for band in bands:
    sntemp = np.asarray([cf1['CalibratorResponseSN'][name] for name in cf1['CalibratorResponseSN'].keys() if str(np.int(bp[name].band/10.)) == band and name != '005.5.1.3.2698'])
    bnames = np.asarray([name for name in cf1['CalibratorResponseSN'].keys() if str(np.int(bp[name].band/10.)) == band and name != '005.5.1.3.2698'])
    ssntemp = np.argsort(sntemp)
    gb[band] = bnames[ssntemp[-10:]]

# get noise data & convert to watts
files3 = glob.glob('/spt/data/bolodata/downsampled/ra0hdec-57.5/'+obsid_str+'/0*.g3')
files3.sort()
noisedict1 = {}
noisedict3 = {}
for band in bands:
    dtemp = {}
    for name in gb[band]:
        dtemp[name] = []
    noisedict1[band] = dtemp
    noisedict3[band] = np.zeros([400,4800])
iscan = 0
for file3 in files3:
    converter = dfmux.unittransforms.ConvertTimestreamUnits(Input='RawTimestreams_I')
    f3 = core.G3File(file3)
    for frame in f3:
        if frame.type is core.G3FrameType.Wiring:
            converter(frame)
        if frame.type is core.G3FrameType.Scan:
            converter(frame)
            for band in bands:
                for name in noisedict1[band].keys():
                    noisedict1[band][name].append(frame['CalTimestreams'][name])
            if 'Turnaround' not in frame:
                for band in bands:
                    for name in noisedict1[band].keys():
                        noisedict3[band][iscan,:] += frame['CalTimestreams'][name][0:4800]
                iscan += 1

noisedict2 = {}
for band in bands:
    noisedict2[band] = {}
    for name in noisedict1[band].keys():
        noisedict2[band][name] = tctools.list_to_array(noisedict1[band][name])
