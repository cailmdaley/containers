from spt3g import core, dfmux, std_processing
#from spt3g.scratch.tcrawfor import tctools
from spt3g.util import tctools
import numpy as np
import pickle
from scipy import ndimage

obsid_str = '15799409'

f4 = core.G3File('/spt/data/bolodata/downsampled/RCW38-pixelraster/'+obsid_str+'/offline_calibration.g3')
cframe = f4.next()
names = []
for key in cframe['CalibratorResponseSN'].keys():
    if cframe['CalibratorResponseSN'][key] > 10. and cframe['ElnodSNSlopes'][key] < -30.:
        names.append(key)
nbolo = len(names)

# !!!
names = []
file1='/spt/user/benderan/noise_test/rcw38_maps/rcw38_15799409_full.g3'
f1 = core.G3File(file1)
for frame in f1:
    if frame.type is core.G3FrameType.Map:
        if frame['Id'] != 'bsmap':
            names.append(frame['Id'])
names = np.asarray(names)
nbolo = len(names)
# !!!


file2 = '/spt/data/bolodata/downsampled/RCW38-pixelraster/'+obsid_str+'/nominal_online_cal.g3'
f2 = core.G3File(file2)
bp = f2.next()['NominalBolometerProperties']
#names = bp.keys()
##names = names[0:2500]
#nbolo = len(names)

# get noise data & convert to watts
file3 = '/spt/data/bolodata/downsampled/RCW38-pixelraster/'+obsid_str+'/0000.g3'
f3 = core.G3File(file3)
noisedict1 = {}
for i in np.arange(nbolo):
    noisedict1[names[i]] = []
converter = dfmux.unittransforms.ConvertTimestreamUnits(Input='RawTimestreams_I')
i = 0
nfnoise = 20
while i < nfnoise:
    frame = f3.next()
    if frame.type is core.G3FrameType.Wiring:
        converter(frame)
    if frame.type is core.G3FrameType.Scan:
        converter(frame)
        for name in names:
            try:
                noisedict1[name].append(frame['CalTimestreams'][name])
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

# W/K for a perfect single-pol bolo with assumed bandwidths
kb = 1.38e-23
wpk_perfect = {}
wpk_perfect[90] = kb*22e9
wpk_perfect[150] = kb*38e9
wpk_perfect[220] = kb*48e9

# integrate maps & do NET calc
netdict = {}
opteffdict = {}
wndict = {}
peakdict = {}
mapdict = {}
#file1 = '/spt/user/production/calibration/RCW38-pixelraster/maps/'+obsid_str+'.g3'
f1 = core.G3File(file1)
for frame in f1:
    if frame.type is core.G3FrameType.Map:
        if 'Wunpol' in frame and frame['Id'] != 'bs':
            twt = np.asarray(frame['Wunpol'].TT)
            twt[np.where(twt == 0.)] = 1.
        if frame['Id'] in names:
            map_weighted = np.asarray(frame['T'])
            if np.max(np.abs(map_weighted)) > 0.:
                map_unw = map_weighted.copy()
                map_unw /= twt
# kludge!!!
                map_unw[:,300:] = 0.
                map_unw[:,0:50] = 0.
                map_unw[190:,:] = 0.
                map_unw[0:10,:] = 0.
# kludge!!!
                band = np.int(bp[frame['Id']].band/10.)
                net, w_per_k, whitelevel = \
                    tctools.quick_net_from_rcw38(map_unw,noisedict[frame['Id']], 
                                                 invert_map=True,band=band,rate=76.3)
                netdict[frame['Id']] = net
                opteffdict[frame['Id']] = w_per_k/wpk_perfect[band]
                wndict[frame['Id']] = whitelevel
                amap_sm = ndimage.gaussian_filter(-map_unw,4)
                ycenter, xcenter = np.unravel_index(np.argmax(amap_sm),[360,360])
                peakdict[frame['Id']] = map_unw[ycenter,xcenter]
                pdtemp = tctools.quick_pspec(noisedict[frame['Id']],rate=76.3)
                psd_dict[frame['Id']] = pdtemp['psd']
                freqs = pdtemp['freqs']
                mapdict[frame['Id']] = map_unw

netnames = netdict.keys()
bands = np.asarray([np.int(bp[name].band/10.) for name in netnames])
nets = np.asarray([netdict[name] for name in netnames])
opteffs = np.asarray([opteffdict[name] for name in netnames])
wlevels = np.asarray([wndict[name] for name in netnames])
peaks = np.asarray([peakdict[name] for name in netnames])


