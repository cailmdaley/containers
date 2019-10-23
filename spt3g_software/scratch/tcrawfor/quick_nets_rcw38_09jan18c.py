from spt3g import core, dfmux, std_processing
#from spt3g.scratch.tcrawfor import tctools
from spt3g.util import tctools
import numpy as np
import pickle
from scipy import ndimage

obsid_str = '31980196'
#obsid_str = '21167083'
#obsid_str = '20029961'

file2 = '/spt/data/bolodata/downsampled/RCW38-pixelraster/'+obsid_str+'/nominal_online_cal.g3'
f2 = core.G3File(file2)
bp = f2.next()['NominalBolometerProperties']
names = bp.keys()
nbolo = len(names)

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
#wpk_perfect[90] = kb*22e9
#wpk_perfect[150] = kb*38e9
#wpk_perfect[220] = kb*48e9
wpk_perfect[90] = kb*28.2e9 # new from Brad
wpk_perfect[150] = kb*42.6e9
wpk_perfect[220] = kb*51.9e9

# integrate maps & do NET calc
netdict = {}
netdict2 = {}
opteffdict = {}
wsdict = {}
wndict = {}
peakdict = {}
rmsdict = {}
mapdict = {}
cmapdict = {}
file1 = '/spt/user/production/calibration/RCW38-pixelraster/maps/'+obsid_str+'.g3'
f1 = core.G3File(file1)
for frame in f1:
    if frame.type is core.G3FrameType.Map:
        name = frame['Id']
        if 'Wunpol' in frame and name != 'bs':
            twt = np.asarray(frame['Wunpol'].TT)
            twt_orig = twt.copy()
            twt[np.where(twt == 0.)] = 1.
            twt1d = np.reshape(twt,ny*nx)
            medwt = np.median(twt1d[np.where(twt1d > 1.)])
            whmed = np.where(twt1d == medwt)
        if name in names:
            map_weighted = np.asarray(frame['T'])
            reso_arcmin = (frame['T']).res/core.G3Units.arcmin
            if np.max(np.abs(map_weighted)) > 0.:
                map_unw = map_weighted.copy()
                map_unw /= twt
                mapdict[name] = map_unw
                ny,nx = np.shape(map_unw)

                band = np.int(bp[name].band/10.)
                net, w_per_k, whitelevel, watts_sr = \
                    tctools.quick_net_from_rcw38(map_unw,noisedict[name], 
                                                 invert_map=True,band=band,rate=76.3)
                netdict[name] = net
                wsdict[name] = watts_sr
                opteffdict[name] = w_per_k/wpk_perfect[band]
                wndict[name] = whitelevel
                amap_sm = ndimage.gaussian_filter(-map_unw,2./reso_arcmin)
                ycenter, xcenter = np.unravel_index(np.argmax(amap_sm),[ny,nx])
                peakdict[name] = map_unw[ycenter,xcenter]
                pdtemp = tctools.quick_pspec(noisedict[name],rate=76.3)
                psd_dict[name] = pdtemp['psd']
                freqs = pdtemp['freqs']

                map1d = np.reshape(map_unw,ny*nx)
                map1d_whmed = map1d[whmed]
                sdtemp = np.std(map1d_whmed)
                whkeep = np.where(np.abs(map1d_whmed) < 5.*sdtemp)
                sdtemp2 = np.std(map1d_whmed[whkeep])
                rmsdict[name] = sdtemp2
                nsec = medwt/76.3 # assume 76.3 Hz sample rate
                netdict2[name] = sdtemp2*np.sqrt(nsec)/w_per_k

netnames = netdict.keys()
bands = np.asarray([np.int(bp[name].band/10.) for name in netnames])
nets = np.asarray([netdict[name] for name in netnames])
nets2 = np.asarray([netdict2[name] for name in netnames])
opteffs = np.asarray([opteffdict[name] for name in netnames])
wlevels = np.asarray([wndict[name] for name in netnames])
peaks = np.asarray([peakdict[name] for name in netnames])
watts_srs=np.asarray([wsdict[name] for name in netnames])
