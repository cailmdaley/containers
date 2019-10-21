from spt3g import core, dfmux, mapmaker
from spt3g.scratch.tcrawfor import tctools
import numpy as np
import pickle
import scipy
from scipy import ndimage

elnod_dir = pickle.load(open('elnod_dir_03feb17.pkl'))

names = []
for name in elnod_dir.keys():
#    if elnod_dir[name]['elnod_sn'] > 3000:
    if elnod_dir[name]['elnod_sn'] > 50:
        names.append(name)

nbolo = len(names)

file2 = '/spt_data/bolodata/fullrate/RCW38/2653414/nominal_online_cal.g3'
f2 = core.G3File(file2)
bp = f2.next()['NominalBolometerProperties']

#file3 = '/spt_data/bolodata/fullrate/RCW38/2653414/0001.g3'
file3 = '/spt_data/bolodata/fullrate/saturn/2815784/0000.g3'
f3 = core.G3File(file3)
data1 = []
converter = dfmux.unittransforms.ConvertTimestreamUnits(Input='RawTimestreams_I')
i = 0
nfnoise = 10
bdict = {}
for name in names:
    bdict[name] = []
while i < nfnoise:
    frame = f3.next()
    if frame.type is core.G3FrameType.Wiring:
        converter(frame)
    if frame.type is core.G3FrameType.Scan:
        if 'Turnaround' not in frame:
            converter(frame)
            for name in names:
                bdict[name].append(frame['CalTimestreams'][name])
            i += 1

psd_dict = {}
#npts_tot = 0
#for i in np.arange(nfnoise):
#    npts_tot += len(data1[nbolo*i])
#bdata = np.zeros([nbolo,npts_tot])
#npts = 0
#for i in np.arange(nfnoise):
#    tnpts = len(data1[nbolo*i])
#    for j in np.arange(nbolo):
#        bdata[j,npts:npts+tnpts] = data1[nbolo*i+j]
#    npts += tnpts
noisedict = {}
#for i in np.arange(nbolo):
#    noisedict[names[i]] = bdata[i,:]
for name in names:
    noisedict[name] = tctools.list_to_array(bdict[name])
wpk_perfect = {}
kb = 1.38e-23
wpk_perfect[90] = kb*22e9
wpk_perfect[150] = kb*38e9
wpk_perfect[220] = kb*48e9

netdict = {}
opteffdict = {}
wndict = {}
peakdict = {}
mapdir = {}
file1 = '/home/nwhitehorn/saturnmap.g3'
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
                mapdir[frame['Id']] = map_unw
                band = np.int(bp[frame['Id']].band/10.)
                net, w_per_k, whitelevel = \
                    tctools.quick_net_from_rcw38(map_unw,
                                                 noisedict[frame['Id']], 
                                                 invert_map=True,band=band,
                                                 saturn_switch=True)
                netdict[frame['Id']] = net
                opteffdict[frame['Id']] = w_per_k/wpk_perfect[band]
                wndict[frame['Id']] = whitelevel
                amap_sm = ndimage.gaussian_filter(-map_unw,4)
                ycenter, xcenter = np.unravel_index(np.argmax(amap_sm),[360,360])
                peakdict[frame['Id']] = map_unw[ycenter,xcenter]
                pdtemp = tctools.quick_pspec(noisedict[frame['Id']])
                psd_dict[frame['Id']] = pdtemp['psd']
                freqs = pdtemp['freqs']

gnames = netdict.keys()
bands = np.asarray([np.int(bp[name].band/10.) for name in gnames])
nets = np.asarray([netdict[name] for name in gnames])
opteff = np.asarray([opteffdict[name] for name in gnames])
wlevel = np.asarray([wndict[name] for name in gnames])
peaks = np.asarray([peakdict[name] for name in gnames])

