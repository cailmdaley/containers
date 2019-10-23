from spt3g import core, dfmux
from spt3g.scratch.tcrawfor import tctools
import numpy as np
import pickle

file2 = '/spt_data/bolodata/fullrate/RCW38/2653414/nominal_online_cal.g3'
f2 = core.G3File(file2)
bp = f2.next()['NominalBolometerProperties']

whitenoise_w = pickle.load(open('whitenoise_w_02feb17.pkl'))
names = pickle.load(open('bolo_names_02feb17.pkl'))
nbolo = len(names)
wndict = {}
for i in np.arange(nbolo):
    wndict[names[i]] = whitenoise_w[i]

file3 = '/spt_data/bolodata/fullrate/RCW38/2653414/0001.g3'
f3 = core.G3File(file3)
data1 = []
converter = dfmux.unittransforms.ConvertTimestreamUnits(Input='RawTimestreams_I')
i = 0
nfnoise = 10
while i < nfnoise:
    frame = f3.next()
    if frame.type is core.G3FrameType.Wiring:
        converter(frame)
    if frame.type is core.G3FrameType.Scan:
        if 'Turnaround' not in frame:
            converter(frame)
            for name in names:
                data1.append(frame['CalTimestreams'][name])
            i += 1

psd_dict = {}
npts_tot = 0
for i in np.arange(nfnoise):
    npts_tot += len(data1[nbolo*i])
bdata = np.zeros([nbolo,npts_tot])
npts = 0
for i in np.arange(nfnoise):
    tnpts = len(data1[nbolo*i])
    for j in np.arange(nbolo):
        bdata[j,npts:npts+tnpts] = data1[nbolo*i+j]
    npts += tnpts
noisedict = {}
for i in np.arange(nbolo):
    noisedict[names[i]] = bdata[i,:]
wpk_perfect = {}
wpk_perfect[90] = kb*22e9
wpk_perfect[150] = kb*38e9
wpk_perfect[220] = kb*48e9


netdict = {}
opteffdict = {}
wndict = {}
peakdict = {}
file1 = '/home/nwhitehorn/rcw38-secondlight-2.g3'
f1 = core.G3File(file1)
for frame in f1:
    if frame.type is core.G3FrameType.Map:
        if 'Wunpol' in frame and frame['Id'] != 'bs':
            twt = np.asarray(frame['Wunpol'].TT)
            twt[np.where(twt == 0.)] = 1.
        if frame['Id'] in names:
            map_weighted = np.asarray(frame['T'])
            map_unw = map_weighted.copy()
            map_unw /= twt
            band = np.int(bp[frame['Id']].band/10.)
            net, w_per_k, whitelevel = \
                tctools.quick_net_from_rcw38(map_unw,noisedict[frame['Id']], 
                                             invert_map=True,band=band)
            netdict[frame['Id']] = net
            opteffdict[frame['Id']] = w_per_k/wpk_perfect[band]
            wndict[frame['Id']] = whitelevel
            amap_sm = ndimage.gaussian_filter(-map_unw,4)
            ycenter, xcenter = np.unravel_index(np.argmax(amap_sm),[360,360])
            peakdict[frame['Id']] = map_unw[ycenter,xcenter]
            pdtemp = tctools.quick_pspec(noisedict[frame['Id']])
            psd_dict[frame['Id']] = pdtemp['psd']
            freqs = pdtemp['freqs']

bands = np.asarray([np.int(bp[name].band/10.) for name in names])
nets = np.asarray([netdict[name] for name in names])
opteff = np.asarray([opteffdict[name] for name in names])
wlevel = np.asarray([wndict[name] for name in names])


