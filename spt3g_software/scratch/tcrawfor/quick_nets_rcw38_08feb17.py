#from spt3g import core, dfmux
#from spt3g.scratch.tcrawfor import tctools
#import numpy as np
#import pickle
#
#names = pickle.load(open('elnod_gnames_08feb17.pkl'))
#
#file3 = '/spt_data/bolodata/fullrate/RCW38/3265922/0000.g3'
#f3 = core.G3File(file3)
#bdict = {}
#for name in names:
#    bdict[name] = []
#
#converter = dfmux.unittransforms.ConvertTimestreamUnits(Input='RawTimestreams_I')
#i = 0
#nfnoise = 10
#while i < nfnoise:
#    frame = f3.next()
#    if frame.type is core.G3FrameType.Wiring:
#        converter(frame)
#    if frame.type is core.G3FrameType.Scan:
#        if 'Turnaround' not in frame:
#            converter(frame)
#            for name in names:
#                try:
#                    bdict[name].append(frame['CalTimestreams'][name])
#                except:
#                    pass
#            i += 1

noisedict = {}
for name in names:
    if len(bdict[name]) > 1:
        noisedict[name] = tctools.list_to_array(bdict[name])
wpk_perfect = {}
wpk_perfect[90] = kb*22e9
wpk_perfect[150] = kb*38e9
wpk_perfect[220] = kb*48e9

netdict = {}
opteffdict = {}
wndict = {}
peakdict = {}
file1 = '/poleanalysis/cal_results/RCW38/3265922.bolomaps.g3.gz'
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

netnames = netdict.keys()

bands = np.asarray([np.int(bp[name].band/10.) for name in netnames])
nets = np.asarray([netdict[name] for name in netnames])
opteff = np.asarray([opteffdict[name] for name in netnames])
wlevel = np.asarray([wndict[name] for name in netnames])


