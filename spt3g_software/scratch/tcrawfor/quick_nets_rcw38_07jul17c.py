from spt3g import core, dfmux, std_processing
#from spt3g.scratch.tcrawfor import tctools
from spt3g.util import tctools
import numpy as np
import pickle
from scipy import ndimage

#obsid_str = '16508393'
obsid_str = '2653414'

file2 = '/spt/data/bolodata/downsampled/RCW38/'+obsid_str+'/nominal_online_cal.g3'
#file2 = '/spt/data/bolodata/downsampled/RCW38-pixelraster/'+obsid_str+'/nominal_online_cal.g3'
f2 = core.G3File(file2)
bp = f2.next()['NominalBolometerProperties']
names = bp.keys()
nbolo = len(names)

# get noise data & convert to watts
file3 = '/spt/data/bolodata/downsampled/RCW38/'+obsid_str+'/0000.g3'
#file3 = '/spt/data/bolodata/downsampled/RCW38-pixelraster/'+obsid_str+'/0000.g3'
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
    noisedict[name] = noisedict[name][npts_running/4:3*npts_running/4]

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
opteffdict = {}
wndict = {}
wndict2 = {}
peakdict = {}
rmsdict = {}
mapdict = {}
cmapdict = {}
file1 = '/spt/user/production/calibration/RCW38/maps/'+obsid_str+'.g3'
#file1 = '/spt/user/production/calibration/RCW38-pixelraster/maps/'+obsid_str+'.g3'
f1 = core.G3File(file1)
for frame in f1:
    if frame.type is core.G3FrameType.Map:
        name = frame['Id']
        if 'Wunpol' in frame and name != 'bs':
            twt = np.asarray(frame['Wunpol'].TT)
            twt_orig = twt.copy()
            twt[np.where(twt == 0.)] = 1.
        if name in names:
            map_weighted = np.asarray(frame['T'])
            reso_arcmin = (frame['T']).res/core.G3Units.arcmin
            if np.max(np.abs(map_weighted)) > 0.:
                map_unw = map_weighted.copy()
                map_unw /= twt
                mapdict[name] = map_unw
# do some outlier cutting. if I don't do this, hot pixels get
# identified as rcw38 center. first cut everything to 10 sigma, then
# smooth and look for source, then mask source and cut to 5
# sigma. this relies on source being many >~10-sigma pixels together.
                ny,nx = np.shape(map_unw)
                map1d = np.reshape(map_unw.copy(),ny*nx)
                twt1d = np.reshape(twt_orig,ny*nx)
                whn0 = np.where(twt1d > 0.)
                sdtemp = np.std(map1d[whn0])
                out_thresh = 1000.
                whout = np.where(np.abs(map1d)/sdtemp > out_thresh)
                while len(whout[0]) > 0:
                    map1d[whout[0]] = 0.
                    sdtemp = np.std(map1d[whn0])
                    whout = np.where(np.abs(map1d)/sdtemp > out_thresh)
                map_unw_temp = np.reshape(map1d,[ny,nx])
                amap_sm_temp = ndimage.gaussian_filter(-map_unw_temp,2./reso_arcmin)
                yctemp, xctemp = np.unravel_index(np.argmax(amap_sm_temp),[ny,nx])
                npmask = np.int(np.round(3./reso_arcmin))
                try:
                    mask = np.ones([ny,nx])
                    mask[yctemp-npmask:yctemp+npmask,xctemp-npmask:xctemp+npmask] = 0.
                except:
                    mask = np.ones([ny,nx])
                map1d_1 = map1d.copy()
                whn0_1 = (np.asarray(whn0[0])).copy()
                map1d = np.reshape(map_unw.copy(),ny*nx)
                mask1d = np.reshape(mask,ny*nx)
                whn0 = np.where(twt1d*mask1d > 0.)
                sdtemp = np.std(map1d[whn0])
                out_thresh = 500.
                whout = np.where(np.abs(map1d*mask1d)/sdtemp > out_thresh)
                while len(whout[0]) > 0:
                    map1d[whout] = 0.
                    sdtemp = np.std(map1d[whn0])
                    whout = np.where(np.abs(map1d*mask1d)/sdtemp > out_thresh)
                map_unw = np.reshape(map1d,[ny,nx])

                band = np.int(bp[name].band/10.)
                net, w_per_k, whitelevel = \
                    tctools.quick_net_from_rcw38(map_unw,noisedict[name], 
                                                 invert_map=True,band=band, 
                                                 rate=76.3,frange=[21.,24.])
                netdict[name] = net
                opteffdict[name] = w_per_k/wpk_perfect[band]
                wndict[name] = whitelevel
                wndict2[name] = np.std(noisedict[name])/np.sqrt(76.3/2.)
                amap_sm = ndimage.gaussian_filter(-map_unw,2./reso_arcmin)
                ycenter, xcenter = np.unravel_index(np.argmax(amap_sm),[ny,nx])
                peakdict[name] = map_unw[ycenter,xcenter]
                rmsdict[name] = sdtemp
                pdtemp = tctools.quick_pspec(noisedict[name],rate=76.3)
                psd_dict[name] = pdtemp['psd']
                freqs = pdtemp['freqs']
                cmapdict[name] = map_unw
#                if name == 'W148/2017.W148.2.78.1.4070':
#                print notavariable

netnames = netdict.keys()
bands = np.asarray([np.int(bp[name].band/10.) for name in netnames])
nets = np.asarray([netdict[name] for name in netnames])
opteffs = np.asarray([opteffdict[name] for name in netnames])
wlevels = np.asarray([wndict[name] for name in netnames])
wlevels2 = np.asarray([wndict2[name] for name in netnames])
peaks = np.asarray([peakdict[name] for name in netnames])
rmss = np.asarray([rmsdict[name] for name in netnames])


