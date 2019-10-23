#!/usr/bin/env python

from spt3g import core, dfmux, std_processing, todfilter
from spt3g.util import tctools
import numpy as np
import pickle
from scipy import ndimage
import sys

# Usage python quick_nets_rcw38_tc.py obsid output_file [fminstr fmaxstr data_dir calresults_dir]

args = sys.argv
nargs = len(args)
if nargs < 2:
    print('Usage: python quick_nets_rcw38_tc.py obsid output_file [fminstr fmaxstr data_dir calresults_dir]')
    print('\n')
obsid_str = sys.argv[1]
outfile = sys.argv[2]
fminstr = '10'
if nargs > 3:
    fminstr = sys.argv[3]
fmaxstr = '20'
if nargs > 4:
    fmaxstr = sys.argv[4]
fmin = np.float(fminstr)
fmax = np.float(fmaxstr)
#data_dir = '/poleanalysis/tcrawfor/'
data_dir = '/spt_data'
if nargs > 5:
    data_dir = sys.argv[5]
#calresults_dir = '/poleanalysis/tcrawfor/'
calresults_dir = '/poleanalysis/sptdaq'
if nargs > 6:
    calresults_dir = sys.argv[6]

file2 = data_dir+'/bolodata/fullrate/RCW38-pixelraster/'+obsid_str+'/nominal_online_cal.g3'
f2 = core.G3File(file2)
bp = f2.next()['NominalBolometerProperties']
names = bp.keys()
nbolo = len(names)

# get noise data & convert to watts
file3 = data_dir+'/bolodata/fullrate/RCW38-pixelraster/'+obsid_str+'/0000.g3'
f3 = core.G3File(file3)
converter = dfmux.unittransforms.ConvertTimestreamUnits(Input='RawTimestreams_I')
i = 0
nfnoise = 20
wndict = {}
while i < nfnoise:
    frame = f3.next()
    if frame.type is core.G3FrameType.Wiring:
        converter(frame)
    if frame.type is core.G3FrameType.Scan:
        if 'Turnaround' not in frame:
            if i == 0:
                for name in frame['RawTimestreams_I'].keys(): 
                    wndict[name] = 0.
            converter(frame)
            filtered_data = todfilter.polyutils.poly_filter_g3_timestream_map(frame['CalTimestreams'], 4)
            psd_dict_temp, freqs_temp =  todfilter.dftutils.get_psd_of_ts_map(filtered_data)
            fhz = freqs_temp/core.G3Units.Hz
            inds = np.where(np.logical_and(fhz >= fmin,fhz <= fmax))[0]
            for name in filtered_data.keys():
                wndict[name] += np.median(np.asarray(psd_dict_temp[name])[inds])/np.float(nfnoise)*core.G3Units.Hz
            i += 1
for name in wndict.keys():
    wndict[name] = np.sqrt(wndict[name])


# W/K for a perfect single-pol bolo with assumed bandwidths
kb = 1.38e-23
wpk_perfect = {}
#wpk_perfect[90] = kb*22e9
#wpk_perfect[150] = kb*38e9
#wpk_perfect[220] = kb*48e9
# note, the following are designed bandwidths (from microstrip
# simulations). comparing to these to get optical efficiency means
# that differences in bandwidth get lumped into optical efficiency.
wpk_perfect[90] = kb*28.2e9 # new from Brad
wpk_perfect[150] = kb*42.6e9
wpk_perfect[220] = kb*51.9e9

# integrate maps & do NET calc
netdict = {}
netdict2 = {}
opteffdict = {}
wsdict = {}
peakdict = {}
rmsdict = {}
mapdict = {}
cmapdict = {}
file1 = calresults_dir+'/calresult/calibration/RCW38-pixelraster/maps/'+obsid_str+'.g3'
f1 = core.G3File(file1)
for frame in f1:
    if frame.type is core.G3FrameType.Map:
        name = frame['Id']
        if 'Wunpol' in frame and name != 'bs':
            twt = np.asarray(frame['Wunpol'].TT)
            twt_orig = twt.copy()
            ny,nx = np.shape(twt_orig)
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
                    tctools.quick_net_from_rcw38(map_unw,wndict[name], 
                                                 invert_map=True,band=band,rate=152.6,frange=[fmin,fmax],data_is_wn=True)
                netdict[name] = net
                wsdict[name] = watts_sr
                opteffdict[name] = w_per_k/wpk_perfect[band]
                wndict[name] = whitelevel
                amap_sm = ndimage.gaussian_filter(-map_unw,2./reso_arcmin)
                ycenter, xcenter = np.unravel_index(np.argmax(amap_sm),[ny,nx])
                peakdict[name] = map_unw[ycenter,xcenter]

netnames = netdict.keys()
bands = np.asarray([np.int(bp[name].band/10.) for name in netnames])
nets = np.asarray([netdict[name] for name in netnames])
opteffs = np.asarray([opteffdict[name] for name in netnames])
wlevels = np.asarray([wndict[name] for name in netnames])
peaks = np.asarray([peakdict[name] for name in netnames])
watts_srs = np.asarray([wsdict[name] for name in netnames])

cfs = np.asarray([1.23,1.73,3.07])
nets_cmb = nets.copy()
wh90 = np.where(bands == 90)
wh150 = np.where(bands == 150)
wh220 = np.where(bands == 220)
nets_cmb[wh90] *= cfs[0]
nets_cmb[wh150] *= cfs[1]
nets_cmb[wh220] *= cfs[2]
netdict_rj = netdict.copy()
netdict = {}
banddict = {}
for i in np.arange(len(netnames)):
    netdict[netnames[i]] = nets_cmb[i]
    banddict[netnames[i]] = bands[i]

outdict = {}
outdict['NET'] = netdict
outdict['NET_RJ'] = netdict_rj
outdict['opteff'] = opteffdict
outdict['whitelevel'] = wndict
outdict['band'] = banddict
outdict['response_in_watts_sr'] = wsdict

pickle.dump(outdict,open(outfile,'w'))
