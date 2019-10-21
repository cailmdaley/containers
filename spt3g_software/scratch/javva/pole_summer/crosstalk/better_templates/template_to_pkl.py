#!/usr/bin/env python
import numpy, sys, os
import scipy.ndimage, scipy.interpolate, scipy.optimize
from spt3g import core, mapmaker, calibration, dfmux
import argparse as ap
import copy
import numpy as np 
import pickle

# Usage: 'python make_rcw38_template.py /spt/user/production/calibration/calframe/RCW38-pixelraster/%s.g3 /spt/user/production/calibration/RCW38-pixelraster/maps/%s.g3 /spt/data/bolodata/downsampled/RCW38-pixelraster/%s/0000.g3 -o /spt/user/javva/crosstalk/templates/%s.pkl' %(obsnum,obsnum,obsnum,obsnum)
#
# Makes a template of rcw38 from all of the single bolo maps

P = ap.ArgumentParser(description='Pointing and calibration off of a point source',
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input', action='store', nargs='+', default=[], help='Inpu\
t files') #Input files should be a calframe, by-wafer coadd maps, one scan frame
P.add_argument('-o', '--output', action='store', default='output.pkl',
               help='Output filename')
P.add_argument('-v', '--verbose', action='store_true', default=False)
args = P.parse_args()

# Load source maps from disk
map_extractor = mapmaker.mapmakerutils.ExtractTheMaps()
bp = None
dfm = None
calresponse = None
calresponsesn = None
t_f = None
def extract_bp(fr):
    global dfm, t_f
    if 'DfMuxHousekeeping' in fr and dfm is None:
        dfm = fr['DfMuxHousekeeping']
        print('found dfm')
    if 'DfMuxTransferFunction' in fr and t_f is None:
        t_f = fr['DfMuxTransferFunction']

    global bp, calresponse, calresponsesn,wm
    if 'BolometerProperties' in fr and bp is None:
        bp = fr['BolometerProperties']
    if 'CalibratorResponse' in fr:
        calresponse = fr['CalibratorResponse']
        calresponsesn = fr['CalibratorResponseSN']
    if 'WiringMap' in fr:
        wm = fr['WiringMap']
        print('found the wm')
pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename=args.input)
pipe.Add(extract_bp)
# Cut maps that are all zero
pipe.Add(lambda fr: 'Id' not in fr or 'Wunpol' in fr or (numpy.asarray(fr['T']) != 0).any())
pipe.Add(map_extractor)
pipe.Run()
data = map_extractor.maps
templates_dir = {}
templates_dir['final_temp'] = {}
mapres = list(data.values())[0]['T'].res
print('The map resolution is %s'%mapres)
for d in data.keys():
    templ = copy.copy(np.asarray(data[d]['T']/np.max(np.abs(data[d]['T']))))
    templ[:templ.shape[0]//2 - 30,:] = 0
    templ[templ.shape[0]//2 + 30:,:] = 0
    templ[:,:templ.shape[1]//2 - 30] = 0
    templ[:,templ.shape[1]//2 + 30:] = 0
    templates_dir['final_temp'][d] = templ

#build up info about all of the bolometers best fit values, as well as previously saved values
off_params = {}
for bolo in bp.keys():
    off_params[bolo] = {}
    off_params[bolo]['pixel_id'] = bp[bolo].physical_name
    off_params[bolo]['x_offset'] = (bp[bolo].x_offset)/mapres
    off_params[bolo]['y_offset'] = (bp[bolo].y_offset)/mapres
n = 0

#Save all of the bolometers carrier frequencies
carrier_freqs = {}
for bolo in bp.keys():
    try:
        carrier_freqs[bolo] = dfmux.HousekeepingForBolo(dfm,wm, bolo).carrier_frequency 
    except:
        n += 1
print('could not find wiring for %s out of %s bolos'%(n, len(list(bp.keys()))))
templates_dir['wiring_info'] = {}

#Save the information about pixel mapping
for bolo in bp.keys():
    try:
        templates_dir['wiring_info'][bolo]= dfmux.HardwareMapTools.PathStringForBolo(wm, bolo)
    except:
        templates_dir['wiring_info'][bolo]= 'bad'

templates_dir['params'] = {}
templates_dir['params']['fit_params'] = 'none for new templtes'
templates_dir['params']['official_params'] = off_params
templates_dir['carrier_freqs'] = carrier_freqs

with open(args.output, 'wb') as handle:
    pickle.dump(templates_dir, handle, protocol=pickle.HIGHEST_PROTOCOL)

