#!/usr/bin/env python
import numpy, sys, os
import scipy.ndimage, scipy.interpolate, scipy.optimize
from spt3g import core, mapmaker, calibration, util
import argparse as ap

# Usage: fit_fluxandpointing.py <files.g3> -o output.g3
#
# Computes best-fit relative detector pointing and flux calibration
# for all detectors in the focal plane. Flux results are normalized to
# a 4 arcminute x 4 arcminute box centered on the brightest point
# in the input maps.

P = ap.ArgumentParser(description='Pointing and calibration off of a point source',
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('-ip','--input', default=['/spt_data/bolodata/downsampled/RCW38-pixelraster/57241701/offline_calibration.g3', '/poleanalysis/sptdaq/calresult/calibration/RCW38-pixelraster/maps/56951575.g3'], help = 'Input Files')#'/spt_data/bolodata/downsampled/RCW38-pixelraster/57241701/0000.g3','/spt_data/bolodata/downsampled/RCW38-pixelraster/57241701/0001.g3','/spt_data/bolodata/downsampled/RCW38-pixelraster/57241701/0002.g3','/spt_data/bolodata/downsampled/RCW38-pixelraster/57241701/0003.g3','/spt_data/bolodata/downsampled/RCW38-pixelraster/57241701/0004.g3','/spt_data/bolodata/downsampled/RCW38-pixelraster/57241701/0005.g3','/spt_data/bolodata/downsampled/RCW38-pixelraster/57241701/0006.g3'],
               #help='Input files')
P.add_argument('-s', '--source', action='store', default='RCW38',
               help='Source name for output data')
P.add_argument('-o', '--output', action='store', default='output.g3',
               help='Output filename')
P.add_argument('-v', '--verbose', action='store_true', default=False,
               help='Print internal likelihoods')

args = P.parse_args()

# Load source maps from disk
import glob
im = sorted(glob.glob('/poleanalysis/sptdaq/calresult/calibration/RCW38-pixelraster/maps/*'))[-10:]
import numpy as np
input_files = np.append('/spt_data/bolodata/downsampled/RCW38-pixelraster/57241701/offline_calibration.g3', im)
map_extractor = mapmaker.mapmakerutils.ExtractTheMaps()
bp = None
calresponse = None
calresponsesn = None

pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename=input_files[2])
pipe.Add(map_extractor)
pipe.Run()
data = map_extractor.maps


def extract_bp(fr):
    global bp, calresponse, calresponsesn
    if 'Id' in fr:
        if False in np.isfinite(np.asarray(fr['T'])):
            #print('baddie')
            return False
        if fr['Id'] != sample_id:
            #print(fr)
            return False
        if fr['Id'] == sample_id:
            print(fr)
    if 'BolometerProperties' in fr and bp is None:
        bp = fr['BolometerProperties']
    if 'CalibratorResponse' in fr:
        calresponse = fr['CalibratorResponse']
        calresponsesn = fr['CalibratorResponseSN']
for sample_id in sorted(data.keys()):
    print('DOING BOLOMETER %s'%sample_id)
    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=input_files)
    pipe.Add(extract_bp)
    # Cut maps that are all zero
    pipe.Add(lambda fr: 'Id' not in fr or 'Wunpol' in fr or (numpy.asarray(fr['T']) != 0).any())
    pipe.Add(util.framecombiner.MapFrameCombiner, fr_id = sample_id)
    pipe.Add(lambda fr: fr.type == core.G3FrameType.Map) # Drop TOD         
    pipe.Add(core.G3Writer, filename = '/poleanalysis/javva/rcw38_singlebolo_coadds/%s_coadd.g3'%sample_id)
    pipe.Run()
