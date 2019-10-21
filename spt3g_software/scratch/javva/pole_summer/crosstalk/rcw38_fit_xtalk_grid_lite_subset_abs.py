#!/usr/bin/env python
import numpy, sys, os
import scipy.ndimage, scipy.interpolate, scipy.optimize
from spt3g import core, mapmaker, calibration
import scipy
import argparse as ap
import numpy as np 
import pickle
# Usage: fit_fluxandpointing.py <files.g3> -o output.g3
#
# Computes best-fit relative detector pointing and flux calibration
# for all detectors in the focal plane. Flux results are normalized to
# a 4 arcminute x 4 arcminute box centered on the brightest point
# in the input maps.

P = ap.ArgumentParser(description='Pointing and calibration off of a point source',
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
#P.add_argument('-ip','--input', default=['/spt/data/bolodata/downsampled/RCW38-pixelraster/64584651/offline_calibration.g3', '/spt/user/production/calibration/RCW38-pixelraster/maps/64584651.g3'], help = 'Input Files')
P.add_argument('input_files', action='store', nargs='+', default=[], help='Inpu\
t files')
P.add_argument('-mb','--mb', default ='None', help = 'number corresponding to observation you want to do the fit for')
P.add_argument('-o', '--output', action='store', default='output.pkl',
               help='Output filename')
P.add_argument('-lstsq', '--lstsq', action='store', default=False,
               help='Output filename')
P.add_argument('-v', '--verbose', action='store_true', default=False)
args = P.parse_args()
print(args.input_files)
templates = pickle.load(open(args.input_files[0],'rb'))['final_temp']
params = pickle.load(open(args.input_files[0],'rb'))['params']
rawmap = pickle.load(open(args.input_files[1],'rb'))['data']
wiring_info = pickle.load(open(args.input_files[0],'rb'))['wiring_info']
carrier_freqs = pickle.load(open(args.input_files[0],'rb'))['carrier_freqs']
w = pickle.load(open(args.input_files[1],'rb'))['w']
pix_info = params['official_params']
print(params.keys())
print(templates.keys())
#do least squares fit 

main_bolo = args.mb
other_places = {}
for group in ['90.0','150.0','220.0']:
    other_places[group] = {}
    print(main_bolo)
    unwmap = rawmap/w
    unwmap[numpy.logical_not(numpy.isfinite(unwmap))] = 0
    invnoise = w**0.5
    invnoise /= numpy.std(unwmap[unwmap != 0]*invnoise[unwmap != 0])
    params_amps = {**params['fit_params']['90.0'],**params['fit_params']['150.0'], **params['fit_params']['220.0']}
    print(params_amps)
    template_spl = scipy.interpolate.RectBivariateSpline(numpy.arange(templates[group].shape[0]),numpy.arange(templates[group].shape[1]),templates[group])
    print('doing lstsq on %s bolos'%len(params_amps.keys()))
    bif = 0
    len_empty = [bolo for bolo in (sorted(params_amps.keys())) if '/'.join(wiring_info[main_bolo].split('/')[:4]) == '/'.join(wiring_info[bolo].split('/')[:4])] 

    #only leave one bolo from each pixel -- the one closest in frequency
    pix_group = {}
    for bolo in len_empty:
        pixel = pix_info[bolo]['pixel_id'].split('.')[0]
        if pixel not in pix_group.keys():
            pix_group[pixel] = [b for b in params_amps.keys() if pix_info[b]['pixel_id'].split('.')[0] == pixel]
    main_frequency = carrier_freqs[main_bolo]
    use_bolos = []
    
    for pixel in pix_group.keys():
        delta_prime = 100 
        for bolo in pix_group[pixel]:
            delta = np.abs(main_frequency - carrier_freqs[bolo])
            print(main_frequency, carrier_freqs[bolo])
            if delta < delta_prime:
                delta_prime = delta
                bolo_prime = bolo
        use_bolos = np.append(use_bolos, bolo_prime)
    len_empty = list(use_bolos)
    a1 = np.zeros([len(len_empty), len(np.ravel(template_spl(np.arange(templates[group].shape[0]) - params_amps[main_bolo][1], numpy.arange(templates[group].shape[1]) - params_amps[main_bolo][2])))])
    for i, bolo in enumerate(sorted(len_empty)):
        if '/'.join(wiring_info[main_bolo].split('/')[:4]) == '/'.join(wiring_info[bolo].split('/')[:4]):
            a1[bif] = np.ravel(template_spl(np.arange(templates[group].shape[0]) - params_amps[bolo][1], numpy.arange(templates[group].shape[1]) - params_amps[bolo][2])*invnoise)
            bif +=1
    print('Starting fit on %s bolos'%bif)
#    co = scipy.optimize.nnls(a1.T,np.abs(np.ravel(unwmap*invnoise)))
    co = np.linalg.lstsq((a1/1e15).T,np.ravel(unwmap*invnoise))
    bif = 0
    a1_noise = np.zeros([len(len_empty), len(np.ravel(template_spl(np.arange(templates[group].shape[0]) - params_amps[main_bolo][1], numpy.arange(templates[group].shape[1]) - params_amps[main_bolo][2])))])
    for i, bolo in enumerate(sorted(len_empty)):
        if '/'.join(wiring_info[main_bolo].split('/')[:4]) == '/'.join(wiring_info[bolo].split('/')[:4]):
            a1_noise[bif] = np.ravel(template_spl(np.arange(templates[group].shape[0]) + params_amps[bolo][1], numpy.arange(templates[group].shape[1]) + params_amps[bolo][2])*invnoise)
            bif +=1
    co_noise = np.linalg.lstsq((a1_noise/1e15).T,np.ravel(unwmap*invnoise))

    other_places[group]['on_source'] = co
    other_places[group]['off_source'] = co_noise
    other_places[group]['other_bolos_order'] = sorted(len_empty)
    other_places[group]['location'] = params_amps
with open(args.output, 'wb') as handle:
    pickle.dump(other_places, handle, protocol=pickle.HIGHEST_PROTOCOL)

