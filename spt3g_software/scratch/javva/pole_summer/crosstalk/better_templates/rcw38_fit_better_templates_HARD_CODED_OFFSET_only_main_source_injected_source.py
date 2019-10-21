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

#grab some information from the input files
templates = pickle.load(open(args.input_files[0],'rb'))['final_temp']
params = pickle.load(open(args.input_files[0],'rb'))['params']
rawmap = pickle.load(open(args.input_files[1],'rb'))['data']
wiring_info = pickle.load(open(args.input_files[0],'rb'))['wiring_info']
carrier_freqs = pickle.load(open(args.input_files[0],'rb'))['carrier_freqs']
w = pickle.load(open(args.input_files[1],'rb'))['w']
pix_info = params['official_params']
print(params.keys())
print(templates.keys())

main_bolo = args.mb
mb = main_bolo
print('Doing the fit for bolometer', main_bolo)
unwmap = rawmap/w
unwmap[numpy.logical_not(numpy.isfinite(unwmap))] = 0
invnoise = w**0.5
invnoise /= numpy.std(unwmap[unwmap != 0]*invnoise[unwmap != 0])
params_amps = params['official_params']
mb_temp = 'rcw38-'+pix_info[mb]['pixel_id'].split('.')[1]+'GHzW'+pix_info[mb]['pixel_id'].split('w')[1].split('_')[0]


print('doing lstsq on %s bolos'%len(params_amps.keys()))
bif = 0
#grab the bolos for the fit that share a comb with the main bolo
len_empty = [bolo for bolo in (sorted(params_amps.keys())) if '/'.join(wiring_info[main_bolo].split('/')[:4]) == '/'.join(wiring_info[bolo].split('/')[:4])] 
#group those bolometers by pixel
pix_group = {}
for bolo in sorted(len_empty):
    pixel = pix_info[bolo]['pixel_id'].split('.')[0]
    if pixel not in pix_group.keys():
        if 'resistor' in pixel:
            continue
        pix_group[pixel] = [b for b in params_amps.keys() if pix_info[b]['pixel_id'].split('.')[0] == pixel and '/'.join(wiring_info[main_bolo].split('/')[:4]) == '/'.join(wiring_info[b].split('/')[:4]) and np.isfinite(pix_info[b]['x_offset'])]
        if pix_group[pixel] == []:
            print('found an empty set')
            del pix_group[pixel]
main_frequency = carrier_freqs[main_bolo]
use_bolos = []

for pixel in pix_group.keys():
    delta_prime = 100 #just a big number, no significance
    print('starting', pixel)
    print(pix_group[pixel])
    skp = 0 
    for bolo in pix_group[pixel]:
        if bolo not in carrier_freqs.keys():
            print('no dfmux info on %s'%bolo)
            if len(list(pix_group[pixel])) == 1:
                skp = 1
            continue
        if 'resistor' in pix_info[bolo]['pixel_id']:
            if len([k for k in pix_group[pixel] if 'resistor' in pix_info[k]['pixel_id']]) == len(pix_group[pixel]):
                skp = 1
            print(pix_info[bolo]['pixel_id'], 'resis')
            continue
        delta = np.abs(main_frequency - carrier_freqs[bolo])
        print(bolo,main_frequency, carrier_freqs[bolo], delta)
        if delta < delta_prime:
            delta_prime = delta
            bolo_prime = bolo
            print('used bolo', bolo)
    if skp > 0:
        continue
    print('actually used', bolo_prime)
    use_bolos = np.append(use_bolos, bolo_prime)
print('use bolos,', use_bolos)
template_spl = scipy.interpolate.RectBivariateSpline(numpy.arange(templates[mb_temp].shape[0]),numpy.arange(templates[mb_temp].shape[1]),templates[mb_temp]) 
len_empty = list(use_bolos)
a1 = np.zeros([len(len_empty), len(np.ravel(template_spl(np.arange(templates[mb_temp].shape[0]) - pix_info[mb]['y_offset'], numpy.arange(templates[mb_temp].shape[1]) - pix_info[mb]['x_offset'])))])
print('This is the order', sorted((len_empty)))
for i, bolo in enumerate(sorted(len_empty)):
    b_temp = 'rcw38-'+pix_info[bolo]['pixel_id'].split('.')[1]+'GHzW'+pix_info[bolo]['pixel_id'].split('w')[1].split('_')[0]
    template_spl = scipy.interpolate.RectBivariateSpline(numpy.arange(templates[b_temp].shape[0]),numpy.arange(templates[b_temp].shape[1]),templates[b_temp])
    a1[bif] = np.ravel(template_spl(np.arange(templates[b_temp].shape[0]) - pix_info[bolo]['y_offset'], numpy.arange(templates[b_temp].shape[1]) - pix_info[bolo]['x_offset'])*invnoise)
    bif +=1
print('Starting fit on %s bolos'%bif)
for i in a1:
    print(i)
co = np.linalg.lstsq((a1/1e10).T,np.ravel(unwmap*invnoise))

#Do the fit to the xtalk region, but only with the main bolometer
bif = 0
template_spl = scipy.interpolate.RectBivariateSpline(numpy.arange(templates[mb_temp].shape[0]),numpy.arange(templates[mb_temp].shape[1]),templates[mb_temp])
len_empty = list(use_bolos)
a1 = np.zeros([1, len(np.ravel(template_spl(np.arange(templates[mb_temp].shape[0]) - pix_info[mb]['y_offset'], numpy.arange(templates[mb_temp].shape[1]) - pix_info[mb]['x_offset'])))])
print('This is the order', sorted((len_empty)))
for i, bolo in enumerate(sorted(len_empty)):
    if bolo != mb:
        continue
    b_temp = 'rcw38-'+pix_info[bolo]['pixel_id'].split('.')[1]+'GHzW'+pix_info[bolo]['pixel_id'].split('w')[1].split('_')[0]
    template_spl = scipy.interpolate.RectBivariateSpline(numpy.arange(templates[b_temp].shape[0]),numpy.arange(templates[b_temp].shape[1]),templates[b_temp])
    a1[bif] = np.ravel(template_spl(np.arange(templates[b_temp].shape[0]) - pix_info[bolo]['y_offset'], numpy.arange(templates[b_temp].shape[1]) - pix_info[bolo]['x_offset'])*invnoise)
    bif +=1
print('Starting fit on %s bolos'%bif)
co_main_bolo = np.linalg.lstsq((a1/1e10).T,np.ravel(unwmap*invnoise))


#Starting the fit on the noise region (all components)

bif = 0
a1_noise = np.zeros([len(len_empty), len(np.ravel(template_spl(np.arange(templates[mb_temp].shape[0]) - pix_info[mb]['y_offset'], numpy.arange(templates[mb_temp].shape[1]) - pix_info[mb]['x_offset'])))])
for i, bolo in enumerate(sorted(len_empty)):
    b_temp = 'rcw38-'+pix_info[bolo]['pixel_id'].split('.')[1]+'GHzW'+pix_info[bolo]['pixel_id'].split('w')[1].split('_')[0]
    template_spl = scipy.interpolate.RectBivariateSpline(numpy.arange(templates[b_temp].shape[0]),numpy.arange(templates[b_temp].shape[1]),templates[b_temp])
    a1_noise[bif] = np.ravel(template_spl(np.arange(templates[b_temp].shape[0]) + pix_info[bolo]['y_offset'], numpy.arange(templates[b_temp].shape[1]) + pix_info[bolo]['x_offset'])*invnoise)
    bif +=1
print('Starting fit on %s bolos'%bif)
co_noise = np.linalg.lstsq((a1_noise/1e10).T,np.ravel(unwmap*invnoise))

#Starting the fit on the noise region (only main bolometer)
bif = 0
template_spl = scipy.interpolate.RectBivariateSpline(numpy.arange(templates[mb_temp].shape[0]),numpy.arange(templates[mb_temp].shape[1]),templates[mb_temp])
len_empty = list(use_bolos)
a1 = np.zeros([1, len(np.ravel(template_spl(np.arange(templates[mb_temp].shape[0]) - pix_info[mb]['y_offset'], numpy.arange(templates[mb_temp].shape[1]) - pix_info[mb]['x_offset'])))])
print('This is the order', sorted((len_empty)))
for i, bolo in enumerate(sorted(len_empty)):
    if bolo != mb:
        continue
    b_temp = 'rcw38-'+pix_info[bolo]['pixel_id'].split('.')[1]+'GHzW'+pix_info[bolo]['pixel_id'].split('w')[1].split('_')[0]
    template_spl = scipy.interpolate.RectBivariateSpline(numpy.arange(templates[b_temp].shape[0]),numpy.arange(templates[b_temp].shape[1]),templates[b_temp])
    a1[bif] = np.ravel(template_spl(np.arange(templates[b_temp].shape[0]) + pix_info[bolo]['y_offset'], numpy.arange(templates[b_temp].shape[1]) + pix_info[bolo]['x_offset'])*invnoise)
    bif +=1
print('Starting fit on %s bolos'%bif)
co_noise_main_bolo = np.linalg.lstsq((a1/1e10).T,np.ravel(unwmap*invnoise))

#make the fake source
a1_neg = np.zeros([360,360])
bif = 0
for i, bolo in enumerate(sorted(len_empty)):
    if bolo != mb:
        continue
    b_temp = 'rcw38-'+pix_info[bolo]['pixel_id'].split('.')[1]+'GHzW'+pix_info[bolo]['pixel_id'].split('w')[1].split('_')[0]
    template_spl = scipy.interpolate.RectBivariateSpline(numpy.arange(templates[b_temp].shape[0]),numpy.arange(templates[b_temp].shape[1]),templates[b_temp])
    a1_neg += co[0][i]*(template_spl(np.arange(templates[b_temp].shape[0]) + pix_info[bolo]['y_offset'], numpy.arange(templates[b_temp].shape[1]) + pix_info[bolo]['x_offset'])*invnoise)
    bif +=1


#starting fit to the noise region with injected source
bif = 0
template_spl = scipy.interpolate.RectBivariateSpline(numpy.arange(templates[mb_temp].shape[0]),numpy.arange(templates[mb_temp].shape[1]),templates[mb_temp])
len_empty = list(use_bolos)
a1 = np.zeros([len(len_empty), len(np.ravel(template_spl(np.arange(templates[mb_temp].shape[0]) - pix_info[mb]['y_offset'], numpy.arange(templates[mb_temp].shape[1]) - pix_info[mb]['x_offset'])))])
print('This is the order', sorted((len_empty)))
for i, bolo in enumerate(sorted(len_empty)):
    b_temp = 'rcw38-'+pix_info[bolo]['pixel_id'].split('.')[1]+'GHzW'+pix_info[bolo]['pixel_id'].split('w')[1].split('_')[0]
    template_spl = scipy.interpolate.RectBivariateSpline(numpy.arange(templates[b_temp].shape[0]),numpy.arange(templates[b_temp].shape[1]),templates[b_temp])
    a1[bif] = np.ravel(template_spl(np.arange(templates[b_temp].shape[0]) + pix_info[bolo]['y_offset'], numpy.arange(templates[b_temp].shape[1]) + pix_info[bolo]['x_offset'])*invnoise)
    bif +=1
#print('Starting fit on %s bolos'%bif)                                                                  
for i in a1:
    print(i)
co_noise_is = np.linalg.lstsq((a1/1e10).T,(np.ravel(unwmap*invnoise)+np.ravel(a1_neg/1e10)))



#starting fit to the noise region with the injected source (only main bolometer)
bif = 0
template_spl = scipy.interpolate.RectBivariateSpline(numpy.arange(templates[mb_temp].shape[0]),numpy.arange(templates[mb_temp].shape[1]),templates[mb_temp])
len_empty = list(use_bolos)
a1 = np.zeros([1, len(np.ravel(template_spl(np.arange(templates[mb_temp].shape[0]) - pix_info[mb]['y_offset'], numpy.arange(templates[mb_temp].shape[1]) - pix_info[mb]['x_offset'])))])
print('This is the order', sorted((len_empty)))
for i, bolo in enumerate(sorted(len_empty)):
    if bolo != mb:
        continue
    b_temp = 'rcw38-'+pix_info[bolo]['pixel_id'].split('.')[1]+'GHzW'+pix_info[bolo]['pixel_id'].split('w')[1].split('_')[0]
    template_spl = scipy.interpolate.RectBivariateSpline(numpy.arange(templates[b_temp].shape[0]),numpy.arange(templates[b_temp].shape[1]),templates[b_temp])
    a1[bif] = np.ravel(template_spl(np.arange(templates[b_temp].shape[0]) + pix_info[bolo]['y_offset'], numpy.arange(templates[b_temp].shape[1]) + pix_info[bolo]['x_offset'])*invnoise)
    bif +=1
print('Starting fit on %s bolos'%bif)
co_noise_main_bolo_is = np.linalg.lstsq((a1/1e10).T,np.ravel(unwmap*invnoise)+np.ravel(a1_neg/1e10))

other_places = {}
other_places['on_source'] = co
other_places['on_source_only_main'] = co_main_bolo
other_places['off_source'] = co_noise
other_places['off_source_only_main'] = co_noise_main_bolo
other_places['off_source_only_main_injected_source'] = co_noise_main_bolo_is
other_places['off_source_injected_source'] = co_noise_is
other_places['other_bolos_order'] = sorted(len_empty)
other_places['location'] = params_amps

'''
print(co[0])
print(sorted(len_empty))
a = other_places
for i, bolo in enumerate(a['other_bolos_order']):
                lab = (100*(a['on_source'][0][i]/np.max(np.abs(a['on_source'][0]))))
                print(bolo)
                print(params_amps[bolo]['x_offset'],params_amps[bolo]['y_offset'], lab)
'''                
with open(args.output, 'wb') as handle:
    pickle.dump(other_places, handle, protocol=pickle.HIGHEST_PROTOCOL)

