#!/usr/bin/env python
import numpy, sys, os
import scipy.ndimage, scipy.interpolate, scipy.optimize
from spt3g import core, mapmaker, calibration
import argparse as ap
import numpy as np 
import pickle
import matplotlib.pyplot as plt

P = ap.ArgumentParser(description='Pointing and calibration off of a point source',
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
#P.add_argument('-ip','--input', default=['/spt/data/bolodata/downsampled/RCW38-pixelraster/64584651/offline_calibration.g3', '/spt/user/production/calibration/RCW38-pixelraster/maps/64584651.g3'], help = 'Input Files')

P.add_argument('input_files', action='store', nargs='+', default=[], help='Inpu\
t files')
P.add_argument('--noise', action = 'store_true')
args = P.parse_args()                                                          

a = pickle.load(open(args.input_files[2], 'rb'))

old= args.input_files[2]
old_old = old.split('_fit_noise')[0]+old.split('_fit_noise')[1]+old.split('_fit_noise')[2]
a_real_fit = pickle.load(open(old_old, 'rb'))
obsnum = args.input_files[0].split('/')[-1].split('.')[0]

templates = pickle.load(open(args.input_files[0],'rb'))['final_temp']           
params = pickle.load(open(args.input_files[0],'rb'))['params']                  
rawmap = pickle.load(open(args.input_files[1],'rb'))['data']                    
print(args.input_files[1])
wiring_info = pickle.load(open(args.input_files[0],'rb'))['wiring_info']
w = pickle.load(open(args.input_files[1],'rb'))['w']

main_bolo = sorted(a['150.0'].keys())[0]
group = '150.0'
template_spl = scipy.interpolate.RectBivariateSpline(numpy.arange(templates[group].shape[0]),numpy.arange(templates[group].shape[1]),templates[group])
params_amps = {**params['90.0'],**params['150.0'], **params['220.0']}
invnoise = w**0.5

#check solution
a_check = np.zeros([360,360])
for i, bolo in enumerate(a[group]['other_bolos_order']):
                a_check+= a[group][main_bolo][0][i]*(template_spl(np.arange(templates[group].shape[0])- params_amps[bolo][1], numpy.arange(templates[group].shape[1]) - params_amps[bolo][2])*invnoise)
                
                                
a_splash = np.zeros([360,360])
for i, bolo in enumerate(a[group]['other_bolos_order']):
                a_splash+= (template_spl(np.arange(templates[group].shape[0])- params_amps[bolo][1], numpy.arange(templates[group].shape[1]) - params_amps[bolo][2])*invnoise)

if args.noise:
                print('doing for a noise fit')
                a_check = np.zeros([360,360])
                for i, bolo in enumerate(a[group]['other_bolos_order']):
                                a_check+= a[group][main_bolo][0][i]*(template_spl(np.arange(templates[group].shape[0])+ params_amps[bolo][1], numpy.arange(templates[group].shape[1]) + params_amps[bolo][2])*invnoise)
                                
                                
                a_splash = np.zeros([360,360])
                for i, bolo in enumerate(a[group]['other_bolos_order']):
                                    a_splash+= (template_spl(np.arange(templates[group].shape[0])+ params_amps[bolo][1], numpy.arange(templates[group].shape[1]) + params_amps[bolo][2])*invnoise)

            
coeffs = [i for i in a[group][main_bolo][0] if i>0]
if args.noise:
                coeffs_max = np.max([i for i in a_real_fit[group][main_bolo][0] if i>0])
else:
                coeffs_max = np.max(coeffs)
                print('not a noise run')
print('found the max %s'%len(coeffs))
try:
                b = pickle.load(open('full_notnoise.pkl','rb'))
except:
                b = {}
b[args.input_files[2]] = coeffs/coeffs_max
print(coeffs/coeffs_max)
import pickle
with open('full_notnoise.pkl', 'wb') as handle:
                pickle.dump(b, handle, protocol=pickle.HIGHEST_PROTOCOL)

if args.noise:
                                    obsnum = 'noise_'+obsnum

'''
print('making histogram')
plt.hist((coeffs/coeffs_max), bins = 100)
plt.ylabel('Number of Bolometers')
plt.xlabel('Percent Crosstalk')
plt.title('%s'%main_bolo) 
print('trying to save')
plt.xlim(0,1)
plt.savefig('/spt/user/javva/crosstalk/plots/histograms/%s_%s.png'%(obsnum,main_bolo))
plt.xlim(0,.1)
plt.title('zoom %s'%main_bolo) 
plt.savefig('/spt/user/javva/crosstalk/plots/histograms/zoom/%s_%s.png'%(obsnum,main_bolo))
plt.clf()
print('checking fit place')
plt.imshow(np.abs(a_check)/np.max(np.abs(a_check)))
plt.colorbar()
plt.clim(0,.1)
plt.title('Fit Coefficients * Shifted Template normalized to maximum')
plt.savefig('/spt/user/javva/crosstalk/plots/plot/%s_%s.png'%(obsnum,main_bolo))
plt.clf()
plt.imshow(a_splash)
plt.title('Splash of where the fit happens - %s bolos'%len(a[group]['other_bolos_order']))
plt.savefig('/spt/user/javva/crosstalk/plots/splash/%s_%s.png'%(obsnum,main_bolo))
plt.clf()
rawmap[numpy.logical_not(numpy.isfinite(rawmap))] = 0
plt.imshow(np.abs(rawmap)/np.max(np.abs(rawmap)))
print(np.where(rawmap>0))
print(np.max(np.abs(rawmap)))
plt.title('Map (data) normalized to map max %s'%main_bolo)
plt.colorbar()
plt.clim(0,.1)
#plt.savefig('/spt/user/javva/crosstalk/plots/data2/%s_%s.png'%(obsnum,main_bolo))
num_over = 0

if args.noise:
    print('doing scatter with noise params')
    for i, bolo in enumerate(a[group]['other_bolos_order']):
                    if bolo == main_bolo:
                                     plt.scatter((360./2)-params_amps[bolo][2],(360./2)-params_amps[bolo][1], marker = 'x', c = 'k')

                    if a[group][main_bolo][0][i]/coeffs_max > 0.15:
                                    plt.scatter((360./2)-params_amps[bolo][2],(360./2)-params_amps[bolo][1], marker = 'x', c = 'g')
                                    num_over += 1
else:
    for i, bolo in enumerate(a[group]['other_bolos_order']):
                    if bolo == main_bolo:
                                     plt.scatter((360./2)+params_amps[bolo][2],(360./2)+params_amps[bolo][1], marker = 'x', c = 'k')

                    if a[group][main_bolo][0][i]/np.max(coeffs) > 0.15:
                                    plt.scatter((360./2)+params_amps[bolo][2],(360./2)+params_amps[bolo][1], marker = 'x', c = 'g')
                                    num_over += 1

plt.title('# over 15 percent = %s for bolo %s'%(num_over, main_bolo))
plt.savefig('/spt/user/javva/crosstalk/plots/data2/x_overlay_more_15percent%s_%s.png'%(obsnum,main_bolo))
plt.clf()

plt.imshow((np.abs(rawmap)-np.abs(a_check))/(np.max(np.abs(rawmap))))
plt.title('Residual -- (map - fit)')
plt.colorbar()
plt.savefig('/spt/user/javva/crosstalk/plots/residual/%s_%s.png'%(obsnum,main_bolo))
'''

