# bolo_matching.py
#
# This module is a collection of functions for verifying the hardware map
# and identification of bolometers.
# Most used are  plot_bolo_offsets() and plot_netanals_and_freq_match()
#
# DPD 2018

import os
from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cPickle as pkl
from spt3g import core,dfmux,gcp,calibration
from pydfmux.spt3g.hwm_tools import wafer_bolo_info
import pdb


# number in offset file is in rad
angle_per_mm = 4.186 * core.G3Units.deg / 1000 

code_repo_dir = '/home/aharkehosemann/repos/' # Change this as needed


def find_bad_pixels_pointing_offset(obs_type,obs_id,plot_focal_plane=False,
				    return_all=False,custom_nom=None):
	'''
	Helps identify mis-matched bolometers by comparing their nominal and 
	offline pointing-derived offsets.
	Given a source observation with offline pointing offsets, this function
	will return a dictionary of information on pixels that have supposed
        member bolometers mapped more than one pixel spacing away.
	
        This function is good at finding pixels with mis-matched bolometers,
        but fails at identifying whole pixels that are mapped to the wrong
        location, and does not by itself identify which bolos are the bad ones.
	
	Parameters
	----------
	obs_type: str
		'RCW38-pixelraster', 'saturn', etc.
	obs_id: str or int
		The id of the observation of interest. e.g. 32009821
	plot_focal_plane: bool
		Set to True to plot the X,Y pointing offsets of the bolos 
		in the obs.
		Use to judge quality of observation
	return_all: bool
		If False, return data will be bad_pixel dict
		If True, return data will be good_pixels dict, bad_pixels dict
		
	Returns
	-------
	bad_pixels: dict
		Dictionary of pixels (wafer_id+'_'+pixel_id) with at least one
		member bolo showing up at a different pixel's location. 
		{pixel:{bolo physical name:{'offsets':(x,y),
		       'nominal':(x,y),'bolo':bolo permanent name}}}
	good_pixels: dict
		Only returned if return_all = True.
		Same form as above, but remaining good pixels.

        addendum 06/2018 AHHH : added functionality with a custom nominal bolometer 
        properties dictionary to find bad pixels after a hardware map rebuild.

	'''
	if not isinstance(obs_id,str):
		obs_id = str(obs_id)
	obs_id = obs_id.replace('.g3','')
	
	if obs_type == 'mars-pixelraster':
		pointing = glob('/spt/user/production/calibration/calframe/mars-pixelraster/'
                                   +obs_id+'.g3')
	else:
		pointing = glob('/spt/user/production/calibration/' +obs_type+ '/'
				   +obs_id+'.g3')
	
	if len(pointing) == 0:
		print("Couldn't find requisite obs data. Exiting")
		return
	
	# Grab nominal bolo physical names, store in dict

	if custom_nom:
		nom = {}
		nom['NominalBolometerProperties'] = custom_nom
		
		phys_names = dict()
		for bolo in nom['NominalBolometerProperties'].keys():
			bolo_name=nom['NominalBolometerProperties'][bolo]['physical_name'].split('/')[-1]
			phys_names[bolo] = bolo_name
		
			
	else:
		nominal_online_cal = glob('/spt/data/bolodata/downsampled/'+obs_type+ 
				  '/'+obs_id+'/nominal_online_cal.g3')
		nom = list(core.G3File(nominal_online_cal[0]))[0]
		
		phys_names = dict()
		for bolo in nom['NominalBolometerProperties'].keys():
			bolo_name=nom['NominalBolometerProperties'][bolo].physical_name
			phys_names[bolo] = bolo_name
			

	# Grab the pointing offsets
	# x,y lists will speed up plotting
	off = list(core.G3File(pointing[0]))[0]
	good_pixels = dict()
	x_offsets = []
	y_offsets = []
	new_bolos = []
	
	print(str(len(off['PointingOffsetY']))+' bolos in this observation')
	for bolo in off['PointingOffsetY'].keys():
		try: 
			pixel = phys_names[bolo].split('.')[0]
			if pixel not in good_pixels.keys():
				good_pixels[pixel] = dict()
			good_pixels[pixel][phys_names[bolo]] = dict()
			# Put pointing offsets into the dictionary
			good_pixels[pixel][phys_names[bolo]]['offsets'] = (
				off['PointingOffsetX'][bolo]/angle_per_mm, 
				off['PointingOffsetY'][bolo]/angle_per_mm)
			# Put nominal offsets into the dictionary
			if custom_nom:
				good_pixels[pixel][phys_names[bolo]]['nominal'] = (
					nom['NominalBolometerProperties'][bolo]['x_offset']/angle_per_mm, 
					nom['NominalBolometerProperties'][bolo]['y_offset']/angle_per_mm)
				good_pixels[pixel][phys_names[bolo]]['lc_channel'] = nom['NominalBolometerProperties'][bolo]['lc_ch']

			else:
				good_pixels[pixel][phys_names[bolo]]['nominal'] = (
					nom['NominalBolometerProperties'][bolo].x_offset/angle_per_mm, 
					nom['NominalBolometerProperties'][bolo].y_offset/angle_per_mm)
			# Add the bolo 'permanent' name to the dictionary, useful for 
			# comparing to metaHWM.csv
			good_pixels[pixel][phys_names[bolo]]['bolo'] = bolo
			x_offsets.append(off['PointingOffsetX'][bolo]/angle_per_mm)
			y_offsets.append(off['PointingOffsetY'][bolo]/angle_per_mm)

		except KeyError:
			new_bolos.append(bolo)

	if plot_focal_plane:
		# Flip the y-offsets to get a more familiar view of focal plane
		plt.plot(x_offsets,-1*np.array(y_offsets),'b.')
		plt.title('Y-Flipped to match view from back port')
		
	# Find bolos not with their nominal pixel
	bad_pixels=dict()
	for pixel in good_pixels.keys():
		# Use first bolo in a pixel to be the ref coords for that pixel
		ref=good_pixels[pixel][good_pixels[pixel].keys()[0]]['offsets']
		# Loop through all bolos in a pixel, stop if too far away
		for bolo in good_pixels[pixel].keys():
			data = good_pixels[pixel][bolo]
			x_diff = data['offsets'][0]-ref[0]
			y_diff = data['offsets'][1]-ref[1]
			dist = np.sqrt(x_diff**2+y_diff**2)
			# Nom pixel separation is 0.0005, wafer radius is 0.005
			if dist > 0.00025/angle_per_mm and dist < 0.004/angle_per_mm:
				# These bolos are believably mismapped
				bad_pixels[pixel]=good_pixels.pop(pixel)
				break
			elif dist > 0.004/angle_per_mm:
				# These bolos are way off, don't deserve to be
				# in good pixels, but not important enough for
				# bad_pixels
				good_pixels.pop(pixel)
				break
	
	#if len(new_bolos) != 0:
	#	print 'New bolos: ', new_bolos

	if return_all:
		return good_pixels, bad_pixels
	return bad_pixels

def plot_bolo_offsets(obs_type,obs_id,label_pix=False,label_ch=False,plot_coord_fits=False,custom_nom=None,scale_cap=0.0041/angle_per_mm, date=None):
	'''
	Makes scatter plot of the focal plane with bolometers colored by
	the difference between their predicted and pointing-derived positions.

	Parameters
	----------
	obs_type: str
	    'RCW38-pixelraster' , etc.
	obs_id: str or int
	    The id of the observation of interest, e.g. 32009821
	label_pix: bool
	    If True, will label the pixels in the scatter plot that do not have
	    any member bolometers located greater than a pixel-spacing apart. 
	plot_coord_fits: bool
	    If True, will generate two additional scatter plots showing the 
	    relation between nominal and pointing-derived X and Y coords
	    for all bolometers.

        addendum 06/2018 AHHH : added channel labelling functionality, and now
        accepts a custom nominal bolometer properties dictionary for post-hwm 
        rebuild viewing pleasure.

	'''

	if not isinstance(obs_id, str):
		obs_id = str(obs_id)

	# Gee, why don't you just include this plot in
	# find_bad_pixels_pointing_offset() ?
	# 1) The whole process was slower
	# 2) The plot didn't end up looking right.
	# Confusing, but if it ain't broke, don't fix it.
	good_pixels, bad_pixels = find_bad_pixels_pointing_offset(
		                     obs_type, obs_id,
                                     plot_focal_plane = False, return_all=True, custom_nom=custom_nom)

	# Lot of plotting, speed it up with these lists
	x_off = []
	y_off = []
	x_nom = []
	y_nom = []
	diff = []

	pixel_coords=dict()
	plot_pix_labels=True

	plt.figure(1)
	plt.clf()
	for pixel in good_pixels.keys():
		
		ref=good_pixels[pixel][good_pixels[pixel].keys()[0]]['offsets']
		pixel_coords[pixel] = ref
		if label_pix:
			plt.figure(1)
			plt.text(ref[0],-1*ref[1],pixel.split('_')[-1],
				 fontsize=8)
		if label_ch:
			plt.figure(1)
			for gbolo in good_pixels[pixel].keys():
				xx_off, yy_off=good_pixels[pixel][gbolo]['offsets'] 
				plt.text(xx_off,-1*yy_off-0.04,good_pixels[pixel][gbolo]['lc_channel'],
                                 fontsize=8)
		for bolo, data in good_pixels[pixel].iteritems():
			x_off.append(data['offsets'][0])
			y_off.append(data['offsets'][1])
			x_nom.append(data['nominal'][0])
			y_nom.append(data['nominal'][1])
			

	for pixel in bad_pixels.keys():
		for bolo, data in bad_pixels[pixel].iteritems():
			x_off.append(data['offsets'][0])
			y_off.append(data['offsets'][1])
			x_nom.append(data['nominal'][0])
			y_nom.append(data['nominal'][1])

		ref=bad_pixels[pixel][bad_pixels[pixel].keys()[0]]['offsets']
                pixel_coords[pixel] = ref
                if label_pix:
                        plt.figure(1)
                        plt.text(ref[0],-1*ref[1],pixel.split('_')[-1],
                                 fontsize=8)
   
		if label_ch:
			plt.figure(1)
			for bbolo in bad_pixels[pixel].keys():
				xx_off, yy_off=bad_pixels[pixel][bbolo]['offsets'] 
				plt.text(xx_off,-1*yy_off-0.04,bad_pixels[pixel][bbolo]['lc_channel'],
                                 fontsize=8)
        
	x_fit = np.polyfit(x_nom,x_off,1)
	y_fit = np.polyfit(y_nom,y_off,1)

	if plot_coord_fits:
		plt.figure(2)
		plt.clf()
		plt.plot(x_nom,x_off,'b.')
		x_hat = np.linspace(-.025,.025,100)
		plt.plot(x_hat,x_fit[0]*x_hat+x_fit[1],'k-')
		plt.title('X')
		plt.xlabel('Nominal offsets')
		plt.ylabel('Pointing-derived offsets')

		plt.figure(3)
		plt.clf()
		plt.plot(y_nom,y_off,'b.')
		y_hat = np.linspace(-.025,.025,100)
		plt.plot(y_hat,y_fit[0]*y_hat+y_fit[1],'k-')
		plt.title('Y')
		plt.xlabel('Nominal offsets')
		plt.ylabel('Pointing-derived offsets')

	for ind in range(len(x_off)):
		x_diff = x_off[ind] - (x_nom[ind]*x_fit[0]+x_fit[1])
		y_diff = y_off[ind] - (y_nom[ind]*y_fit[0]+y_fit[1])
		dist = np.sqrt(x_diff**2+y_diff**2)
		diff.append(dist)
		

	plt.figure(1)
	plt.scatter(x_off,-1*np.array(y_off),c=diff,vmin=0,vmax=scale_cap,
		    cmap='inferno_r',alpha=.5,edgecolor='None')
	plt.colorbar(shrink=0.5,label='Dist. from predicted location [mm]')
	plt.ylim(-.02/angle_per_mm,.02/angle_per_mm)
	plt.xlim(-.025/angle_per_mm,.025/angle_per_mm)
	plt.xlabel('X-Offset [mm]')
	plt.ylabel('Y-Offset [mm]')

	if date == None:
		title = obs_type+' '+obs_id
	else:
		title = obs_type+' '+date
	plt.title(title)

	"""
	### find one channel
	for pixel in bad_pixels.keys():
		if pixel == 'w181_42':
			print 'found pixel'
			for bolo in good_pixels[pixel].keys():
				ax = plt.gca()
				offsets = good_pixels[pixel][bolo]['offsets']
				ax.plot(offsets[0],-1*offsets[1], 'r.')
	"""

def find_certain_correct_channels(good_pixels,bad_pixels):
	'''
	Compares the pointing-derived offset of possibly bad bolos
	with 'perfect' pixels. If a bolo is located within a different pixel's
	footprint, the frequency schedule and metahwm are consulted to see how 
	such a mix-up could have occurred.

	Takes as input the results of find_bad_pixels_pointing_offset,
	with return_all=True

	Conclusions are printed to screen, and returned in a dictionary
	Not comprehensive or fool-proof
	'''
	
	pixel_coords = dict()
	bad_bolos = dict()

	for pixel in good_pixels.keys():
		ref=good_pixels[pixel][good_pixels[pixel].keys()[0]]['offsets']
		pixel_coords[pixel] = ref

	for pixel in bad_pixels.keys():
		for bolo, data in bad_pixels[pixel].items():
			for pix, coords in pixel_coords.items():
				x_diff = data['offsets'][0] - coords[0]
				y_diff = data['offsets'][1] - coords[1]
				dist = np.sqrt(x_diff**2+y_diff**2)
				if dist < 0.00025/angle_per_mm:
					bad_bolos[bolo] = dict()
					bad_bolos[bolo]['bolo'] = bad_pixels[pixel][bolo]['bolo']
					bad_bolos[bolo]['offset_pix'] = pix
	
	mapping = wafer_bolo_info()
	unclear=[]

	for phys_name, data in bad_bolos.items():
		bad_bolos[phys_name]['module']=data['bolo'].rpartition('.')[0]
		name = phys_name.split('_')[-1]
		ch=int(mapping.loc[mapping['physical_name']==name]['lc_ind']+1)
		right_pix = int(data['offset_pix'].split('_')[-1])
		try:
			right_ch=int(mapping.loc[(mapping['pixel'] == right_pix) & 
						 (abs(mapping['lc_ind']+1-ch)== 1)]['lc_ind']+1)
		except TypeError:
		#	print(phys_name+' is just confused')
			unclear.append(phys_name)
			continue
		bad_bolos[phys_name]['ch'] = ch
		bad_bolos[phys_name]['correct_ch'] = right_ch

	metahwm_path = (code_repo_dir+'hardware_maps_southpole/2018/2018_global/metaHWM.csv')
	metahwm = pd.read_csv(metahwm_path, delimiter = '\t')

	for bolo in unclear:
		bad_bolos.pop(bolo)

	for bolo in bad_bolos.keys():
		mod = bad_bolos[bolo]['module']
		# bolo permanent names changed
		if '/' in mod:
			mod = mod.split('/')[-1]
		crate, slot, mezz, mod = mod.split('.')
		lc_chip = metahwm.loc[(metahwm['crate']==int(crate)) & 
				      (metahwm['slot']==int(slot)) &
				      (metahwm['mezzanine']==int(mezz)) & 
				      (metahwm['module']==int(mod))]['lc_chip']
		lc_chip = lc_chip.values[0]
		bad_bolos[bolo]['lc_chip'] = lc_chip

	for bolo in bad_bolos.keys():
		print(bolo +' should be in '+bad_bolos[bolo]['offset_pix']+'. '
		      'Currently Ch'+str(bad_bolos[bolo]['ch'])+' on '+
		      bad_bolos[bolo]['lc_chip'] +
		      ', should be Ch'+str(bad_bolos[bolo]['correct_ch']))

	return bad_bolos

def plot_netanals_and_freq_match(chip=None,pstring=None):
	'''
	Overplots cold and warm (if extant) netanals with the lc channels
	identified by the lc matching algorithm.

	Can either pass the LC chip name, e.g. 'LC68.v1.b3.c4'
	or the module pstring, e.g. '005/5/2/3'
	'''

	if not (chip or pstring):
		raise ValueError('Invalid number of arguments')
	
	if pstring:
		metahwm_path = (code_repo_dir+'hardware_maps_southpole/'
				'2018/2018_global/metaHWM.csv')
		metahwm = pd.read_csv(metahwm_path, delimiter = '\t')
		crate,slot,mezz,mod = pstring.split('/')
                chip = metahwm.loc[(metahwm['crate']==int(crate)) &
				      (metahwm['slot']==int(slot)) &
				      (metahwm['mezzanine']==int(mezz)) & 
				      (metahwm['module']==int(mod))]['lc_chip']
		
		assert(len(chip.values)>0),'LC Chip/Module not in metaHWM'
		chip = chip.values[0]
	

	netanal_dir = '/big_scratch/poleanalysis/receiver_commissioning_2017_2018/netanals/'
	warm_pkls = glob(netanal_dir+'warm_good/*OUTPUT.pkl')
	cold_pkls = glob(netanal_dir+'cold_good/*OUTPUT.pkl')

	warm_mods = np.array([path.rsplit('/')[-1] for path in warm_pkls])
	cold_mods = np.array([path.rsplit('/')[-1] for path in cold_pkls])

	f = open(code_repo_dir + '/hardware_maps_southpole/2018/2018_global/3g_lc_matching.pkl','r')
	lc_data = pkl.load(f)
	f.close()

	mod = lc_data[chip]['na_data'].split('/')[-1]

	f = open(cold_pkls[np.where(cold_mods==mod)[0]],'r')
	cold_na = pkl.load(f)
	f.close()

	plt.figure(figsize=(10,4))
	plt.clf()

	amp = cold_na['carrier_NA']['amp']/cold_na['nuller_NA']['amp']
	plt.plot(cold_na['carrier_NA']['freq'], amp, 'b',
		 label = 'Cold NA')
	for i in range(len(lc_data[chip]['optimal_freqs'])):
		label = '%d' % (lc_data[chip]['lc_ind'][i]+1)
		plt.axvline(lc_data[chip]['optimal_freqs'][i], color='k', 
			    linestyle = '--')
		plt.text(lc_data[chip]['optimal_freqs'][i]+1e3, 40, label,
			 fontsize=10,rotation=90, color = 'k')

	if mod in warm_mods:
		f = open(warm_pkls[np.where(warm_mods==mod)[0]])
		warm_na = pkl.load(f)
		f.close()
		amp = warm_na['carrier_NA']['amp']/warm_na['nuller_NA']['amp']
		plt.plot(warm_na['carrier_NA']['freq'], amp, 'r', 
			 label='Warm NA')

	plt.title(chip+', '+mod.replace('_OUTPUT.pkl',''))
	plt.suptitle('Channels numbered 1-68')
	plt.xlim(1.5e6, 6.3e6)
	plt.legend(loc=1,fontsize=10)

def get_pix_info(wafer = None, pixel= None):
	'''
	wafer: e.g. 'w172'
	pixel: e.g. 243
	'''
	if wafer==None or pixel==None:
		print(' Please specify a wafer name, e.g. ''w172"'
		      'and pixel number, e.g. 151')

	metahwm_path = (code_repo_dir+'hardware_maps_southpole/2018/'
			'2018_global/metaHWM.csv')
	metahwm = pd.read_csv(metahwm_path, delimiter = '\t')
	mapping = wafer_bolo_info()

	side=int(mapping.loc[(mapping['pixel']==pixel)]['side'].values[0])
	flex=mapping.loc[(mapping['pixel']==pixel)]['flex_cable'].values[0]
	lc_ch = (mapping.loc[(mapping['pixel']==pixel)]['lc_ind'].values)+1
	# Some pixels split across flex cables, always on same squid
	if flex in [1,2]: pair = '1,2'
	if flex in [3,4]: pair = '3,4'
	if flex in [5,6]: pair = '5,6'
	if flex in [7,8]: pair = '7,8'
	mod = metahwm.loc[(metahwm['wafer'] == wafer) &
			   (metahwm['side'] == side) &
			   (metahwm['flex_cable'] == pair)]

	print('Pixel '+str(pixel)+' has bolos at LC Ch[1-68]:'+ str(lc_ch))
	print('MetaHWM info:')
	print mod
