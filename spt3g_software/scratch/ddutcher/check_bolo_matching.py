# check_bolo_matching.py
#
# This module is a collection of functions for verifying the hardware map
# and identification of bolos using measured bolometer pointing offsets.
#
# DPD 2018

import os
from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from spt3g import core,dfmux,gcp,calibration
from spt3g.util import files
from pydfmux.spt3g.hwm_tools import wafer_bolo_info

metahwm_path = os.path.join('/home/ddutcher/code',  # Change this as needed
                            'hardware_maps_southpole/2019/',
                            'global/metaHWM_2019.csv')

def check_bolo_matching(obsid, minnum_bolos=8000, verbose = True,
                        find_correct_ch=False, hwm_dir = None,
                        boloprops = None, return_data=False, 
                        plot_pointing_offsets=True, **plot_kwargs):
    '''
    Identify mismatched bolometers by comparing their nominal and 
    offline pointing-derived offsets.
    Given a source observation with offline pointing offsets, this function
    will return a dictionary of information on pixels that have supposed
    member bolometers mapped more than one pixel spacing away.
    
    This function is good at finding pixels with mis-matched bolometers,
    but fails at identifying whole pixels that are mapped to the wrong
    location.
    
    Parameters:
    -----------
    obsid: str or int
        The observation id of the pixelraster observation of interest.
    minnum_bolos: int
        Don't consider pixelrasters with fewer bolometers than this.
    find_correct_ch: bool
        If a bolo's offsets locate it within a different pixel's
        footprint, the frequency schedule and hwm are consulted to see how 
        such a mix-up could have occurred, and findings are printed to screen.
            * Check these results for yourself!!
        Will also check your metahwm for squid swaps.
    hwm_dir: str
        The path to the relevant hardware map directory. Required for 
        find_correct_ch to function.
    boloprops: G3File, G3Frame, or BolometerPropertiesMap
            If not None, the script will use these Nominal Bolometer
            Properties instead of the on-disk nominal_online_cal
    return_data: bool
        Return data or not
    plot_pointing_offsets: bool
        Plot the pointing offsets for each bolometer, colored by distance
        from their predicted location.
        
    **plot_kwargs: all further arguments passed to plot_bolo_offsets
        
    Returns:
    --------
    Dictionary with the following keys:
        bad_pixels: dict
            Dictionary of pixels (wafer_id+'_'+pixel_id) with at least one
            member bolo showing up at a different pixel's location. 
            {pixel:{bolo0 physical name:{'offsets':(x,y),
                   'nominal':(x,y),'bolo':bolo0 permanent name} ,
                   bolo1 physical name:{'offsets':(x,y),
                   'nominal':(x,y),'bolo':bolo1 permanent name} , ... }}
        good_pixels: dict
            Same form as above, but remaining good pixels.
        bolo_fixes: dict
            Included if find_correct_ch=True. Contains info on fixable bolos.
    '''    
    datapath = '/spt/data/bolodata/downsampled/'
    calpath = '/spt/user/production/calibration/'
    if not isinstance(obsid,str):
        obsid = str(obsid)
    obsid = os.path.basename(obsid)
    obsid = obsid.replace('.g3','')
    
    test = glob(os.path.join(datapath,'*',obsid))
    if len(test)==0:
        print('Could not find any observation with obsid %s'%obsid)
        return
    else:
        source = test[0].split('/')[-2]
        if '-pixelraster' not in source:
            raise TypeError("'-pixelraster' observation required.")
        
    if not os.path.isfile(os.path.join(calpath,source,obsid+'.g3')):
        print('Could not find pointing fits for %s %s. Did autoprocessing run?'
             %(source, obsid))
        return
    pointing = os.path.join(calpath,source,obsid+'.g3')
    ######################################################
    ### Section dealing with possible custom boloprops ###
    nbp = None
    if boloprops is None:
        nominal_online_cal = os.path.join(datapath,source,obsid,
                                          'nominal_online_cal.g3')    
        nbp = list(core.G3File(nominal_online_cal))[0]
        nbp = nbp['NominalBolometerProperties']
    elif isinstance(boloprops, str):
        if boloprops.split('.')[-1] != 'g3':
            raise OSError('boloprops file must be a .g3 file')
        tmp = core.G3File(boloprops)
        for fr in tmp:
            if 'NominalBolometerProperties' in fr:
                nbp = fr['NominalBolometerProperties']
        assert nbp is not None
    elif isinstance(boloprops, core.G3Frame):
        nbp = fr['NominalBolometerProperties']
    elif isinstance(boloprops, calibration.BolometerPropertiesMap):
        nbp = boloprops
    else:
        raise ValueError('Custom bolometer properties format not recognized')
    #######################################################
    
    # Grab nominal bolo physical names, store in dict
    phys_names = dict()
    
    for bolo in nbp.keys():
        bolo_name=nbp[bolo].physical_name
        phys_names[bolo] = bolo_name
            
    # Grab the pointing offsets
    # x,y lists will speed up plotting
    off = list(core.G3File(pointing))[0]
    if verbose:
        print(str(len(off['PointingOffsetY']))+' bolos in this observation')
    if len(off['PointingOffsetY']) < minnum_bolos:
        print('Fewer than minnum_bolos=%s in observation. Quitting...'
             %minnum_bolos)
        return
    
    good_pixels = dict()
    x_offsets = []
    y_offsets = []
    
    for bolo in off['PointingOffsetY'].keys():
        pixel = phys_names[bolo].split('.')[0]
        if pixel not in good_pixels.keys():
            good_pixels[pixel] = dict()
        good_pixels[pixel][phys_names[bolo]] = dict()
        # Put pointing offsets into the dictionary
        good_pixels[pixel][phys_names[bolo]]['offsets'] = (
            off['PointingOffsetX'][bolo], 
            off['PointingOffsetY'][bolo])
        # Put nominal offsets into the dictionary
        good_pixels[pixel][phys_names[bolo]]['nominal'] = (
            nbp[bolo].x_offset, nbp[bolo].y_offset)
        # Add the bolo 'permanent' name to the dictionary, useful for 
        # comparing to metaHWM.csv
        good_pixels[pixel][phys_names[bolo]]['bolo'] = bolo
        x_offsets.append(off['PointingOffsetX'][bolo])
        y_offsets.append(off['PointingOffsetY'][bolo])
        
    if plot_pointing_offsets:
        plot_bolo_offsets(obsid, boloprops = nbp, **plot_kwargs)
        
    # Find bolos not with their nominal pixel
    bad_pixels=dict()
    for pixel in list(good_pixels.keys()):
        # Use first bolo in a pixel to be the ref coords for that pixel
        ref=good_pixels[pixel][list(good_pixels[pixel].keys())[0]]['offsets']
        # Loop through all bolos in a pixel, stop if too far away
        for phys_name in list(good_pixels[pixel].keys()):
            data = good_pixels[pixel][phys_name]
            x_diff = data['offsets'][0]-ref[0]
            y_diff = data['offsets'][1]-ref[1]
            dist = np.sqrt(x_diff**2+y_diff**2)
            # Nom pixel separation is 0.0005, wafer radius is 0.005
            if dist > 0.00025 and dist < 0.005:
                # These bolos are believably mismapped
                bad_pixels[pixel]=good_pixels.pop(pixel)
                break
            elif dist > 0.005:
                # These bolos are way off, don't deserve to be
                # in good pixels, but not important enough for
                # bad_pixels
                good_pixels.pop(pixel)
                break
    return_d = {'good_pixels':good_pixels,'bad_pixels':bad_pixels}

    if find_correct_ch:
        if hwm_dir is None:
            hwm_dir = os.path.dirname(os.path.dirname(metahwm_path))
            hwm_dir = os.path.join(hwm_dir,'hwm_pole')
        if os.path.exists(hwm_dir):
            if verbose:
                print('Using hwm at %s to attempt to find corrections'%hwm_dir)
            bolo_fixes = find_correct_channels(
                good_pixels, bad_pixels, hwm_dir, 
                return_data = True, verbose = verbose)
            return_d['bolo_fixes'] = bolo_fixes
            
    if return_data:
        return return_d
    
def plot_bolo_offsets(obsid, label_bolos = False, wafers = None,
                      print_fit=False, boloprops=None,
                      vmax = 12):
    '''
    Takes the pointing offsets from '-pixelraster' observations of sources
    and compares them to their predicted values* from the hardware map.
    The pointing offsets for each bolometer are then plotted,
    colored by distance from their predicted location.
    
    *A linear relationship is fit for when comparing the pointing offsets and 
     the hardware map values for bolometer positions.
    
    Paramters:
    ----------
    obsid: str or int
            The observation ID of the pixelraster
    verbose: bool
            prints more info to screen
    label_bolos: boolean
            labels each bolometer in the focal plane plot
    wafers: str or list of str
            If label_bolos=True then only the wafers included in `wafers`
            will have their bolometers labelled.
    boloprops: G3File, G3Frame, or BolometerPropertiesMap
            If not None, the script will use these Nominal Bolometer
            Properties instead of the on-disk nominal_online_cal
    vmax [12]: int
            The upper limit of the colorscale in the plot, in arcminutes.   
    '''
    datapath = '/spt/data/bolodata/downsampled/'
    calpath = '/spt/user/production/calibration/'
    if not isinstance(obsid,str):
        obsid = str(obsid)
    obsid = os.path.basename(obsid)
    obsid = obsid.replace('.g3','')    
    
    test = glob(os.path.join(datapath,'*',obsid))
    if len(test)==0:
        print('Could not find any observation with obsid %s'%obsid)
        return
    else:
        source = test[0].split('/')[-2]
        if '-pixelraster' not in source:
            raise TypeError("'-pixelraster' observation required.")
        
    if not os.path.isfile(os.path.join(calpath,source,obsid+'.g3')):
        print('Could not find pointing fits for %s %s. Did autoprocessing run?'
             %(source, obsid))
        return

    ######################################################
    ### Section dealing with possible custom boloprops ###
    nbp = None
    if boloprops is None:
        nominal_online_cal = os.path.join(datapath,source,obsid,
                                          'nominal_online_cal.g3')    
        nbp = list(core.G3File(nominal_online_cal))[0]
        nbp = nbp['NominalBolometerProperties']
    elif isinstance(boloprops, str):
        if boloprops.split('.')[-1] != 'g3':
            raise OSError('boloprops file must be a .g3 file')
        tmp = core.G3File(boloprops)
        for fr in tmp:
            if 'NominalBolometerProperties' in fr:
                nbp = fr['NominalBolometerProperties']
        assert nbp is not None
    elif isinstance(boloprops, core.G3Frame):
        nbp = fr['NominalBolometerProperties']
    elif isinstance(boloprops, calibration.BolometerPropertiesMap):
        nbp = boloprops
    else:
        raise ValueError('Custom bolometer properties format not recognized')
    #######################################################
    pointing = list(core.G3File(os.path.join(calpath,source,obsid+'.g3')))[0]
    
    hwm_x,hwm_y,pointing_x,pointing_y, diff = [],[],[],[],[]
    physical_names = []
    
    for bolo in pointing['PointingOffsetX'].keys():
        if bolo not in nbp:
            print('%s has no bolometer properties'%bolo)
            continue
        if not np.isfinite(pointing['PointingOffsetX'][bolo]):
            continue
        if not np.isfinite(pointing['PointingOffsetY'][bolo]):
            continue
        pointing_x.append(pointing['PointingOffsetX'][bolo])
        pointing_y.append(pointing['PointingOffsetY'][bolo])
        
        hwm_x.append(nbp[bolo].x_offset)
        hwm_y.append(nbp[bolo].y_offset)
        physical_names.append(nbp[bolo].physical_name)     
        
    # Obtain a mapping from hwm coordinates -> pointing offsets
    # A linear fit is sufficient
    # Don't fit nans
    hwm_x = np.array(hwm_x)
    hwm_y = np.array(hwm_y)
    pointing_x = np.array(pointing_x)
    pointing_y = np.array(pointing_y)
    
    x_fit = np.polyfit(hwm_x, pointing_x, 1)
    y_fit = np.polyfit(hwm_y, pointing_y, 1)
    
    if print_fit:
        print('X fit: m = %.2e, b = %.2e'%(x_fit[0],x_fit[1]))
        print('Y fit: m = %.2e, b = %.2e'%(y_fit[0],y_fit[1]))
    for ind in range(len(pointing_x)):
        x_diff = pointing_x[ind] - (hwm_x[ind]*x_fit[0]+x_fit[1])
        y_diff = pointing_y[ind] - (hwm_y[ind]*y_fit[0]+y_fit[1])
        mag = np.sqrt(x_diff**2 + y_diff**2)
        diff.append(mag)    
    
    # I flip the y-units to match my preferred focal plane orientation
    pointing_y = -1*np.array(pointing_y)/core.G3Units.arcmin
    pointing_x = np.array(pointing_x)/core.G3Units.arcmin
    diff = np.array(diff)/core.G3Units.arcmin
    
    # sorting is crucial to have the bad bolos stand out
    ind = np.argsort(diff)
    
    plt.figure()
    if label_bolos:
        for idx, bolo in enumerate(physical_names):
            if wafers is not None:
                if bolo.split('_')[0] in wafers:
                    plt.text(pointing_x[idx],pointing_y[idx],bolo, fontsize=8,
                            alpha=0.7)
            else:
                plt.text(pointing_x[idx], pointing_y[idx], bolo,
                         fontsize=8, alpha=0.7)
    plt.scatter(pointing_x[ind], pointing_y[ind],c=diff[ind], s=10,
                vmin=0,vmax=vmax, cmap='inferno_r',alpha=.5,edgecolor='None')
    plt.colorbar(shrink=0.5,label='Dist. from predicted [arcmin]')
    plt.xlabel('arcmin')
    plt.ylabel('arcmin')
    plt.title(source+' '+obsid)
    
def find_correct_channels(good_pixels, bad_pixels, hwm_dir,
                          return_data = False, verbose=True):
    '''
    Compares the pointing-derived offset of possibly bad bolos
    with 'perfect' pixels. If a bolo is located within a different pixel's
    footprint, the frequency schedule and hwm are consulted to see how 
    such a mix-up could have occurred.

    Takes as input the results of find_bad_bolos

    Conclusions are printed to screen, and returned in a dictionary
    Not comprehensive or fool-proof
    '''
    check_for_sq_swaps(metahwm_path)
    
    metahwm = pd.read_csv(metahwm_path, delimiter = '\t')
    if '.yaml' in hwm_dir:
        hwm_dir = os.path.dirname(hwm_dir)
        
    mapping = wafer_bolo_info()
    
    pixel_coords = dict()
    bad_bolos = dict()

    for pixel in good_pixels.keys():
        ref=good_pixels[pixel][list(good_pixels[pixel].keys())[0]]['offsets']
        pixel_coords[pixel] = ref

    for pixel in bad_pixels.keys():
        for phys_name, data in bad_pixels[pixel].items():
            for pix, coords in pixel_coords.items():
                x_diff = data['offsets'][0] - coords[0]
                y_diff = data['offsets'][1] - coords[1]
                dist = np.sqrt(x_diff**2+y_diff**2)
                if dist < 0.00025:
                    bad_bolos[phys_name] = dict()
                    bad_bolos[phys_name]['bolo'] = \
                        bad_pixels[pixel][phys_name]['bolo']
                    bad_bolos[phys_name]['offset_pix'] = pix
    
    unclear=[]
    if verbose:
        print('Recommended changes:')
        print('(Check that these are repeatable before implementing!)')
    for phys_name in list(sorted(bad_bolos.keys())):
        wafer, name = phys_name.split('_')
        try:
            mapping_csv = glob(os.path.join(hwm_dir,'mappings',
                                            'mapping_'+wafer+'.csv'))[0]
        except IndexError:
            print('Could not find mapping file for %s within %s/mappings'
                  %(wafer, hwm_dir))
            continue
        waf_mapping = pd.read_csv(mapping_csv, delimiter='\t')
        row = waf_mapping.loc[
            waf_mapping['bolometer']==wafer+'/'+bad_bolos[phys_name]['bolo']]
        lc_chip, ch = (row['lc_path'].values[0]).split('/')
        correct_pix = bad_bolos[phys_name]['offset_pix']
        try:
            correct_info = mapping.loc[
                (mapping['pixel'] == int(correct_pix.split('_')[-1])) & 
                (abs(mapping['lc_ind']+1-int(ch)) == 1)]
            correct_ch = int(correct_info['lc_ind'])+1
            band = str(int(correct_info['observing_band']))
            xy = correct_info['pol_xy'].values[0]
        except TypeError:
        #    print(phys_name+' is just confused')
            unclear.append(phys_name)
            bad_bolos.pop(phys_name)
            continue
        bad_bolos[phys_name]['lc_chip'] = lc_chip
        bad_bolos[phys_name]['ch'] = int(ch)
        bad_bolos[phys_name]['correct_ch'] = correct_ch
        bad_bolos[phys_name]['correct_physname'] = \
            str(correct_pix)+'.'+band+'.'+xy
        if verbose:
            print(phys_name+' ---> '+bad_bolos[phys_name]['correct_physname'],
                  '\t Currently Ch'+ch+' on '+lc_chip+
                  ', should be Ch'+str(correct_ch))
        
    if return_data:
        return bad_bolos

def plot_netanals_and_freq_match(chip = None, pstring = None,
                                 lc_matching = None,
                                 cold_data_dir = None,
                                 warm_data_dir = None):
    '''
    Overplots cold and warm (if extant) netanals with the lc channels
    identified by the lc matching algorithm.
    Can either pass the LC chip name or the module pstring.
    
    Parameters:
    ----------
    chip: The lc chip name, e.g. 'LC68.v1.b3.c4'
    pstring: Module pstring, e.g. '005/5/2/3'
    lc_matching: Location of the lc_matching pkl used for building
        the hardware map
    cold_data_dir: directory containing the cold NetAnal pkls
    warm_data_dir: directory containing the warm NetAnal pkls
    '''

    if not (chip or pstring):
        raise ValueError("Must assign one of 'chip' or 'pstring' kwargs")
        
    if cold_data_dir is None:
        raise ValueError("Missing required 'cold_data_dir' kwarg")

    if pstring:
        metahwm = pd.read_csv(metahwm_path, delimiter = '\t')
        crate,slot,mezz,mod = pstring.split('/')
        chip = metahwm.loc[(metahwm['crate']==int(crate)) &
                           (metahwm['slot']==int(slot)) &
                           (metahwm['mezzanine']==int(mezz)) & 
                           (metahwm['module']==int(mod))]['lc_chip']
        if len(chip.values)==0:
            raise KeyError('LC Chip/Module not in metaHWM')
        chip = chip.values[0]
        
    cold_pkls = glob(os.path.join(cold_data_dir,'*OUTPUT.pkl'))
    cold_mods = np.array([path.rsplit('/')[-1] for path in cold_pkls])
    
    if warm_data_dir is not None:
        warm_pkls = glob(os.path.join(warm_data_dir,'*OUTPUT.pkl'))
        warm_mods = np.array([path.rsplit('/')[-1] for path in warm_pkls])
    else:
        warm_mods=[]
    if lc_matching is None:
        global_dir = os.path.dirname(metahwm_path)
        lc_matching = glob(global_dir+'/*lc_matching*.pkl')
        if len(lc_matching)!=1:
            print("Could not find lc_matching pickle in expected place.")
            print("You'll probably need to edit this script. . .")
            return
        lc_matching = lc_matching[0]
    lc_data = files.load_pickle(lc_matching)

    mod = lc_data[chip]['na_data'].split('/')[-1]

    cold_na = files.load_pickle(cold_pkls[np.where(cold_mods==mod)[0][0]])

    plt.figure(figsize=(14,5))
    plt.clf()

    amp = cold_na['carrier_NA']['amp']/cold_na['nuller_NA']['amp']
    plt.plot(cold_na['carrier_NA']['freq'], amp, 'b',
         label = 'Cold NA')
    for i in range(len(lc_data[chip]['optimal_freqs'])):
        label = '%d' % (lc_data[chip]['lc_ind'][i]+1)
        plt.axvline(lc_data[chip]['optimal_freqs'][i], color='gray', 
                linestyle = '--')
        plt.text(lc_data[chip]['optimal_freqs'][i]+1e3, 40, label,
             fontsize=10,rotation=90, color = 'k')

    if mod in warm_mods:
        warm_na=files.load_pickle(warm_pkls[np.where(warm_mods==mod)[0]])
        amp = warm_na['carrier_NA']['amp']/warm_na['nuller_NA']['amp']
        plt.plot(warm_na['carrier_NA']['freq'], amp, 'r', 
             label='Warm NA')

    plt.title(chip+', '+mod.replace('_OUTPUT.pkl',''))
    plt.suptitle('Channels numbered 1-68')
    plt.xlim(1.5e6, 6.3e6)
    plt.legend(loc=1,fontsize=10)

def get_pix_info(wafer = None, pixel= None):
    '''
    Prints the lc channels for the bolometers in the 
    specified pixel and the relevant metahwm entry.
    
    wafer: e.g. 'w172'
    pixel: e.g. 243
    '''
    if wafer==None or pixel==None:
        print(' Please specify a wafer name, e.g. ''w172"'
              'and pixel number, e.g. 151')
    if not isinstance(pixel,int):
        pixel=int(pixel)
        
    metahwm = pd.read_csv(metahwm_path, delimiter = '\t')
    mapping = wafer_bolo_info()

    side=int(mapping.loc[(mapping['pixel']==pixel)]['side'].values[0])
    flex=int(mapping.loc[(mapping['pixel']==pixel)]['flex_cable'].values[0])
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
    print(mod)

def check_for_sq_swaps(metahwm_path):
    '''
    Checks for mislabelled squids in the metahardware map
    
    A common mistake is to mix up the sides of an LC board when
    recording the squid-LC chip matching:
        Sq1<->Sq4, Sq2<->Sq3, Sq5<->Sq6, Sq7<->Sq8
    One way to check for this is to compare the flex cable numbers
    connected to each squid:
        Flex cables 1,2 and 7,8 always plug in to SideB, 
            which is always connected to an odd-numbered squid.
        Flex cables 3,4 and 5,6 always plug in SideA,
            which is always connected to an even-numbered squid.
    '''
    metahwm = pd.read_csv(metahwm_path, delimiter = '\t')
    bad_rows = []
    for i in np.arange(len(metahwm['flex_cable'])):
        if ((metahwm['flex_cable'][i] in ['3,4','5,6'] and 
             metahwm['squid'][i]%2 == 1) or 
            (metahwm['flex_cable'][i] in ['1,2','7,8'] and
             metahwm['squid'][i]%2 == 0)):
            bad_rows.append(i)
    if len(bad_rows)>0:
        print('Checking your meta hardware map for potential squid swaps...')
        print('Found %.0d swapped pairs of squids'%np.ceil(len(bad_rows)/2))
        print('Incorrect lines from metahwm printed below:')
        print([k for k in metahwm.keys()])
        for idx,line in enumerate(bad_rows):
            print([i for i in metahwm.iloc[line]])
            if idx%2 == 1:
                print('\n')
