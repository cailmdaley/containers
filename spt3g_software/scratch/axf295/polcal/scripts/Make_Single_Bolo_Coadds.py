from spt3g import core
import sys
sys.path.append('/home/axf295/2019/code/spt3g_software/scratch/axf295/polcal/python/')
import CenA_Map_Utils as CMU
import General_Utils  as GU
import os
import pickle as pk
import numpy as np


## Static data locations
Pointing = 'Offline'
singlebolomap_loc = '/spt/user/axf295/Condor_Return/CenA/%s_Pointing/'%Pointing #'/spt/user/production/calibration/CenA-pixelraster/singlebolomaps/'
coaddmap_loc      = '/spt/user/axf295/Condor_Return/CenA/%s_Pointing/'%Pointing #'/spt/user/production/calibration/CenA-pixelraster/coaddmaps/'
## line 91 edited to use new maps

calframe_loc = '/spt/user/production/calibration/boloproperties/60000000.g3'

SN_dir = '/spt/user/axf295/Condor_Return/CenA/%s_Pointing/'%Pointing

analyzed_data_dir = SN_dir+'both_fit/'

mask_dir = '/home/axf295/2019/code/spt3g_software/scratch/axf295/polcal/Masks/0p25arcmin/'


## Define output directory and make it if necessary
output_dir = analyzed_data_dir+'All/'
if not os.path.exists(output_dir):
    GU.make_directory(output_dir)

analyzed_obs_file = output_dir+'Analyzed_Obs.txt'
analyzed_obs = np.array([])
if os.path.exists(analyzed_obs_file):
    analyzed_obs = np.loadtxt(analyzed_obs_file)
else:
    pass
    
log = open(output_dir+'logfile.txt','w')

## Load all observations
obsids = GU.get_cena_obsids(years=[2019])
failed_obs = {}

## Get latest BoloProperties
boloprops = {}    
for frame in core.G3File(calframe_loc):
    for bolo in frame['BolometerProperties'].keys():
        boloprops[bolo] =  frame['BolometerProperties'][bolo]

## Make Dictionaries for IV-weighted coadds and their weight sums
coadd_map_coadds = {}
single_bolo_coadd_maps = {}
Noise_Mask = {}
for obs in obsids:
    print(obs)
    log.write('Attempting to load data from %s\n'%obs)
    
    ## Make sure data_dir exists
    data_dir = analyzed_data_dir+'%s/'%obs
    if not os.path.exists(data_dir):
        failed_obs[obs] = 'Data path, %s doesn\'t exist'%data_dir
        continue
        
    ## Load the good bolos (as definined in Fit_CenA_Polarization)
    goodbololist = data_dir+'%s_GoodBolos.pkl'%obs
    if os.path.exists(goodbololist):
        good_bolos = pk.load(open(goodbololist,'rb'))
    else:
        print('Couldnt load %s'%(goodbololist))
        failed_obs[obs] = 'Couldnt load %s'%goodbololist
        continue
    
    ## Load the bolo map variance
    mapvarfile = SN_dir+'%s/'%obs + 'Bolo_Noise_Dictionary.pkl'
    if os.path.exists(mapvarfile):
        bolo_map_noise = pk.load(open(mapvarfile,'rb'))
    else:
        print('Couldnt load %s'%(mapvarfile))
        failed_obs[obs] = 'Couldnt load %s'%mapvarfile
        continue         
        
    ## Load the SN calculated in a previous script
    ## Load the center locations
    centerlocfile = SN_dir+'%s/'%obs +'CenA_center_location_Dictionary.pkl'
    try:
        center_locs = pk.load(open(centerlocfile,'rb'))
    except Exception:
        print('Couldnt load %s'%(centerlocfile))
        failed_obs[obs] = 'Couldnt load %s'%(centerlocfile)
        continue
    ## Set the center location of the source.    
    xo = int(np.nanmedian(center_locs['x']))
    yo = int(np.nanmedian(center_locs['y']))
    
    print('Center of obs %s'%obs,xo,yo)
    ## Load autoprocessed TQU maps ; per band coadds only
    coaddfile = coaddmap_loc+'%s_CenA_Coadd.g3'%obs
    if os.path.exists(coaddfile):
        template_tqu_maps = CMU.load_perband_coaddmaps(coaddfile)#,xo=xo,yo=yo)
    else:
        print('Couldnt load %s'%coaddfile)
        failed_obs[obs] = 'Couldnt load %s'%(coaddfile)
        continue
   
            
    log.write('Coadd Maps Loaded!\n')
    
    ## Load up individual bolo maps 
    singlebolofile = singlebolomap_loc+'%s_CenA_SingleBolomaps.g3'%(obs)
    if os.path.exists(singlebolofile):
        ind_map_data = core.G3File(singlebolofile)
    else:
        print('Couldnt load %s'%(singlebolofile))
        failed_obs[obs] = 'Couldnt load %s'%(singlebolofile)
        continue
              
    ## Coadd tqu coadds by IV weighting
    for band in template_tqu_maps:
        if band not in Noise_Mask:
            Noise_Mask[band] = 1.-np.loadtxt(mask_dir+'%sGHz_Amp_Mask.txt'%band)[::-1]
        map_var = CMU.calc_source_masked_variance(template_tqu_maps[band].maps['T'],Noise_Mask = Noise_Mask[band])
        if map_var == -1:
            continue
        if band not in coadd_map_coadds:
            coadd_map_coadds[band] = template_tqu_maps[band]
            for mt in ['T','Q','U']:
                coadd_map_coadds[band].maps[mt]/=np.copy(map_var)
            coadd_map_coadds[band].noise  = np.copy(1./map_var)
            coadd_map_coadds[band].numobs = 1

        else:
            for mt in ['T','Q','U']:
                coadd_map_coadds[band].maps[mt]+=np.copy(template_tqu_maps[band].maps[mt]/map_var)
            coadd_map_coadds[band].noise  += np.copy(1./map_var)
            coadd_map_coadds[band].numobs += 1
    ## Loop through frames and analyze individual maps
    ## if not map frame, continue
    ## if bolo not in goodbolos, continue
    log.write('Coadding Individual Bolo Maps\n')
    while True:
        try:
            fr = ind_map_data.next()
            if 'Wunpol' in fr:
                b = fr['Id']
                if b == 'map':
                    continue
                
                ## Do not include dark bolos, or resistors,or noisy bolos in coadds
                if b not in good_bolos:
                    continue
                if b not in boloprops:
                    continue
                
                band = str(int(boloprops[b].band/core.G3Units.GHz))
                if band == '-1':
                    continue
                    
                if bolo_map_noise[b] == -1:
                    continue
                bolo_map = CMU.CenAMap(fr)
                bolo_map.get_brightest_pixel()
                bolo_map.center_map()#xo=xo,yo=yo)
                
                wafer  = boloprops[b].physical_name.split('_')[0].upper()
                nomang = np.floor(np.rad2deg(boloprops[b].pol_angle))
                bolo_map.band   = band
                bolo_map.wafer  = wafer
                bolo_map.nomang = nomang
                #mapvar = bolo_map_noise[b] **2
                mapvar = CMU.calc_source_masked_variance(bolo_map.maps['T'],Noise_Mask = Noise_Mask[band])
                
                bs = np.shape(bolo_map.maps['T'])
                if bs[0] != bs[1]:
                    continue
                
                ## coadd single bolo maps using IV weighting
                if b not in single_bolo_coadd_maps:
                    single_bolo_coadd_maps[b] = bolo_map
                    single_bolo_coadd_maps[b].maps['T']/= np.copy(mapvar)
                    single_bolo_coadd_maps[b].noise  = np.copy(1./mapvar)
                    single_bolo_coadd_maps[b].numobs = 1
                    
                else:
                    single_bolo_coadd_maps[b].maps['T'] += np.copy(bolo_map.maps['T']/mapvar)
                    single_bolo_coadd_maps[b].noise  += np.copy(1./mapvar)
                    single_bolo_coadd_maps[b].numobs += 1
                    
                
        except StopIteration:
            break
            
        if obs not in analyzed_obs:
            np.append(analyzed_obs,obs)
            
        np.savetxt(analyzed_obs_file,np.asarray(analyzed_obs))
    log.write('Done obs %s\n'%obs)
    log.write('---------------------------\n')

    
for b in single_bolo_coadd_maps:
    single_bolo_coadd_maps[b].maps['T']/=single_bolo_coadd_maps[b].noise
    
for band in coadd_map_coadds:
    for mt in coadd_map_coadds[band].maps:
        coadd_map_coadds[band].maps[mt]/=coadd_map_coadds[band].noise
    
    
with open(output_dir+'All_Obs_Single_Bolo_IVweightedCoadds.pkl', 'wb') as handle:
            pk.dump(single_bolo_coadd_maps, handle, protocol=pk.HIGHEST_PROTOCOL)
        
with open(output_dir+'All_Obs_CoaddTQU_IVweightedCoadds.pkl', 'wb') as handle:
            pk.dump(coadd_map_coadds, handle, protocol=pk.HIGHEST_PROTOCOL)
        
        
log.write('Failed Observations:\n')
failedobs = []
for obs in sorted(failed_obs.keys()):
    if obs not in analyzed_obs:
        failedobs.append(obs)
    log.write(obs+' : ' +failed_obs[obs]+'\n')    
np.append(analyzed_obs,failedobs)
np.savetxt(analyzed_obs_file,np.asarray(analyzed_obs))
log.write('Done!')
log.close()

    
    
