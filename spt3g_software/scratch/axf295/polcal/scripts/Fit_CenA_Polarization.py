import glob
import os
import argparse as ap
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pylab as pl
import numpy as np
import pickle as pk
import sys
sys.path.append('/home/axf295/2019/code/spt3g_software/scratch/axf295/polcal/python/')
import CenA_Map_Utils as CMU
import General_Utils  as GU
import Polarization_Fitting as PF
import Plotting_Utils as PU
from spt3g import core,std_processing



# Usage: Fit_CenA_Polarization.py -obs [ObsIDs] -o /path/to/outputdir/ -m /path/to/sourcemask/ -s CenA -P plot? -op /path/to/plot_directory/
P = ap.ArgumentParser(description='Fits the polarization angle by comparing single bolo T maps to TQU coadds',
              formatter_class=ap.ArgumentDefaultsHelpFormatter)

P.add_argument('input_files', action='store', nargs='+', default=[], 
               help='Input files; [calframe, singlebolomaps, coaddmaps]')

P.add_argument('-obs','--ObsIDs', action='store', nargs='*', default=[],
                help='CenA Observation IDs')


P.add_argument('-o', '--output_dir', action='store', default='./', 
               help='Output Directory')

P.add_argument('-m', '--mask_dir', action='store',
               default='/home/axf295/2019/code/spt3g_software/polcal/Masks/0p25arcmin/',
               help='source mask location')

P.add_argument('-s', '--source', action='store', default='CenA', help='name of source')

P.add_argument('-P', '--PLOT', action='store', default=0, help='To Plot, or not to Plot')

P.add_argument('-op', '--plot_dir', action='store', default='./tmp_plots/', help='plot directory')

P.add_argument('-r', '--resolution', action='store', default=.25, help='Map Resolution in arcmin')

P.add_argument('-V', '--variable', action='store', default='p', help='Polarization variable to fit (G,p, or both)')

P.add_argument('-n','--Nuclear_Mask',action='store',default = 0, help = 'Do you want to mask the nucleus?')

#P.add_argument('-f','--fft_filt',action='store',default = 0, help = 'Do you want to FFT-filter the bolo maps?')

P.add_argument('-ow','--overwrite',action='store',default = 1, help = 'Do you want to overwrite this observation?')
               

## Set wafers to analyze by adding/removing from this list
wafer_colors_18 = {'W188':'c', 'W174':'b', 'W177':'purple', 'W176':'y', 'W172':'k', 'W180':'orange', 'W181':'brown','W203':'r'}## Removed w2xx #'w203':'r', 'w201':'g',  'w187':'darkblue'

wafer_colors_19 = {'W172':'r', 'W174':'g', 'W176':'c', 'W177':'b', 'W180':'purple', 'W181':'y', 'W188':'k', 'W203':'orange', 'W204':'brown', 'W206':'darkblue'}

## Constants for plotting, etc.
bandmarkers = {'90':'s','150':'o','220':'*'}        
bandcolors  = {'90':'b','150':'g','220':'r'}  
bandfig     = {'90':1,'150':2,'220':3}
bands       = ['90','150','220']
maptypes    = ['T','Q','U']



input_file_list = []

## Load inputs
args     = P.parse_args()
source   = args.source
obsids   = args.ObsIDs
PLOT_DIR = args.plot_dir
PLOT     = args.PLOT
res      = args.resolution
overwrite= args.overwrite
#fftfilt  = args.fft_filt
input_file_list = args.input_files
masknucleus = args.Nuclear_Mask


if len(input_file_list) != 3:
    input_file_list = []
    ## Static data locations
    singlebolomap_loc = '/spt/user/production/calibration/CenA-pixelraster/singlebolomaps/'
    coaddmap_loc      = '/spt/user/production/calibration/CenA-pixelraster/coaddmaps/'
    calframe_loc      = '/spt/user/production/calibration/calframe/CenA-pixelraster/' 
    
print(input_file_list)
    
## HARDCODED!! WARNING!! IF MAPSIZE CHANGES, THIS NEEDS TO CHANGE!    
ps = int(60) ## proper mapsize
proper_size = (ps,ps) ## Size of individual bolo maps  
mapsize = ps//2 ## Half length of map side, in pixels

## Load CenA Mask
mask = {}
nucleusR = {'90':8,'150':6,'220':6}
for band in bands:
    mask[band] = np.loadtxt(args.mask_dir+'%sGHz_Amp_Mask.txt'%band)
    #mask[band] = np.loadtxt(args.mask_dir+'Upper_Lobe_Mask.txt')+ np.loadtxt(args.mask_dir+'Lower_Lobe_Mask.txt')
    masksize = np.shape(mask[band])[0]//2
    ## Mask Nucleus
    if masknucleus:
        #print("Masking the Nucleus!")
        for i in np.arange(0,masksize*2,1):
            for j in np.arange(0,masksize*2,1):
                if np.sqrt((masksize-i)**2+(masksize-j)**2)<nucleusR[band]:
                    mask[band][i][j] *= 0.
    
## Function to get all obsIDs from 2018-2019
if obsids == []:
    obsids = GU.get_cena_obsids(years=[2019])
elif obsids == 'All':
    singlebolomap_loc = '/big_scratch/axf295/2019/CenA/Analyzed_Data/All/'
    coaddmap_loc      = '/big_scratch/axf295/2019/CenA/Analyzed_Data/All/'
    calframe_loc      =  '/spt/user/production/calibration/boloproperties/60000000.g3'
## IF running remotely; need to define filenames 
if len(input_file_list)>0:
    ## if a list was passed, it will stuff that list in another list
    ## this gives us the list we want
    if len(input_file_list)==1:
        input_file_list = input_file_list[0] 
    calframe_loc      = input_file_list[0]
    singlebolomap_loc = input_file_list[1]
    coaddmap_loc      = input_file_list[2]
    

## Cut on low SN
## SN_cut[Nucleus,Lobes]
## Nuclear cuts are just there to cut out problem bolos.
sn_cuts = {}
sn_cuts['90'] = [10.,1.]
sn_cuts['150'] = [10.,1.]
sn_cuts['220'] = [2.,1.]

## By-eye cuts looking at noise/pixel in Noise_Mask Region for good observations
## set stupidly high because S/N cuts usually do the job
noise_cuts = {'90':100, '150':100, '220':200}    
    
## Arbitrary cut based on CenA Brightness, for filtering out a noise spike
peakbrightnesscut = 300
numgoodboloscut   = 3000

failed_obs = {}
for ob in sorted(obsids):
    output_dir = args.output_dir+'%s/'%ob
    try:
        os.mkdir(output_dir)
    except FileExistsError:
        pass
    
    if not overwrite:
        if os.path.exists(output_dir+'%s_Polarziation_Fit_logfile.txt'%(ob)):
            print('skipping %s'%ob)
            failed_obs[ob] = 'Skipped Because Not Overwriting.'
            continue
    print('Analyzing %s'%ob)
    ## Start Logger
    log = open(output_dir+'%s_Polarziation_Fit_logfile.txt'%(ob),'w')
    log.write('Running PolAngle Fitting Script on observation %s\n'%ob)
    for key in vars(args):
        log.write(key+' : '+str(vars(args)[key])+'\n')
    for key in sn_cuts:
        log.write('sn_cut['+key+'] : '+str(sn_cuts[key])+'\n')
        log.write('noise_cuts['+key+'] : '+str(noise_cuts[key])+'\n')
    log.write('numgoodboloscut : '+str(numgoodboloscut)+'\n')
    log.write('peakbrightnesscut : '+str(peakbrightnesscut)+'\n')
    log.write('\n ----------------------------\n')
    
    ## Choose wafers
    if ob == 'All':
        usecoupling = True
        wafers = sorted(wafer_colors_19)
    elif int(ob)>60000000:
        usecoupling = True
        wafers = sorted(wafer_colors_19)
    else:
        usecoupling = False
        wafers = sorted(wafer_colors_18)

    log.write('Using Wafers :' +str(wafers)+'\n')
    
    
    if ob != 'All':
        ## assume this for now... it's accurate enough for 2019 data.
        numframes = 14000
        num_bolos = numframes-3
        
        if len(input_file_list) == 0:
            try:
                ind_map_data = core.G3File(singlebolomap_loc+'%s.g3'%(ob))
            except Exception:
                failed_obs.update({ob:'Single Bolo map file does not exist.'})
                continue
        else:
            try:
                ind_map_data = core.G3File(singlebolomap_loc)
            except Exception:
                failed_obs.update({ob:'Single Bolo map file does not exist.'})
                continue
        
        ## Load bolometer properties and
        ##  autoprocessed TQU maps ; per band coadds only
        ## Can uncomment this to use the online
        if len(input_file_list) == 0:
            cal_file    = calframe_loc+'%s.g3'%ob
            auxfile_loc = output_dir
            coaddfile   = coaddmap_loc+'%s.g3'%ob
        else:
            cal_file    = calframe_loc
            auxfile_loc = './'
            coaddfile   = coaddmap_loc
            
        bolo2band,pname,AB,nominal_angles,unique_nominal_angles,coupling = GU.load_bolometer_properties(cal_file)
        log.write('%s Maps and BoloProps Loaded!\n'%num_bolos)

        
        ## Load the SN calculated in a previous script
        ## Load the center locations
        try:
            bolo_nuclear_sn = pk.load(open(auxfile_loc +'Source_Signal_Dictionary.pkl','rb'))
            bolo_noise      = pk.load(open(auxfile_loc +'Bolo_Noise_Dictionary.pkl','rb'))
            center_locs     = pk.load(open(auxfile_loc +'CenA_center_location_Dictionary.pkl','rb'))
            log.write('Bolo SN files loaded\n')
        except Exception:
            print('Couldnt load a file')
            print(auxfile_loc)
            continue
        ## Set the center location of the source.    
        xo = int(np.nanmedian(center_locs['x']))
        yo = int(np.nanmedian(center_locs['y']))

        template_tqu_maps    = CMU.load_perband_coaddmaps(coaddfile)#,xo=xo,yo=yo)
    
        log.write('Coadd Maps Loaded!\n')

        crapsignalbolos = []
        noiseybolos     = []
        bolo_maps       = {}
        maxpixval       = {}
        improperly_shaped_maps = []
        prob_zero_or_nansinmap = []
        is_resistor            = []
        is_dark                = []
        background = {}
        map_sn     = {}
        for band in bands:
            maxpixval[band]  = []
            background[band] = []
            map_sn[band]     = []
        ## Loop through frames and analyze individual maps
        ## if not map frame, continue
        ## if problematic/noisy bolo, continue
        while True:
            try:
                fr = ind_map_data.next()
                if 'Wunpol' in fr:
                    b = fr['Id']
                    if b == 'map':
                        continue

                    ## Do not include dark bolos, or resistors in bolo_maps
                    if b not in bolo2band or b not in coupling:
                        continue
                    if bolo2band[b] == -1:
                        is_resistor.append(b)
                        continue
                    if usecoupling:
                        if str(coupling[b]) != 'Optical':
                            is_dark.append(b)
                            continue
                    else:
                        if np.isnan(bolo2band[b]):
                            is_dark.append(b)
                            continue
                    bolo_maps[b] = CMU.CenAMap(fr)
                    bolo_maps[b].get_brightest_pixel()
                    bolo_maps[b].center_map()#xo=xo,yo=yo)
                    '''
                    ## filter map in fourier space to remove high-freq stuff
                    if fftfilt:
                        bolo_maps[b].fft_filter_map()
                    '''
                    bmt = bolo_maps[b].maps['T']
                    
                    
                    if b in bolo_nuclear_sn and b in bolo_noise and ~np.any(np.isnan(bmt)):
                        bolo_maps[b].band   = str(bolo2band[b])
                        bolo_maps[b].nomang = int(nominal_angles[b])
                        bolo_maps[b].wafer  = pname[b].split('_')[0].upper()
                        noise  = bolo_noise[b]
                        ## This would work to calc noise if we don't have it already
                        ## may consider a conditional case.
                        #CMU.calc_source_masked_variance(bolo_maps[b].maps['T'],1.-mask[bolo_maps[b].band])
                        
                        if np.isnan(noise) or noise <= 0. :
                            prob_zero_or_nansinmap.append(b)
                            continue
                        elif not CMU.map_is_square(bmt):
                            improperly_shaped_maps.append(b)
                            continue
                        elif noise > noise_cuts[band] :
                            noiseybolos.append(b)
                            continue

                        bnsn = bolo_nuclear_sn[b][0]/bolo_nuclear_sn[b][1]
                        if bnsn < sn_cuts[bolo_maps[b].band][0]  or np.amax(bmt) > peakbrightnesscut:
                            crapsignalbolos.append(b)
                        else:
                            bolo_maps[b].noise = noise
                            bolo_maps[b].sn    = bnsn
                            map_sn[bolo_maps[b].band].append(bolo_maps[b].sn)
                            background[bolo_maps[b].band].append(noise)
                            maxpixval[bolo_maps[b].band].append(np.amax(bmt))

                    else:    
                        prob_zero_or_nansinmap.append(b)
            except StopIteration:
                break
        log.write('%s dark bolos\n'%len(is_dark))
        log.write('%s not optically coupled \n'%len(is_resistor))
        log.write('%s crap signal bolos\n'%len(crapsignalbolos))
        log.write('%s improperly shaped maps\n'%len(improperly_shaped_maps))
        log.write('%s probably zero or nans in map\n'%len(prob_zero_or_nansinmap))

        mapvar3sigmaupper = {}

        log.write('Calculating Per-band Background Noise mean+-std\n')
        
        pl.figure()
        for band in background:
            #print(background[band])
            cleanedbackground = GU.remove_outliers(background[band])
            mapvar3sigmaupper[band] = np.nanmean(cleanedbackground)+3*np.nanstd(cleanedbackground)
            log.write(band+' : '+str(np.nanmean(cleanedbackground))+' +- '+str(np.nanstd(cleanedbackground))+'\n')
            pl.hist(cleanedbackground,bins=np.linspace(0,mapvar3sigmaupper[band],201),alpha=.5,label=band)

        
        pl.legend()
        pl.title('Map background histogram')
        pl.savefig(output_dir+'%s_MapBackHist.png'%(ob))
        pl.close('all')

        if PLOT:
            pl.figure()
            for band in maxpixval:
                pl.hist(maxpixval[band],bins=np.linspace(0,400,201),alpha=.5,label=band)

            pl.legend()
            pl.title('Map max pixelval histogram')
            pl.savefig(output_dir+'%s_MapMaxPixValHist.png'%(ob))
            pl.close('all')
        goodbolos = []
        ## Build a list of good bolos which pass the cuts
        for b in bolo_maps:
            if  bolo_maps[b].noise > mapvar3sigmaupper[band]:
                noiseybolos.append(b)
                continue
            if b not in crapsignalbolos and b not in improperly_shaped_maps and b not in prob_zero_or_nansinmap and b not in noiseybolos and b not in is_dark and b not in is_resistor:
                    goodbolos.append(b)

        log.write(str(len(noiseybolos)) +' noisy bolos\n')
        log.write(str(len(goodbolos))   +' good bolos\n')
        ## Save good bolo list
        with open(output_dir+'%s_GoodBolos.pkl'%(ob), 'wb') as handle:
                pk.dump(goodbolos, handle, protocol=pk.HIGHEST_PROTOCOL)    

        if len(goodbolos)<numgoodboloscut:
            log.write('obs %s Failed Good Bolos Cut; continuing\n'%ob)
            failed_obs.update({ob:'Failed Good Bolos Cut'})
            log.close()
            continue

        #template_tqu_maps = PF.normalize_TQU_template_maps(template_tqu_maps, perbandcoadds,savedir=output_dir,obs=ob, debug = True)

    else:
        bolomap_file  = output_dir+'All_Obs_Single_Bolo_IVweightedCoadds.pkl'
        coaddmap_file = output_dir+'All_Obs_CoaddTQU_IVweightedCoadds.pkl'
        bolo_maps = pk.load(open(bolomap_file,'rb'))
        template_tqu_maps = pk.load(open(coaddmap_file,'rb'))
        ## set arbitrary cut that bolos should have more than 6 observations
        ## in which they are included in the coadd. -- gold standard would be 8
        ## as of 6/16/2019
        goodbolos = []
        for b in bolo_maps:
            if bolo_maps[b].numobs > 6:
                goodbolos.append(b)
            
            
    log.write('Making nominal angle coadd maps\n')
    ## Make nominal angle coadds (split per-wafer and per-band)
    nomang_coadds,numbolos_per_nomang = CMU.make_nominalangle_coadds(bolo_maps,wafers,goodbolos = goodbolos)
    
    
    ## Make per-band T coadds from individual maps
    log.write('Making Per-Band Coadds')
    perbandcoadds,numdetsperband = CMU.make_perband_coadds(bolo_maps,goodbolos = goodbolos)
    PU.plot_perband_coadds(perbandcoadds,numperband = numdetsperband,obs=ob,source=source,savedir=output_dir)
    
    log.write('Saving Nominal Angle Coadd Maps\n')
    with open(output_dir+'%s_NomAngCoaddMaps.pkl'%(ob), 'wb') as handle:
        pk.dump(nomang_coadds, handle, protocol=pk.HIGHEST_PROTOCOL)    
    
    
    
   
    
    
    
    log.write('Calculating Nominal-Angle Best-Fit Angles\n')
    
    nomangfitdata = PF.fit_nominal_angle_coadd_polarization(nomang_coadds,perbandcoadds,template_tqu_maps,
                                                            mask,numbolos_per_nomang,log,
                                                            mapsize=mapsize,variable=args.variable,
                                                            plot_loc=PLOT_DIR,PLOT=PLOT)
    
    '''
    ## get a rotation angle which minimizes chi squared fit
    ## apply that rotation to the nominal angle fits, and store 
    ## that angle to rotate individual map fits.
    rot_ang = PF.get_optimal_nomang_rotation(nomangfitdata['all'])
    nomangfitdata = PF.rotate_nomang_fits(nomangfitdata,rot_ang)
    log.write('Applied a Rotation of %s deg, storing for use in single bolo maps \n'%rot_ang)
    '''
    log.write('Saving Nominal Angle Coadd Fits\n')
    with open(output_dir+'%s_NomAngCoadd_BestFitParams.pkl'%(ob), 'wb') as handle:
        pk.dump(nomangfitdata, handle, protocol=pk.HIGHEST_PROTOCOL)    
    
    
    log.write('Calculating Individual Bolometer Best-Fit Angles\n')
    fitdata = PF.fit_individual_bolometer_polarization(bolo_maps,template_tqu_maps,mask,log=log,
                                                        goodbolos=goodbolos,mapsize=mapsize,
                                                        variable=args.variable,
                                                        plot_loc=PLOT_DIR,PLOT=PLOT)

    print('Saving Individual Bolometer Best-Fit Angles')
    with open(output_dir+'%s_Individual_Bolometer_FitParams.pkl'%(ob), 'wb') as handle:
                pk.dump(fitdata, handle, protocol=pk.HIGHEST_PROTOCOL) 
    
print('The Following Observations failed due to too few good bolos:\n')
for obs in failed_obs:
    print(obs)
