##################################
import sys
import datetime
from glob import glob
import time
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage.filters import gaussian_filter as gf 
from scipy import interpolate
import os
import argparse as ap
from spt3g import core, todfilter, timestreamflagging, gcp, std_processing, mapmaker, dfmux, calibration, xtalk, coordinateutils
from spt3g.coordinateutils.coordsysmodules import FillCoordTransRotations
import pickle as pickle
## To do angle averaging
import CircFun as CF
#######################################


 

    
    
def Make_PerWafer_PerBand_Coadds(input_file_list,Output_dir,obsID,res,source,wafer_lst=[],Overwrite = 0):
    
    outname = Output_dir +'%s_perWafer_coaddmaps_dynamic_noCMfilt.g3'%obsID
    print('Searching for %s'%outname)
    if os.path.exists(outname) and not Overwrite:
        print('%s already exists, skipping because Overwrite is not True'%outname.split('/')[-1])
        return
    
    #maskedpix = Get_Source_Mask(res)
    print('Running PerWafer PerBand Map Maker Script')
    
    start_time = None
    bolos = None
    if not bolos:
        for fname in input_file_list:
            print(fname)
            if (bolos is None) or (start_time is None):
                for frame in core.G3File(fname):
                    
                    if 'CalibratorResponseSN' in frame:
                        bolos  = [k for k in frame['CalibratorResponseSN'].keys()]
                    if 'RawTimestreams_I' in frame:
                        start_time = frame['RawTimestreams_I'].start
                        break
    print(len(bolos),' bolos')
    totnumbolos = len(bolos)
    
    
    starttime = time.time()
    print('Starting at : ',starttime)
    
   
    # Generate map stub
    
    smstub = std_processing.CreateSourceMapStub(
        source.lower(), x_len = 1.5*core.G3Units.deg/res,
        y_len = 1.5*core.G3Units.deg/res, res = res,
        proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
        at_time = start_time)
    
    '''
    psmask = std_processing.CreateSourceMapStub(
    'rcw38', x_len = 1.5*core.G3Units.deg/res,
    y_len = 1.5*core.G3Units.deg/res, res = res,
    proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
    at_time = start_time)
    
    RA  = maskedpix[0]
    Dec = maskedpix[1]
    rad = maskedpix[2]
    mapmaker.fill_point_source_mask_flat_map(RA,Dec,rad,False,False,psmask)
    '''
    pipe = core.G3Pipeline()

    input = input_file_list

    pipe.Add(core.G3Reader, filename=input)
    

    # Combine our various in-progress cal data
    pipe.Add(calibration.build_cal_frames.MergeCalibrationFrames)
    pipe.Add(core.DeduplicateMetadata)
    
    # Cut turnarounds
    pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])

    #pipe.Add(core.Dump)
    
    #only needed for data before feb 2018
    pipe.Add(FillCoordTransRotations,
             transform_store_key = 'OnlineRaDecRotation',
             bs_az_key = 'RawBoresightAz', bs_el_key = 'RawBoresightEl',
             bs_ra_key = 'OnlineBoresightRa', bs_dec_key = 'OnlineBoresightDec',
             do_bad_transform = True)

    pipe.Add(std_processing.flagsegments.FieldFlaggingPreKcmbConversion,ts_key='RawTimestreams_I')
    pipe.Add(std_processing.CalibrateRawTimestreams, units=core.G3TimestreamUnits.Tcmb, output = 'BoloMapTimestreams')
    #pipe.Add(lambda fr: fr.type != core.G3FrameType.Calibration)

    pipe.Add(std_processing.flagsegments.FieldFlaggingPostKcmbConversion,ts_key='BoloMapTimestreams')


    
    pipe.Add(timestreamflagging.RemoveFlagged, 
         input_ts_key = 'BoloMapTimestreams',
         input_flag_key = 'Flags',
         output_ts_key = 'DeflaggedTimestreams')
    
    
    

    pipe.Add(mapmaker.TodFiltering,
         ts_in_key = 'DeflaggedTimestreams', #
         ts_out_key = 'FilteredTs',
         use_dynamic_source_filter = True,
         poly_order = 4,
         #point_source_mask_id = '%sMask'%source,
         point_source_pointing_store_key = 'PixelPointing',
         filter_mask_key = 'FilterMask') #
         
    
    # remove flagged detectors
    pipe.Add(timestreamflagging.RemoveFlagged, 
             input_ts_key = 'FilteredTs',
             input_flag_key = 'Flags',
             output_ts_key = 'BoloMapTimestreamsCMfilt')
    
    #pipe.Add(core.Dump)
    ## No CMFilt for fitting
    #pipe.Add(todfilter.polyutils.CommonModeFilter, in_ts_map_key = 'BoloMapTimestreamsFilt', out_ts_map_key='BoloMapTimestreamsCMfilt',per_wafer=True, per_band=True)
    #pipe.Add(core.Dump)

    pipe.Add(std_processing.weighting.AddMaskedVarWeight, input='BoloMapTimestreamsCMfilt',mask_key='FilterMask', output='TodWeights',store_var = True)#

    #pipe.Add(timestreamflagging.GenerateFlagStats, flag_key='Flags')
    
    pipe.Add(mapmaker.MapInjector,
         map_id = "%sMask"%source,
         maps_lst = [smstub],
         is_stub = True)
    
    pipe.Add(mapmaker.CalculatePointing, map_id='%sMask'%source,
         pointing_store_key = 'PixelPointing',
         ts_map_key = 'DeflaggedTimestreams',
         trans_key='OnlineRaDecRotation')
    
    pipe.Add(core.Delete, keys=['RawTimestreams_I', 'RawTimestreams_Q', 'TimestreamsWatts', 'BoloMapTimestreamsFilt', 'DeflaggedTimestreams', 'FilteredTs', 'BoloMapTimestreams'])
    #pipe.Add(core.Dump)
    
    
    
    pipe.Add(calibration.SplitByWafer, input='BoloMapTimestreamsCMfilt')
    
    if len(wafer_lst) == 0:
        wafer_lst = ['W203', 'W201', 'W188', 'W174', 'W177', 'W176', 'W172', 'W180', 'W181', 'W187']
    # Kick off maps
    for wafer in wafer_lst:
        pipe.Add(mapmaker.MapInjector, map_id = wafer,
                 maps_lst = [smstub], is_stub = True, 
                 make_polarized = True, do_weight = True)

        pipe.Add(mapmaker.BinMap, map_id=wafer,
                 ts_map_key='BoloMapTimestreamsCMfilt'+wafer,
                 pointing_store_key='PixelPointing',
                 timestream_weight_key = 'TodWeights',trans_key = 'OnlineRaDecRotation')
    
        pipe.Add(calibration.SplitByBand, input='BoloMapTimestreamsCMfilt'+wafer,
                 output_root='DfTS_'+wafer+'_')
    
        for band in ['_90GHz', '_150GHz', '_220GHz']:
            pipe.Add(mapmaker.MapInjector, map_id = wafer+band,
                     maps_lst = [smstub], is_stub = True, 
                     make_polarized = True, do_weight = True)

            pipe.Add(mapmaker.BinMap, map_id=wafer+band,
                     ts_map_key='DfTS_'+wafer+band,
                     pointing_store_key='PixelPointing',
                     timestream_weight_key = 'TodWeights',
                     trans_key = 'OnlineRaDecRotation')

    #pipe.Add(core.Dump)
    pipe.Add(mapmaker.MapInjector, map_id = 'SumMap',
             maps_lst = [smstub], is_stub = True, 
             make_polarized = True, do_weight = True)

    pipe.Add(mapmaker.BinMap, map_id='SumMap',
             ts_map_key='BoloMapTimestreamsCMfilt',
             pointing_store_key='PixelPointing',
             timestream_weight_key = 'TodWeights',
             trans_key = 'OnlineRaDecRotation')
    
    pipe.Add(lambda fr: fr.type == core.G3FrameType.Map) # Drop TOD
    #pipe.Add(core.Dump)
    
    pipe.Add(core.G3Writer, filename=outname)
    pipe.Run()
    print('Total Time to Generate Maps: %i seconds'%(time.time()-starttime))
    return


def make_individual_bolo_maps(input_file_list,Output_dir,obsID,res,source,Overwrite=0):
    
    outname = Output_dir +'%s_individual_bolo_maps'%obsID
    #maskedpix = Get_Source_Mask(res)
    
    
    if os.path.exists(Output_dir+'%s_Centered_Individual_Maps.pkl'%obsID) and not Overwrite:
        print('Centered Maps Already Exist, and Overwrite is False. Exiting Individual Bolo Mapmaker')
        return
    
    indbolomapnames = glob(outname+'*.g3')
    for f in indbolomapnames:
        if os.path.exists(f) and not Overwrite:
            print('%s already exists, skipping because Overwrite is not True'%f.split('/')[-1])
            return
    
    print('Running Individual Map Maker Script')

    start_time = None
    bolos = None
    if not bolos:
        for fname in input_file_list:
            print(fname)
            if (bolos is None) or (start_time is None):
                for frame in core.G3File(fname):
                    if 'CalibratorResponseSN' in frame:
                        bolos_all  = [k for k in frame['CalibratorResponseSN'].keys() if frame['CalibratorResponseSN'][k]>=10.]
                    if 'RawTimestreams_I' in frame:
                        start_time = frame['RawTimestreams_I'].start
                        break
    print(len(bolos_all))
    totnumbolos = len(bolos_all)
    bolonums = [0]
    start = 0
    while start < totnumbolos:
        start+=7500
        bolonums.append(start)
    bolonums.append(-1)
    
    starttime = time.time()
    print('Starting at : ',starttime)
    '''
    RA  = maskedpix[0]
    Dec = maskedpix[1]
    rad = maskedpix[2]
    '''
    # Generate map stub
    smstub = std_processing.CreateSourceMapStub(
        source.lower(), x_len = 3.*core.G3Units.deg/res,
        y_len = 3.*core.G3Units.deg/res, res = res,
        proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
        at_time = start_time)
    
    '''
    psmask = std_processing.CreateSourceMapStub(
        'rcw38', x_len = 1.5*core.G3Units.deg/res,
        y_len = 1.5*core.G3Units.deg/res, res = res,
        proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
        at_time = start_time)

    mapmaker.fill_point_source_mask_flat_map(RA,Dec,rad,False,False,psmask)
    '''
    print(bolonums)
    for i in range(len(bolonums)-1):
        print('Running on bolos %s to %s'%(bolonums[i],bolonums[i+1]))
        bolos = bolos_all[bolonums[i]:bolonums[i+1]]
        pipe = core.G3Pipeline()

        input = input_file_list
        print(input)
        pipe.Add(core.G3Reader, filename=input)
        #pipe.Add(core.DeduplicateMetadata)
        
        # Combine our various in-progress cal data
        pipe.Add(calibration.build_cal_frames.MergeCalibrationFrames)
        # Cut turnarounds
        pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])
        pipe.Add(todfilter.util.CutTimestreamsWithoutProperties, input='BoloMapTimestreams')
        
        #only needed for data before feb 2018
        pipe.Add(FillCoordTransRotations,
                 transform_store_key = 'OnlineRaDecRotation',
                 bs_az_key = 'RawBoresightAz', bs_el_key = 'RawBoresightEl',
                 bs_ra_key = 'OnlineBoresightRa', bs_dec_key = 'OnlineBoresightDec',
                 do_bad_transform = False)
        
        
        
        pipe.Add(std_processing.flagsegments.FieldFlaggingPreKcmbConversion,ts_key='RawTimestreams_I')
        pipe.Add(std_processing.CalibrateRawTimestreams, units=core.G3TimestreamUnits.Tcmb, output = 'BoloMapTimestreams')
        #pipe.Add(lambda fr: fr.type != core.G3FrameType.Calibration)

        #pipe.Add(std_processing.flagsegments.FieldFlaggingPostKcmbConversion,ts_key='BoloMapTimestreams')
        
        
        
        pipe.Add(mapmaker.TodFiltering,
             ts_in_key = 'BoloMapTimestreams',
             ts_out_key = 'FilteredTs',
             use_dynamic_source_filter = True,
             poly_order = 4 ,
             #point_source_mask_id = '%sMask',%source
             point_source_pointing_store_key = 'PixelPointing',
             filter_mask_key = 'FilterMask')
        
        
        pipe.Add(mapmaker.MapInjector,
                 map_id = '%sMap'%source,
                 maps_lst = [smstub],
                 is_stub = False)
        
        pipe.Add(mapmaker.CalculatePointing, map_id='%sMap'%source,
             pointing_store_key = 'PixelPointing',
             trans_key='OnlineRaDecRotation',
             ts_map_key = 'FilteredTs')
        
        
        pipe.Add(std_processing.weighting.AddMaskedVarWeight, input='FilteredTs',mask_key='FilterMask', output='TodWeights',store_var = True)#

        #pipe.Add(std_processing.weighting.AddMaskedVarWeight, input='FilteredTs', output='TodWeights')
        
        # do some very nice flagging
        pipe.Add(std_processing.flagsegments.FlagNonResponsive, flag_key = 'Flags')
        pipe.Add(timestreamflagging.FlagNaNs, ts_key='FilteredTs')
        pipe.Add(timestreamflagging.FlagBadHousekeeping,ts_key = 'FilteredTs')
        pipe.Add(timestreamflagging.noiseflagging.FlagUnphysicallyLowVariance, ts_key = 'FilteredTs')
        
        
        # remove flagged detectors
        pipe.Add(timestreamflagging.RemoveFlagged, 
                 input_ts_key = 'FilteredTs',
                 input_flag_key = 'Flags',
                 output_ts_key = 'DeflaggedTimestreams')
        pipe.Add(todfilter.polyutils.CommonModeFilter, in_ts_map_key = 'DeflaggedTimestreams', out_ts_map_key='CMfiltTs',per_wafer=True,per_band=True)

        
        
        #pipe.Add(core.Dump)
        pipe.Add(mapmaker.MapInjector,
                 map_id = '%sMap'%source,
                 maps_lst = [smstub],
                 is_stub = True,
                 make_polarized = False,
                 do_weight = True)

        #creates our weighted histogram
        pipe.Add(mapmaker.BinMap,
                 map_id = '%sMap'%source,
                 ts_map_key = 'CMfiltTs',#'BoloMapTimestreamsFilt',
                 pointing_store_key = 'PixelPointing',
                 timestream_weight_key = 'TodWeights',
                 individual_bolos_to_map = bolos,
                 trans_key = 'OnlineRaDecRotation')
        
        pipe.Add(core.Delete,keys=['DeflaggedTimestreams','FilteredTs','BoloMapTimestreams','BoloMapTimestreamsFilt'])
        
        pipe.Add(lambda fr: fr.type == core.G3FrameType.Map) # Drop TOD
        #pipe.Add(core.Dump)
        
        pipe.Add(core.G3Writer, filename=outname+'_%s.g3'%(i+1))
        pipe.Run()
        
    print('Total Time to Generate Maps: %i seconds'%(time.time()-starttime))
    return





def Center_Maps_on_Brightest_Pixel(input_file_list,obsID,source,res = .25, MapUnits='Power', Output_dir = '',root_folder = '/spt/data/bolodata/downsampled/',Overwrite = 0):
    root_folder+='%s-pixelraster/%s/'%(source,obsID)
    if Output_dir == '':
        Output_dir = '/big_scratch/axf295/2019/%s/Analyzed_Data/'%source
    try:
        cal_file = input_file_list[0]
    except Exception:
        cal_file = '/spt/user/production/calibration/calframe/%s-pixelraster/%s.g3'%(source,obsID)
    
    if not os.path.exists(cal_file):
        print('%s cal file not found; exiting'%obsID)
        return
    
    Bolo_Maps = input_file_list[1:]
    
    print(Bolo_Maps)
    
    if os.path.exists(Output_dir+'%s_Centered_Individual_Maps.pkl'%obsID) and not Overwrite:
        print('Centered Maps Found, and Overwrite is False. Returning from Center_Maps_on_Brightest_Pixel algorithm.')
        return
    
    Source_Centers_Saved = 0
    Save_Individual_Plots = 0
    
    ## Want 30 arcmin maps from 3 deg box
    nominalsize = int(30/res * 3)
    nominal_center = [nominalsize,nominalsize]
    ## Search for individual Center in ~2.8 deg box (trying to avoid edge of obs field)
    centered_image_size =  int(85/res)

    offline_cal = cal_file
    offline_bolo_band = {}
    for frame in core.G3File(offline_cal):
        for bolo in frame['BolometerProperties'].keys():
            if not np.isnan(frame['BolometerProperties'][bolo].band):
                offline_bolo_band[bolo] = int(frame['BolometerProperties'][bolo].band/core.G3Units.GHz)
            else:
                offline_bolo_band[bolo] = frame['BolometerProperties'][bolo].band
    totnumbolos = len(offline_bolo_band.keys())
    print('There are %s bolometers in the offline_cal_file '%totnumbolos)
    
    print('Using offline cal for bolo band mapping')
    bolo_band = offline_bolo_band

    print('Centering bolo_maps')
    
    numbolos = 0
    weight = core.G3Units.mK
    goodbolonum = 0
    Centered_Source_Maps = {}        
    Source_Centers = {}
    plotnums = 0
    bolonum = 0
    percent = 0
    try:
        for i in range(len(Bolo_Maps)):
            datafile = core.G3File(Bolo_Maps[i])
            while True:
                try:
                    frame = datafile.next()
                    #print(frame)
                    if 'T' in frame.keys() and frame['Id']!='bsmap':
                        Id = frame['Id']
                        #print(Id)
                        if 'Wunpol' in frame:
                            w = np.array(frame['Wunpol'].TT)[::-1]
                            if MapUnits == 'Power':
                                ## put into fW
                                weight = -1. *core.G3Units.W / 1e15 #
                            
                            else:
                                ## put into mK
                                weight = 1* core.G3Units.mK

                        if MapUnits == 'Power':
                            outmapsT = np.array(frame['T'])[::-1]/weight
                            
                        elif MapUnits == 'Tcmb':
                            outmapsT = np.array(frame['T'])[::-1]/weight

                        numbolos+=1
                        try:
                            bolo = Id
                            bolonum+=1
                            if int(100*(bolonum/totnumbolos)) > percent:
                                print("Analyzing bolo %s out of %s : %s%% "%(bolonum,totnumbolos,percent))
                                percent+=1
                            band = str(bolo_band[bolo])
                            xlen,ylen = np.shape(outmapsT)
                            xlen//=2
                            ylen//=2
                            sourceT = outmapsT[xlen-centered_image_size:xlen+centered_image_size,ylen-centered_image_size:ylen+centered_image_size]
                            sourceW = w[xlen-centered_image_size:xlen+centered_image_size,ylen-centered_image_size:ylen+centered_image_size]
                            smoothedT = gf(sourceT, 1., order=0,  mode='reflect', cval=0.0, truncate=3.0)
                            maxpixels = np.where(smoothedT== np.amax(np.nan_to_num(smoothedT)))
                            xo = maxpixels[1][0]
                            yo = maxpixels[0][0]
                            #print(xo,yo)
                            dx = int(15/res)
                            dy = int(15/res)
                            if xo >= dx and yo>=dy and not np.all(sourceT == 0.) and not np.all(sourceW==0.): 
                                centered_source = {}
                                sourceT /=sourceW
                                if plotnums<0:
                                    plt.figure(bolo)
                                    plt.title(bolo+' %s GHz'%band)
                                    plt.imshow(sourceT)
                                    plt.colorbar()
                                    plt.scatter(xo,yo,marker='x',s=500,color='w',label=maxpixels)
                                    plt.legend()
                                    plt.savefig('/home/axf295/2019/tmpPlots/%sCenterLoc.png'%(bolo))
                                    plt.close('all')
                                    plotnums+=1
                                #print(np.shape(sourceT))
                                '''
                                if res == .5:
                                    sourceT = sourceT.repeat(2,axis=0).repeat(2,axis=1)
                                    xo*=2
                                    yo*=2
                                    dx*=2
                                    dy*=2
                                #print(np.shape(sourceT))
                                '''
                                moi = np.copy(sourceT[yo-dy:yo+dy,xo-dx:xo+dx])
                                if plotnums<0:
                                    plt.figure(bolo)
                                    plt.title(bolo+' %s GHz'%band)
                                    plt.imshow(sourceT,vmin=0,vmax=10)
                                    plt.colorbar()
                                    plt.scatter(xo,yo,marker='x',s=500,color='w',label=maxpixels)
                                    plt.legend()
                                    plt.savefig('/home/axf295/2019/tmpPlots/%sCentered_Tmap.png'%(bolo))
                                    plt.close('all')
                                
                                
                                
                                if not np.all(np.isfinite(moi)):
                                    try:
                                        centered_source['T'] = Interpolate_Map_NaNs(moi)
                                    except ValueError:
                                        #print('probably no finite values')
                                        centered_source['T'] = np.zeros([2*dx,2*dy])
                                else:
                                    centered_source['T'] = moi
                                
                                if plotnums<0:
                                    plt.figure(bolo)
                                    plt.title(bolo+' %s GHz'%band)
                                    plt.imshow(centered_source['T'])
                                    plt.colorbar()
                                    plt.savefig('/home/axf295/2019/tmpPlots/%sCentered_interped.png'%(bolo))
                                    plt.close('all')
                                    
                                Centered_Source_Maps[bolo] = centered_source
                                Source_Centers[bolo] = (yo+(ylen-centered_image_size),xo+(xlen-centered_image_size))
                                goodbolonum+=1
                                #print(goodbolonum)
                        except (IndexError, KeyError) as e:
                            #print(bolo,e)
                            pass
                except StopIteration:
                    break
                    
                    
    except RuntimeError:
        print('Error on file: %s'%Bolo_Maps[i])
        print('There is either a problem with this file or it is still being created.')
    
    print('There are %s bolos in this map, of which %s were centered'%(numbolos,goodbolonum))
    print('\n ################################################### ')
    print('################################################### \n')
    
    print('%s bolos passed though map centering algorithm'%goodbolonum)

    print('\n Saving Centered Maps')
    
    with open(Output_dir+'%s_Centered_Individual_Maps.pkl'%obsID, 'wb') as handle:
            pickle.dump(Centered_Source_Maps, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(Output_dir+'%s_Center_Locations.pkl'%source, 'wb') as handle:
            pickle.dump(Source_Centers, handle, protocol=pickle.HIGHEST_PROTOCOL)


    print('Done Centering Maps \n')
    print('Done Source Map Centering Algorithm')
    return

def pix2dist(x0,y0,x1,y1,pix2angle):
        pixdist = np.sqrt(abs(x0-x1)**2+abs(y0-y1)**2)
        dist = pixdist*pix2angle
        return dist

def GetFluxAroundPix(xpix,ypix,pix2ang,data,radius=2,xo=0,yo=0,circcolor='r',source=False,debug=False):
    integral_flux = 0.0
    pixinregion = []
    xregion = np.array(np.arange(int(xo)-int(2*radius),int(xo)+int(2*radius),pix2ang)/pix2ang,dtype=int)
    yregion = np.array(np.arange(int(yo)-int(2*radius),int(yo)+int(2*radius),pix2ang)/pix2ang,dtype=int)
    #print xo,yo
    #print(xregion,yregion)
    for i in xregion:
        for j in yregion:
            if i>=0 and j >=0 and i<len(xpix) and j<len(ypix):
                x1 = xpix[i]*pix2ang
                y1 = ypix[j]*pix2ang
                #print x1,y1
                dist = pix2dist(xo,yo,x1,y1,pix2ang)
                if dist<=radius*pix2ang:
                    #print x1,y1
                    #plt.scatter(x1,y1,marker='+',color='y')
                    if source:
                        integral_flux+=data[i][j]
                    else:
                        integral_flux+=data[i][j]**2
                    pixinregion.append(data[i][j])
    if integral_flux == 0.0:
        circcolor='w'
    if not source:
        integral_flux = np.sqrt(integral_flux)
    #print 'Circle going to : ',yo,xo
    if debug:
        plt.figure('Source SN Regions')
        arcmincirc = plt.Circle((yo, xo), radius,facecolor='none', edgecolor=circcolor,lw=2)
        plt.axes().add_artist(arcmincirc)
        plt.annotate('%.2f'%integral_flux, xy=(yo+1, xo+1),color='w',fontsize=12)

    return integral_flux,pixinregion






def calc_signal_to_noise(map_dir,obsID,source,Output_dir = '',search_radius = 1.5,Save_pngs = 0):
    '''
    res is pixel scale in arcmin
    Search radius in arcmin.
    '''
    print('\n################################################### ')
    print('Welcome to the Source-Centered Map S/N Calculating Script.')
    print(' ################################################### \n')
    
    
    if Output_dir=='':
        Output_dir = map_dir
    
    ## as of 5/8/2019 singlebolomaps are CenA-centered
    '''
    if os.path.exists(Output_dir+'%s_Centered_Individual_Maps.pkl'%obsID):
        individual_maps = pickle.load(open(Output_dir+'%s_Centered_Individual_Maps.pkl'%obsID,'rb'))
    else:
        print(Output_dir+'%s_Centered_Individual_Maps.pkl not found; returning.'%obsID)
        return
    '''
    datafile = '/spt/user/production/calibration/CenA-pixelraster/singlebolomaps/%s.g3'%obsID
    if os.path.exists(datafile):
        ind_map_data = core.G3File(datafile)
    else:
        print('%s does not exist; returning.'%datafile)
        return
    
    
    starttime = time.time()
    
    print('Reading in %s'%(datafile))
    
    num_frames = 0
   
    while True:
        try:
            fr = ind_map_data.next()
            num_frames+=1
        except StopIteration:
            break
    dt = time.time()-starttime
    
    print('Read in %s frames in %.2f s'%(num_frames,dt))
    num_bolos = num_frames-3

    left  = 0
    right = -1
    
    bolonum = 0
    
    print('Starting SN Calculations')
    
    startofSNcalc = time.time()
    
    
    analnum = 0
    Source_SN = {}
    Integral_Flux = {}
    ET = 0
    percent = 0
    bolonum = 0
    ind_map_data = core.G3File(datafile)
    while True:
        try:
            fr = ind_map_data.next()
            if 'Wunpol' in fr:
                bolo = fr['Id']
                if bolo == 'map':
                    continue
                mapT = np.asarray(fr['T']/fr['Wunpol'].TT)[::-1]
                res = fr['T'].res/core.G3Units.arcmin
                
                analnum+=1
                SN = -1
                if int(100*(bolonum/num_bolos)) > percent:
                    print("Analyzing bolo %s out of %s : %s%% "%(bolonum,num_bolos,percent))
                    percent+=1
                bolonum+=1
                
                data_im = gf(mapT, .5/res, order=0,  mode='reflect', cval=0.0, truncate=3.0)

                xlen,ylen = np.array(np.shape(data_im))

                if xlen==ylen:
                    pixel_length= res
                    flux_xs = np.arange(5,int((xlen-5)*pixel_length),int(xlen/3.0)*pixel_length)
                    flux_ys = np.arange(5,int((ylen-5)*pixel_length),int(ylen/3.0)*pixel_length)
                    xpix = np.linspace(1,xlen,xlen)
                    ypix = np.linspace(1,ylen,ylen)

                    if Save_pngs:
                        plt.figure('%s SN Regions'%source,figsize=(10,10))
                        plt.imshow(data_im,extent=[0,xlen*pixel_length,  ylen*pixel_length, 0],cmap='bone',interpolation='none')
                        plt.colorbar()
                        plt.grid(color='w')

                    xo,yo = xlen/2*pixel_length,ylen/2*pixel_length

                    dataradius = search_radius 

                    sourceflux,sourcepix = GetFluxAroundPix(xpix,ypix,pixel_length,data_im,radius = dataradius,xo=xo,yo=yo,circcolor='r',source = True,debug=Save_pngs)
                    sourcestd = np.std(sourcepix)

                    noiseflux = []
                    noisestd  = []
                    if Save_pngs:
                        plt.figure('Pixel Raster Plot '+ bolo,figsize=(10,10))
                    for i in range(len(flux_xs)):
                        for j in range(len(flux_ys)):
                            if flux_xs[i]!=flux_ys[j]:
                                flux,pix = GetFluxAroundPix(xpix,ypix,pixel_length,data_im,radius=dataradius, xo=flux_xs[i],yo=flux_ys[j],circcolor='y',debug=Save_pngs)
                                if flux!=0.0:
                                    noiseflux.append(flux)
                                    noisestd.append(np.std(pix))
                                    if Save_pngs:
                                        plt.figure('Pixel Raster Plot '+bolo)
                                        plt.plot(pix,label='Flux around center %s,%s'%(flux_xs[i],flux_ys[j]))

                    if Save_pngs:
                        plt.figure('Pixel Raster Plot '+bolo)            
                        plt.plot(sourcepix,color='k',lw=1,label='Source Flux')
                        plt.legend(fontsize=12)
                        plt.title('Pixel Raster Plot '+bolo,fontsize=20)
                        plt.savefig(Output_dir+bolo+'_SNregionraster.png')
                    if len(noiseflux)>0:
                        meannoise = np.nanmean(noiseflux)
                    else:
                        meannoise='nan'

                    meanstd = np.nanmean(noisestd)

                    if Save_pngs:
                        if meannoise!='nan':
                            print('Source flux and std: %.2f, %.2f'%(sourceflux,sourcestd))
                            print('Mean Noise Region Flux and std: %.2f, %.2f'%(meannoise,meanstd))
                        else:
                            print('Noise is nan')

                    ## Remove ridiculous noise or SN values (set SN to 0)
                    if meannoise == 'nan' or meannoise <= 1e-6 or meannoise > 1e6:
                        SN = 0
                    else:
                        SN = int(np.nan_to_num(sourceflux/meannoise))
                    if SN>1e6 or SN < 0:
                        SN = 0
                    #print('S/N: ', SN)
                    if Save_pngs:
                        plt.figure('Source SN Regions')
                        plt.title('Source '+bolo,fontsize=20)
                        if SN!=0:
                            plt.annotate('Signal: %.3f , Noise: %.3f, S/N: %i '%(sourceflux,meannoise,SN), xy=(0.55, 0.85), xycoords='axes fraction',color='w',fontsize=12)
                        else:
                            plt.annotate('Signal: %.3f , Noise: %s, S/N: %i '%(sourceflux,meannoise,SN), xy=(0.55, 0.85), xycoords='axes fraction',color='w',fontsize=12)

                        plt.savefig(Output_dir+bolo+'_%sSNregion.png'%source)
                        plt.close('all')

                    Source_SN[bolo] = SN
                    Integral_Flux[bolo] = sourceflux
            et = time.time()-startofSNcalc
            if et > 5.+ET:
                ET = et
                print('\n ------------------------------------------------------')
                print('Elapsed time in SN calc: %i s.'%(ET))
                print('Average time per bolo: %.3f s.'%(ET/float(analnum)))
                print('Estimated time remaining: %.3f s.'%(ET/float(analnum)*(num_bolos-analnum)))
                print(' ------------------------------------------------------\n')

        except StopIteration:
            break
    with open(Output_dir+'%s_Core_SN_Dictionary.pkl'%source, 'wb') as handle:
        pickle.dump(Source_SN, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print('Plotting the Nucleus SN histograms')    
    offline_cal = '/spt/data/bolodata/downsampled/%s-pixelraster/%s/offline_calibration.g3'%(source,obsID)
    if not os.path.exists(offline_cal):
        print('Using offline cal file from most recent bolo props')
        offline_cal = '/spt/user/production/calibration/boloproperties/60000000.g3'
    
    bolo_band = {}
    for frame in core.G3File(offline_cal):
        for bolo in frame['BolometerProperties'].keys():
            if not np.isnan(frame['BolometerProperties'][bolo].band):
                bolo_band[bolo] = str(int(frame['BolometerProperties'][bolo].band/core.G3Units.GHz))
            else:
                bolo_band[bolo] = -1
    
    print('plotting')

    bandcolors = {}
    bandcolors['90']  = 'b'
    bandcolors['150'] = 'g'
    bandcolors['220'] = 'r'
    
        
    SourceSN = {}
    SourceSN['90']  = []
    SourceSN['150'] = []
    SourceSN['220'] = []
    
    PerBandFlux = {}
    PerBandFlux['90']  = []
    PerBandFlux['150'] = []
    PerBandFlux['220'] = []

    for bolo in Source_SN:
        if Source_SN[bolo]>0. and bolo_band[bolo]!=-1:
            SourceSN[bolo_band[bolo]].append(Source_SN[bolo])
            if Integral_Flux[bolo] < 1e6 and Source_SN[bolo]>10.:
                PerBandFlux[bolo_band[bolo]].append(Integral_Flux[bolo])

    print('bolos SN and IF added to dict')
    
    
    plt.figure('%s Source SN'%obsID)
    plt.title('%s Source SN Histogram'%obsID)
    SN_vmax = {}
    for band in SourceSN:
        medianSN = np.median(SourceSN[band])
        stdSN = np.std(SourceSN[band])
        SN_vmax[band] = medianSN+2*stdSN
        BINS = np.linspace(medianSN-3*stdSN,medianSN+3*stdSN,101)
        plt.hist(SourceSN[band],bins=BINS,color=bandcolors[band],alpha=.5,label=band+' GHz')
    plt.xlabel('%s Nucleus SN'%source,fontsize=16)
    plt.legend()
    #plt.show()
    plt.savefig(Output_dir+'%s_%sNucleus_SN_Hist_byband.png'%(source,obsID))
    plt.close()

    
    plt.figure('%s Source Integral Flux'%obsID)
    plt.title('%s Source Integral Flux Histogram'%obsID)
    IF_vmax = {}
    for band in PerBandFlux:
        medianIF = np.nanmedian(PerBandFlux[band])
        stdIF = np.nanstd(PerBandFlux[band])
        IF_vmax[band] = medianIF+2*stdIF
        BINS = np.linspace(medianIF-3*stdIF,medianIF+3*stdIF,101)
        plt.hist(PerBandFlux[band],bins=BINS,color=bandcolors[band],alpha=.5,label=band+' GHz: %.2f +- %.2f'%(medianIF,stdIF))
    plt.xlabel('%s Nuclear Flux'%source,fontsize=16)
    plt.legend()
    #plt.show()
    plt.savefig(Output_dir+'%s_%sNucleus_Flux_Hist_byband.png'%(obsID,source))
    plt.close()
    
    
    with open(Output_dir+'%s_PerBand_Nuclear_Flux.pkl'%source, 'wb') as handle:
        pickle.dump(PerBandFlux, handle, protocol=pickle.HIGHEST_PROTOCOL)

    
    #nom = list(core.G3File(nominal_online_file))[0]
    off = core.G3File(offline_cal).next()
    good_pixels = dict()
    x_offsets  = []
    y_offsets  = []
    lost_bolos = []
    for band in SourceSN:
        plt.figure('Source SN, %s GHz'%band,figsize=(10,10))
        numbolos = 0
        VMAX = SN_vmax[band]

        for bolo in off['BolometerProperties'].keys():
            if np.isfinite(off['BolometerProperties'][bolo].pol_angle):
                x = off['BolometerProperties'][bolo].x_offset
                y = off['BolometerProperties'][bolo].y_offset
                NA =int(np.rad2deg(off['BolometerProperties'][bolo].pol_angle%np.pi))
            
                if bolo in bolo_band:
                    if bolo_band[bolo] == band and bolo in Source_SN:
                        bolomark = (2,0,NA)
                        plt.scatter(x,y,marker=bolomark,c=Source_SN[bolo],vmin=0,vmax=VMAX,cmap='plasma',edgecolor='face')
                        numbolos+=1
        print(numbolos, '%s GHz bolos plotted'%band)
        if numbolos!=0:
            plt.colorbar(label='Source SN')
        plt.title('Obs %s, %s GHz'%(obsID,band))
        plt.savefig(Output_dir+'%s_%sNucleus_SN_FP_map_%sGHz.png'%(obsID,source,band))
        plt.close()
    print('\n ################################################### ')
    print('Congrats! You are finished calculating source SN. Elapsed time: %i s.'%(time.time()-starttime))
    print(' ################################################### \n' )
    
    return







## Bolometer Polarization Angle Fitting Functions





def Load_PerWafer_CoaddMaps(coadd_file,wafers,factor = 1.):
    ## Read in Per-wafer Per-band Coadds
    ## 
    wafer_coadds = {}
<<<<<<< HEAD
    try:
        for frame in core.G3File(coadd_file):
=======
    fpath = coadd_file #map_dir+'%s/%s_perWafer_coaddmaps.g3'%(ObsID,ObsID)
    try:
        for frame in core.G3File(fpath):
>>>>>>> 884205b41e8b8e6d95b86f88d73a568bba4a10f3
                waferband = frame['Id'].split('_')
                if len(waferband)==2:
                    wafer = waferband[0]
                    band  = waferband[1]
                    if wafer.upper() in wafers or wafer.lower() in wafers:
                        if wafer not in wafer_coadds:
                            wafer_coadds[wafer] = {}
                            wafer_coadds[wafer]['90GHz']  = {}
                            wafer_coadds[wafer]['150GHz'] = {}
                            wafer_coadds[wafer]['220GHz'] = {}
                        wafer_coadds[wafer][band]['T'] = np.nan_to_num((np.asarray(frame['T'])/frame['Wpol'].TT)[::-1]*factor)
                        wafer_coadds[wafer][band]['Q'] = np.nan_to_num((np.asarray(frame['Q'])/frame['Wpol'].QQ)[::-1]*factor)
                        wafer_coadds[wafer][band]['U'] = np.nan_to_num((np.asarray(frame['U'])/frame['Wpol'].UU)[::-1]*factor)
    except RuntimeError:
            print(wafer_coadds)
    return wafer_coadds

def Plot_PerWafer_Coadds(wafer_coadds,obsID,source, savedir='/.'):
    waferfig = {}
    i = 1
    for w in sorted(wafer_coadds.keys()):
        waferfig[w.upper()] = i
        i+=1
    bandfig = {'90GHz':1,'150GHz':2,'220GHz':3}
    vlims = {}
    vlims['T'] = [-2,2]
    vlims['Q'] = [-1,1]
    vlims['U'] = [-1,1]
    for mt in ['T','Q','U']:
        plt.figure(figsize=(35,50))
    
        for wafer in wafer_coadds:
            for band in wafer_coadds[wafer]:
                smoothed = wafer_coadds[wafer][band][mt]#gf(wafer_coadds[wafer][band][mt], 1., order=0,  mode='reflect', cval=0.0, truncate=3.0)
                dx = np.shape(smoothed)[0]//2
                imsize = 30
                zoom = smoothed[dx-imsize:dx+imsize,dx-imsize:dx+imsize]
                plt.subplot(10,3,3*(waferfig[wafer]-1)+bandfig[band])
                plt.imshow(zoom)#,vmin=vlims[mt][0],vmax=vlims[mt][1])
                plt.colorbar(fraction=.046,pad=.04)
                plt.title('%s, %s '%(band,wafer))
        plt.suptitle('%s %s Coadds, Perwafer,Perband, Obs %s'%(source,mt,obsID))
        plt.savefig(savedir+'%s_PerWafer_PerBand_%sCoaddmaps_%s.png'%(source,mt,obsID))
        plt.close('all')
    return

<<<<<<< HEAD

=======
def Plot_PerNomang_Coadds(nomang_coadds,obsID,source,savedir='/.'):
    waferfig = {}
    i = 1
    for w in sorted(nomang_coadds.keys()):
        waferfig[w.upper()] = i
        i+=1
    bandfig = {90:1,150:2,220:3}
    
    
    for wafer in nomang_coadds:
        plt.figure(figsize=(20,15))
        angnum = 1
        for ang in nomang_coadds[wafer]:
            for band in nomang_coadds[wafer][ang]:
                plt.subplot(4,3,3*(angnum-1)+bandfig[band])
                
                smoothed = nomang_coadds[wafer][ang][band] #gf(nomang_coadds[wafer][ang][band], 2., order=0,  mode='reflect', cval=0.0, truncate=3.0)
                plt.imshow(smoothed)
                plt.colorbar(fraction=.046,pad=.04)
                plt.title('%sGHz, %s %s deg'%(band,wafer,ang))
            angnum+=1
        plt.suptitle('%s %s Nomang Coadds Obs %s'%(source,wafer,obsID))
        plt.savefig(savedir+'%s_%s_NomAng_Coaddmaps_%s.png'%(source,wafer,obsID))
        plt.close('all')
    return
>>>>>>> 884205b41e8b8e6d95b86f88d73a568bba4a10f3




def Center_PerBandMaps(bandcoadds,Noise_Mask,mapsize = 30,Output_Dir = './'):
    ## Center the band summaps
    dx = mapsize
    weightperband = {}
    centeredMaps = {}
    for band in bandcoadds:
        centeredMaps[band] = {}
        
        ## Get Region of interest -- called 'map of interest' or moi
        s = np.shape(bandcoadds[band]['T'])
        moi_T = bandcoadds[band]['T'][(s[0]//2-dx):(s[0]//2+dx),(s[1]//2-dx):(s[1]//2+dx)]
        smoothed_T = gf(moi_T, 1., order=0,  mode='reflect', cval=0.0, truncate=3.0)
        maxT = np.where(smoothed_T == np.amax(smoothed_T))
        offsetx = dx-maxT[0][0]
        offsety = dx-maxT[1][0]
        #print(offsetx,offsety)
        ## Get noise region variance
        mapvar = np.var(bandcoadds[band]['T'][(s[0]//2-dx-offsetx):(s[0]//2+dx-offsetx),(s[1]//2-dx-offsety):(s[1]//2+dx-offsety)][np.nonzero(Noise_Mask)])
        ## Center moi on brightest pixel
        for mt in ['T','Q','U']:
            moi = bandcoadds[band][mt][(s[0]//2-dx-offsetx):(s[0]//2+dx-offsetx),(s[1]//2-dx-offsety):(s[1]//2+dx-offsety)]

            centeredMaps[band][mt] = np.copy(moi)

    
    with open(Output_Dir+'Centered_Coadds.pkl', 'wb') as handle:
            pickle.dump(centeredMaps, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    return centeredMaps



def Center_PerWaferMaps(wafer_coadds,Noise_Mask = [],mapsize = 30,IVweighting=True,Output_Dir = './'):
    ## Create the Sum Map from wafer coadds.
    ## Inverse variance weight using the noise mask if IVweighting
    
    dx = mapsize
    weightperband = {}
    Summap = {} 
    centeredMaps = {}
    for wafer in wafer_coadds:
        centeredMaps[wafer] = {}
        for band in wafer_coadds[wafer]:
            centeredMaps[wafer][band] = {}
            if band not in Summap:
                Summap[band] = {}
                weightperband[band] = 0.
            
            
            ## Get Region of interest -- called 'map of interest' or moi
            s = np.shape(wafer_coadds[wafer][band]['T'])
            #print(wafer,band)
            moi_T = wafer_coadds[wafer][band]['T'][(s[0]//2-dx):(s[0]//2+dx),(s[1]//2-dx):(s[1]//2+dx)]
            maxT = np.where(moi_T == np.amax(moi_T))
            #print(maxT)
            offsetx = dx-maxT[0][0]
            offsety = dx-maxT[1][0]
            ## Get noise region variance
            
            mapvar = np.var(wafer_coadds[wafer][band]['T'][(s[0]//2-dx-offsetx):(s[0]//2+dx-offsetx),(s[1]//2-dx-offsety):(s[1]//2+dx-offsety)][np.nonzero(Noise_Mask)])
            ## Center moi on brightest pixel
            for mt in ['T','Q','U']:
                moi = wafer_coadds[wafer][band][mt][(s[0]//2-dx-offsetx):(s[0]//2+dx-offsetx),(s[1]//2-dx-offsety):(s[1]//2+dx-offsety)]
                
                centeredMaps[wafer][band][mt] = np.copy(moi)
                
                if mt not in Summap[band]:
                    Summap[band][mt] = np.zeros((2*dx,2*dx))
                if IVweighting:
                    moi/=mapvar
                Summap[band][mt] += np.copy(moi)
            if IVweighting:
                weight = 1./mapvar
            else:
                weight = 1.
            weightperband[band] += weight
            
    for band in Summap:
        for mt in Summap[band]:
            Summap[band][mt] /= weightperband[band]
    
    with open(Output_Dir+'Centered_PerWafer_Coadds.pkl', 'wb') as handle:
            pickle.dump(centeredMaps, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open(Output_Dir+'Summap.pkl', 'wb') as handle:
            pickle.dump(Summap, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    return Summap,centeredMaps

def Calc_Map_Background(bolo_map,Noise_Mask):
    ms = np.shape(Noise_Mask)[0]
    dsx = np.shape(bolo_map)[0]
    dsy = np.shape(bolo_map)[1]
    if dsx > ms or dsy > ms:
        bolo_map = bolo_map[dsx//2-ms//2:dsx//2+ms//2, dsy//2-ms//2:dsy//2+ms//2]
        bolonoisemask = Noise_Mask
    elif dsx < ms or dsy < ms:
        bolonoisemask = Noise_Mask[ms//2-dsx//2:ms//2+dsx//2, ms//2-dsy//2:ms//2+dsy//2]
    else:
        bolonoisemask = Noise_Mask
    ms = np.shape(bolonoisemask)[0]     
    numpix = len(bolo_map[np.nonzero(bolonoisemask)])
    try:
        background = np.sum(abs(bolo_map[np.nonzero(bolonoisemask)]))/numpix
    except IndexError:
        background = -1
    return background




def Make_Variance_Map_From_Centered_Maps(map_location,ob,save_dir=''):
    maps = pk.load(open(map_location,'rb'))
    ## Hopefully the first bolo has a properly shaped map
    varmap = np.zeros(np.shape(maps[sorted(maps.keys())[0]]['T']))
    mapshape = np.shape(varmap)
    
    bolo2band,pname,AB,nominal_angles,unique_nominal_angles = Load_Bolometer_Properties(ob,ABcorr=False)
    
    propslists = {}
    bands = [90,150,220]
    for band in bands:
        propslists[band] = {}
    
    for b in bolo2band:
        band = bolo2band[b]
        if band in propslists:
            if nominal_angles[b] not in propslists[band]:
                propslists[band][nominal_angles[b]] = []
            propslists[band][nominal_angles[b]].append(b)
    perbandperangVar = {}
    for band in propslists:
        perbandperangVar[band] = {}
        for nomang in propslists[band]:
            
            bolos = propslists[band][nomang]
            numbolos = 0
            for b in bolos:
                if b in maps:
                    numbolos+=1
            
            print(band,'GHz, ',nomang,'deg: ',numbolos)
            pixelvalarray = np.zeros((numbolos,len(np.ndarray.flatten(varmap))))
            bolonum = 0
            for b in bolos:
                if b in maps:
                    bm = np.ndarray.flatten(np.asarray(maps[b]['T']))
                    for i in range(len(bm)):
                        pixelvalarray[bolonum,i] = bm[i]
                    bolonum+=1
            flatvar = np.zeros(len(np.ndarray.flatten(varmap)))
            for i in range(len(flatvar)):
                pixdata = CAF.remove_outliers(pixelvalarray[:,i],m=3)
                flatvar[i] = np.nanvar(pixdata)
            perbandperangVar[band][nomang] = np.reshape(flatvar,mapshape)
            
    with open(save_dir+'%s_perBandperAngle_VarianceMaps.pkl', 'wb') as handle:
        pickle.dump(perbandperangVar,handle,protocol=pickle.HIGHEST_PROTOCOL)
        
    return   



def Make_NominalAngle_Coadds(bolo_maps,bolo_map_variance,wafers,pname,bolo2band,nominal_angles,goodbolos = [], noiseybolos=[],crapsignalbolos=[],IVweighting = True,proper_size=(120,120)):
    ## Returns nominal angle coadd maps per-wafer per-band
    ## coadds are variance weighted if IVweighting == True
    bands = ['90','150','220']
    if len(goodbolos)==0:
        print('No Good Bolos, returning')
        return {},{},{}
    fitdata = {}
    if IVweighting:
        IVfactor = 1.
    else:
        IVfactor = 0.
    nomangle_maps = {}
    numanglebolos = {}
    varweights = {}
    wafer_coadds = {}
    wafer_coadds['all'] = {}
    all_wafer_coadds = wafer_coadds['all']
    all_wafer_coadd_weight = {}
    for wafer in wafers:
        print(wafer)
        nomangle_maps[wafer] = {}
        numanglebolos[wafer] = {}
        varweights[wafer] = {}
        
        for b in goodbolos:
            if pname[b].split('_')[0].upper() == wafer:
                if b not in crapsignalbolos and b not in noiseybolos:
                    band = str(bolo2band[b])
                    noise = bolo_map_variance[b]
                    try:
                        bmt = bolo_maps[b]['T']
                    except Exception:
                        bmt = bolo_maps[b]
                    if ~np.isnan(noise) and noise!=0.:
                        if nominal_angles[b] in nomangle_maps[wafer] and np.shape(bmt) == proper_size:
                            nomangle_maps[wafer][nominal_angles[b]][bolo2band[b]] += np.copy(bmt)*(1.+(-1+1./noise)*IVfactor) 
                            numanglebolos[wafer][nominal_angles[b]][bolo2band[b]] += 1.
                            varweights[wafer][nominal_angles[b]][bolo2band[b]] += 1./noise

                        elif np.shape(bmt) == proper_size:
                            nomangle_maps[wafer][nominal_angles[b]] = {}
                            numanglebolos[wafer][nominal_angles[b]] = {}
                            varweights[wafer][nominal_angles[b]] = {}
                            for band in bands:
                                band = int(band)
                                if band not in nomangle_maps[wafer][nominal_angles[b]]:
                                    nomangle_maps[wafer][nominal_angles[b]][band] = np.zeros(np.shape(bmt))
                                    numanglebolos[wafer][nominal_angles[b]][band] = 0.
                                    varweights[wafer][nominal_angles[b]][band] = 0.
                            nomangle_maps[wafer][nominal_angles[b]][bolo2band[b]] = np.copy(bmt)*(1.+(-1+1./noise)*IVfactor) 
                            numanglebolos[wafer][nominal_angles[b]][bolo2band[b]] = 1.
                            varweights[wafer][nominal_angles[b]][bolo2band[b]] += 1./noise
                        
                else:
                    pass
        if 'all' not in numanglebolos:
            numanglebolos['all'] = {}
        for ang in nomangle_maps[wafer]:
            if ang not in all_wafer_coadds:
                all_wafer_coadds[ang] = {}
                all_wafer_coadd_weight[ang] = {}
                numanglebolos['all'][ang] = {}
            for band in nomangle_maps[wafer][ang]:
                if varweights[wafer][ang][band] >0. and numanglebolos[wafer][ang][band] >0.:
                    nomangle_maps[wafer][ang][band] /=((1.-IVfactor)*numanglebolos[wafer][ang][band] + IVfactor*varweights[wafer][ang][band])
                    if band not in all_wafer_coadds[ang]:
                        all_wafer_coadds[ang][band] = np.copy(nomangle_maps[wafer][ang][band])
                        all_wafer_coadd_weight[ang][band] = 1.

                        numanglebolos['all'][ang][band] = numanglebolos[wafer][ang][band]
                    else:
                        all_wafer_coadds[ang][band] += np.copy(nomangle_maps[wafer][ang][band])
                        all_wafer_coadd_weight[ang][band] += 1.
                        numanglebolos['all'][ang][band] += numanglebolos[wafer][ang][band]
    for ang in all_wafer_coadds:
        for band in all_wafer_coadds[ang]:
            all_wafer_coadds[ang][band] /= all_wafer_coadd_weight[ang][band]
    
                
    #print(numanglebolos)
    return nomangle_maps, numanglebolos, wafer_coadds



def Split_individual_bolos_per_nomang(bolo_maps,wafers,pname,bolo2band,nominal_angles,goodbolos = [], noiseybolos=[],crapsignalbolos=[]):
    ## Returns bolo names in each per-wafer per-band Nomang
    bands = ['90','150','220']
    if len(goodbolos)==0:
        print('No Good Bolos, defaulting to all bolos')
        goodbolos = bolo_maps.keys()
    nomangle_maps = {}
    for wafer in wafers:
        print(wafer)
        nomangle_maps[wafer] = {}
        
        for b in goodbolos:
            if pname[b].split('_')[0].upper() == wafer:
                if b not in crapsignalbolos and b not in noiseybolos:
                    band = str(bolo2band[b])
                    
                    if nominal_angles[b] in nomangle_maps[wafer]:
                        nomangle_maps[wafer][nominal_angles[b]][bolo2band[b]].append(b)
                        
                    else:
                        nomangle_maps[wafer][nominal_angles[b]] = {}
                        
                        for band in bands:
                            band = int(band)
                            nomangle_maps[wafer][nominal_angles[b]][band] = []
                        nomangle_maps[wafer][nominal_angles[b]][bolo2band[b]].append(b)
                        
                else:
                    pass
        
    
                
    return nomangle_maps

                                                    
def Make_PerWafer_PerBand_Coadds_FromIndBoloMaps(bolo_maps,bolo_map_variance,wafers,pname,bolo2band,goodbolos = [], noiseybolos=[],crapsignalbolos=[],IVweighting = True,proper_size=(120,120)):
    ## Returns coadd maps per-wafer per-band
    ## coadds are variance weighted if IVweighting == True
    bands = ['90','150','220']
    if len(goodbolos)==0:
        goodbolos = bolo_maps.keys()
    fitdata = {}
    if IVweighting:
        IVfactor = 1.
    else:
        IVfactor = 0.
    coadd_maps = {}
    numbolos = {}
    varweights = {}
    dx = proper_size[0]//2
    for wafer in wafers:
        coadd_maps[wafer] = {}
        numbolos[wafer] = {}
        varweights[wafer] = {}
        
        for b in goodbolos:
            if pname[b].split('_')[0].upper() == wafer:
                if b not in crapsignalbolos and b not in noiseybolos and b in bolo_map_variance:
                    band = str(bolo2band[b])
                    noise = bolo_map_variance[b]
                    mapx = np.shape(bolo_maps[b]['T'])[0]//2
                    mapy = np.shape(bolo_maps[b]['T'])[1]//2
                    #print np.shape(bolo_maps[b]['T'])
                    #print proper_size
                    if mapx >= dx and mapy>=dx:
                        bmap = bolo_maps[b]['T'][mapx-dx:mapx+dx,mapy-dx:mapy+dx]
                    
                        if band in coadd_maps[wafer]:
                            coadd_maps[wafer][band] += np.copy(bolo_maps[b]['T'])*(1.+(-1+1./noise)*IVfactor) 
                            numbolos[wafer][band] += 1.
                            varweights[wafer][band] += 1./noise
                        
                        else:
                        
                            coadd_maps[wafer][band] = np.copy(bolo_maps[b]['T'])*(1.+(-1+1./noise)*IVfactor) 
                            numbolos[wafer][band]   = 1.
                            varweights[wafer][band] = 1./noise
                else:
                    pass
        
        for band in coadd_maps[wafer]:
            coadd_maps[wafer][band] /=((1.-IVfactor)*numbolos[wafer][band] + IVfactor*varweights[wafer][band])
    #print numanglebolos
    return coadd_maps, numbolos
   




def Get_Bolo_Params_Faster(t,T,Q,U,mask,plotname=''):
    # Bolo Map t given by: 
    # t = 1/2 G * [(2-p)*T + p * cos(2Theta) * Q + p * sin(2Theta) * U]
    # G is gain, p is pol eff, Theta is pol angle of bolometer
    
    ## Generate inverse mask, and mask out nucleus; for generating noise variance.
    inverse_mask = 1.-mask
    dx = np.shape(inverse_mask)[0]//2
    #centermasksize = int(2/mapres) ## in pixels
    #mask[dx-centermasksize:dx+centermasksize,dx-centermasksize:dx+centermasksize] = 0.
    
    T = np.ndarray.flatten(T)
    Q = np.ndarray.flatten(Q)
    U = np.ndarray.flatten(U)
    t = np.ndarray.flatten(t)
    
    mask = np.ndarray.flatten(mask)
    inverse_mask = np.ndarray.flatten(inverse_mask)
    nonzeronoise = t[np.nonzero(inverse_mask)[0]]
    var = np.nanvar(nonzeronoise)
    dof = len(np.nonzero(mask)[0])
    
    ## Fix initial start conditions; set G=2 bc our cal procedure is supposed to work.
    Grange = [2.] 
    prange = np.arange(.5,1.5,.2)
    Thetarange = np.deg2rad(np.arange(0,180,4))
   
    iternum = 0
    while iternum < 3:
        minres = I_Map_Residual(Grange,prange,Thetarange,t,T,Q,U,mask,var,dof)
        minG = minres[0]
        minp = minres[1]
        minTheta = minres[2]
        Grange = Get_Parameter_Range(minG,Grange)
        prange = Get_Parameter_Range(minp,prange)
        Thetarange = Get_Parameter_Range(minTheta,Thetarange)
        #print(Grange,prange,np.rad2deg(Thetarange))
        iternum+=1
    
    return minres[0],minres[1],np.rad2deg(minres[2])    
    


def QU_Decomp(QU,A,B):
    q,u = QU[:]
    return A*q+B*u
    
def Residual(p,x,y):
    res = 1.-np.sqrt(p[1]**2+p[0]**2)
    return y-QU_Decomp(x,p[0],p[1])-res*100.

def Uncert_theta(s1,s2,ds1,ds2):
    dTheta = .5 * (1./(s1**2+s2**2)) *( s2*ds1 - s1*ds2 )
    return dTheta
    

def get_delta_alphas(qmap,umap,dTmap):
    dalpha = .5*np.arctan2(umap-qmap,qmap+umap)
    
    
def center_on_brightest_pix(Map):
    brightestpix = np.where(Map == np.amax(Map))
    centerx = brightestpix[0][0]
    centery = brightestpix[1][0]
    return [centerx,centery]
    
def Correct_AB_Angle(band,nomang,split):
    ## Fit angles from 2018 data
    corr = {90:{'A':3.24,'B':-3.76},150:{'A':.53,'B':-.96},220:{'A':-6.55,'B':6.5}}
    return nomang-np.round(np.deg2rad(corr[band][split]),decimals=3)


<<<<<<< HEAD
=======
def align_maps(map1,map2,shiftrange = 4):
    ## Shift maps around to find shift which
    ## minimizes the residual. map1 is shifted while map2 remains fixed.
    ## map1 and map2 need to have same shape
    xrange = range(-shiftrange,shiftrange+1,1)
    yrange = range(-shiftrange,shiftrange+1,1)
    s = np.shape(map1)[0]//2
    boundary  = s - shiftrange
    res = np.zeros((len(xrange),len(yrange)))
    x = 0
    for i in xrange:
        y = 0
        for j in yrange:
            #print(i,j,x,y)
            #print(s-boundary+i,s+boundary+i,s-boundary+j,s+boundary+j,s-boundary,s+boundary)
            res[x][y] = np.sum(abs(map1[s-boundary+i:s+boundary+i,s-boundary+j:s+boundary+j]-map2[s-boundary:s+boundary,s-boundary:s+boundary]))
            y+=1
        x+=1
    minvals = np.where(res == np.amin(abs(res)))
    if len(minvals[0])==0:
        minvals = np.where(res == -np.amin(abs(res)))
    fp = [xrange[minvals[0][0]],yrange[minvals[1][0]]]
    map1 = np.pad(map1[s-boundary+fp[0]:s+boundary+fp[0],s-boundary+fp[1]:s+boundary+fp[1]],(shiftrange,shiftrange),'constant')
    
    return map1
    






>>>>>>> 884205b41e8b8e6d95b86f88d73a568bba4a10f3


def Fit_Individual_Bolometer_Polarization_PerWafer(bolo_maps,wafer_coadds,mask,bolo2band,nominal_angles,pname,goodbolos = [],mapsize = 30,variable='gain'):
    ## Fits nominal angle maps per wafer, per band to TQU sum maps.  
    bands = ['90','150','220']
    fitdata = {}
    if len(goodbolos) == 0:
        goodbolos = bolo_maps.keys()
    for band in mask:
        ms = np.shape(mask[band])[0]//2
        if ms <= mapsize:
            mapsize=ms
        else:
            mask[band] = mask[band][ms-mapsize:ms+mapsize,ms-mapsize:ms+mapsize]
    
    numbolos = len(goodbolos)
    startofFit = time.time()
    bolonum = 0.
    percentdone = 1.
    plotnum = 0
    for bolo in goodbolos:  
        bolonum+=1.
        if (bolonum/numbolos)*100.>percentdone:
            dt = (time.time()-startofFit)
            print('\n ------------------------------------------------------')
            print('%i Percent Done!'%percentdone)
            print('Elapsed time in calc: %i s.'%dt)
            print('Average time per bolo: %.3f s.'%((dt/float(bolonum))))
            print('Estimated time remaining: %.3f s.'%(dt/float(bolonum)*(numbolos-bolonum)))
            print(' ------------------------------------------------------\n')
            percentdone+=1.
                  
        fitdata[bolo] = {}
        ang = nominal_angles[bolo]
        band = str(bolo2band[bolo])
        coband = band+'GHz'
        wafer = pname[bolo].split('_')[0].upper()
        if wafer in wafer_coadds:
            summapsize = np.shape(wafer_coadds[wafer][coband]['T'])[0]//2
            T = wafer_coadds[wafer][coband]['T'][summapsize-mapsize:summapsize+mapsize,summapsize-mapsize:summapsize+mapsize]
            Q = wafer_coadds[wafer][coband]['Q'][summapsize-mapsize:summapsize+mapsize,summapsize-mapsize:summapsize+mapsize]
            U = wafer_coadds[wafer][coband]['U'][summapsize-mapsize:summapsize+mapsize,summapsize-mapsize:summapsize+mapsize]
            
            if 'T' in bolo_maps[bolo]:
                bmap = bolo_maps[bolo]['T']
            else:
                bmap = bolo_maps[bolo]
            cmapsize = np.shape(bmap)
            if cmapsize[0]//2>= mapsize and cmapsize[1]//2 >= mapsize :
                xsize = cmapsize[0]//2
                ysize = cmapsize[1]//2
                ti = np.nan_to_num(bmap[xsize-mapsize:xsize+mapsize,ysize-mapsize:ysize+mapsize])
                
                p2 = Get_Bolo_Params(band,ang,ti,T,Q,U,mask[str(band)],Variable=variable,plotname='')
                phi = p2[2][0]
            

                if plotnum < 0:
                    plt.figure(figsize=(35,5))
                    plt.subplot(151)
                    plt.imshow(ti*mask[str(band)])
                    plt.title('%s t_i'%bolo)
                    plt.colorbar()
                    plt.subplot(152)
                    plt.imshow(T*mask[str(band)])
                    plt.colorbar()
                    plt.title('T')
                    plt.subplot(153)
                    plt.imshow((ti-T)*mask[str(band)])
                    plt.colorbar()
                    plt.title('t_i-T')
                    plt.subplot(154)
                    plt.imshow((.5* p2[0] * ((2-p2[1])*T + p2[1] * np.cos(2*np.deg2rad(p2[2])) * Q + p2[1] * np.sin(2*np.deg2rad(p2[2])) * U))*mask[str(band)])
                    plt.title("Best Fit t_i")
                    plt.colorbar()
                    plt.subplot(155)
                    plt.imshow((ti-.5* p2[0] * ((2-p2[1])*T + p2[1] * np.cos(2*np.deg2rad(p2[2])) * Q +
p2[1] * np.sin(2*np.deg2rad(p2[2])) * U))*mask[str(band)]/abs(ti),vmin=-.1,vmax=.1)
                    plt.title("Fit Residual")
                    plt.colorbar()
                    plt.savefig('/home/axf295/2019/tmpPlots/%s_FittingPlot.png'%bolo)
                    plt.close('all')
                    plotnum+=1

                if phi > 0:
                    phi%=180.
                elif phi <0.:
                    phi%=-180.
                nomang =  ang #np.rad2deg(ang)
                a = nomang-phi
                if a > 180:
                    a -= 360 
                elif a < -180:   
                    a += 360 
                if (phi > 100. and nomang< 50.) or (nomang>125. and phi < 50.):
                    phi = 180.-phi
                res = nomang-phi
                fitdata[bolo]['rawPhi'] = p2[2]
                fitdata[bolo]['nom'] = nomang
                fitdata[bolo]['fitPhi'] = phi
                fitdata[bolo]['res'] = a
                fitdata[bolo]['G'] = p2[0]
                fitdata[bolo]['rho'] = p2[1]

    return fitdata



        
        
        
        
        
        
