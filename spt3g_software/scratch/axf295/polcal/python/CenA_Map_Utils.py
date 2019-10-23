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
from spt3g import core
import pickle as pickle
#######################################





class CenAMap(object):
    '''
    Class for CenA Map objects; has attributes with resolution, nominal angle, map noise,etc.
    
    '''
    def __init__(self, mapframe, has_pol = False):
        self.map_id = mapframe['Id']
        self.maps = {}
        self.res  = mapframe['T'].res/core.G3Units.arcmin
        self.shape= mapframe['T'].shape
        self.has_pol = has_pol
        if self.has_pol:
            self.maps['T'] = np.asarray(mapframe['T']/mapframe['Wpol'].TT)[::-1]
            self.maps['Q'] = np.asarray(mapframe['Q']/mapframe['Wpol'].QQ)[::-1]
            self.maps['U'] = np.asarray(mapframe['U']/mapframe['Wpol'].UU)[::-1]
        else:
            self.maps['T'] = np.asarray(mapframe['T']/mapframe['Wunpol'].TT)[::-1]
        self.noise  = -1
        self.sn     = -1
        self.nomang = -1
        self.wafer  = -1
        self.band   = -1
        self.numobs = -1
        self.xo     = -1
        self.yo     = -1
    def center_map(self,xo=-1,yo=-1,mapsize=40):
        if xo<0 or yo<0:
            xo = self.xo
            yo = self.yo
        if xo >= mapsize and yo>=mapsize and not np.all(self.maps['T'] == 0.): 
            self.maps['T']=np.asarray(self.maps['T'])[yo-mapsize:yo+mapsize,xo-mapsize:xo+mapsize]
            if self.has_pol:
                self.maps['Q']=np.asarray(self.maps['Q'])[yo-mapsize:yo+mapsize,xo-mapsize:xo+mapsize]
                self.maps['U']=np.asarray(self.maps['U'])[yo-mapsize:yo+mapsize,xo-mapsize:xo+mapsize]
        else:
            self.maps['T']=np.zeros((2*mapsize,2*mapsize))
            if self.has_pol:
                self.maps['Q']=np.zeros((2*mapsize,2*mapsize))
                self.maps['U']=np.zeros((2*mapsize,2*mapsize))
        self.shape = np.shape(self.maps['T'])
        return 
    
    
    def get_brightest_pixel(self):
        if np.all(~np.isfinite(self.maps['T'])):
            return
        smoothT = gf(np.nan_to_num(self.maps['T']), .5/self.res, order=0,  mode='reflect', cval=0.0, truncate=3.0)
        maxpixels = np.where(smoothT== np.amax(np.nan_to_num(smoothT)))
        #if np.all(~np.isnan(smoothT)):
        try:
            self.xo = maxpixels[1][0]
            self.yo = maxpixels[0][0]
        except IndexError:
            pass
        return
    
    
    def interpolate_map_nans(self,interptype='linear'):
        if self.has_pol:
            pollist = ['T','Q','U']
        else:
            pollist = ['T']
        anyinfinite = 0
        for mt in pollist:
            if np.any(~np.isfinite(self.maps[mt])):
                anyinfinite+=1
        if anyinfinite == 0:
            return
                      
        if self.shape[0] > 0 and self.shape[1] > 0:
            x = np.arange(0,self.shape[1])
            y = np.arange(0.,self.shape[0])
            maskedmap = np.ma.masked_invalid(self.maps['T'])
            xx,yy = np.meshgrid(x,y)
            x1 = xx[~maskedmap.mask]
            y1=yy[~maskedmap.mask]

            for mt in pollist:
                newarr = self.maps[mt][~maskedmap.mask]
                try:
                    self.maps[mt] = interpolate.griddata((x1,y1),newarr.ravel(),(xx,yy),method=interptype)
                except Exception:
                    pass
        return 
    
    
    def fft_filter_map(self,maskradius = 6.):
        for mt in self.maps:
            fft = np.fft.fft2(self.maps[mt])
            fftmask = np.zeros(self.shape)
            #radius  = maskradius/self.res
            ## I don't know what 18 corresponds to in frequency space.
            o = self.shape[0]//2
            for i in range(2*o):
                for j in range(2*o):
                    if np.sqrt((o-i)**2+(o-j)**2)<18.0:
                        fftmask[i,j] = 1.

            map_filt = np.fft.ifft2(np.fft.ifftshift(np.fft.fftshift(fft)*fftmask))
            self.maps[mt] = abs(map_filt)
            
        return

def make_nominalangle_coadds(bolo_maps,wafers,goodbolos = []):
    '''
    Coadds individual bolometer maps per-wafer, per-band and per-nominal angle.
    Uses inverse variance weighting in the coadds.
    Also performs a wafer-averaged coadd.
    
    Arguments
    ---------
    bolo_maps : dict
        Dictionary of CenAMap class objects; these contain information about maps
        such as map noise, bolometer properties, and map sn.
        
    wafers : list
        List of all the wafers to use in this analysis.
           
    goodbolos : list
        List of bolometer names that pass cuts.
    
    Returns
    -------
    nomangle_maps : dict
        per-wafer, per-nominal angle, per-band,  coadds as 2D numpy arrays.
        dict struct --> nomangle_maps[wafer][nomang][band]
        key 'all' contains the wafer-averaged coadd
    
    numanglebolos : dict
        per-wafer, per-nominal angle, per-band counts of the number of bolos in each 
        coadd map
    
    all_wafer_coadds : dict
        wafer-averaged coadds per-angle, per-band
        dict struct --> all_wafer_coadds[nomang][band]
    
    '''
    
    BANDS = ['90','150','220']
    if len(goodbolos)==0:
        print('No Good Bolos, returning')
        return {},{},{}
    fitdata       = {}
    nomangle_maps = {}
    numanglebolos = {}
    varweights    = {}
    all_wafer_coadds       = {}
    all_wafer_coadd_weight = {}
    for wafer in wafers:
        nomangle_maps[wafer] = {}
        numanglebolos[wafer] = {}
        varweights[wafer] = {}
        for b in goodbolos:
            if bolo_maps[b].wafer == wafer:
                band   = bolo_maps[b].band
                nomang = bolo_maps[b].nomang

                if ~np.isnan(bolo_maps[b].noise) and bolo_maps[b].noise!=-1 and bolo_maps[b].noise !=0. and band!= -1:
                    if nomang in nomangle_maps[wafer]:
                        if bolo_maps[b].shape == np.shape(nomangle_maps[wafer][nomang][band]):
                            nomangle_maps[wafer][nomang][band] += np.copy(bolo_maps[b].maps['T'])*(1./bolo_maps[b].noise) 
                            numanglebolos[wafer][nomang][band] += 1.
                            varweights[wafer][nomang][band]    += np.copy(1./bolo_maps[b].noise)


                    elif map_is_square(bolo_maps[b].maps['T']):
                        if nomang not in nomangle_maps[wafer]:
                            nomangle_maps[wafer][nomang] = {}
                            numanglebolos[wafer][nomang] = {}
                            varweights[wafer][nomang]    = {}
                        for band in BANDS:
                            if band not in nomangle_maps[wafer][nomang]:
                                nomangle_maps[wafer][nomang][band] = np.zeros(bolo_maps[b].shape)
                                numanglebolos[wafer][nomang][band] = 0.
                                varweights[wafer][nomang][band]    = 0.
                        nomangle_maps[wafer][nomang][bolo_maps[b].band] = np.copy(bolo_maps[b].maps['T'])*(1./bolo_maps[b].noise)
                        numanglebolos[wafer][nomang][bolo_maps[b].band] = 1.
                        varweights[wafer][nomang][bolo_maps[b].band]    += np.copy(1./bolo_maps[b].noise)

               
        if 'all' not in numanglebolos:
            numanglebolos['all'] = {}
        for ang in nomangle_maps[wafer]:
            if ang not in all_wafer_coadds:
                all_wafer_coadds[ang] = {}
                all_wafer_coadd_weight[ang] = {}
                numanglebolos['all'][ang] = {}
            for band in nomangle_maps[wafer][ang]:
                if varweights[wafer][ang][band] > 0. and numanglebolos[wafer][ang][band] > 0.:
                    ## normalize per-wafer,per-ang,per-band maps

                    if band not in all_wafer_coadds[ang]:
                        all_wafer_coadds[ang][band]       = np.copy(nomangle_maps[wafer][ang][band])
                        all_wafer_coadd_weight[ang][band] = np.copy(varweights[wafer][ang][band])
                        numanglebolos['all'][ang][band]   = np.copy(numanglebolos[wafer][ang][band])
                    else:
                        all_wafer_coadds[ang][band]       += np.copy(nomangle_maps[wafer][ang][band])
                        all_wafer_coadd_weight[ang][band] += np.copy(varweights[wafer][ang][band])
                        numanglebolos['all'][ang][band]   += np.copy(numanglebolos[wafer][ang][band])
                    
                    nomangle_maps[wafer][ang][band] /= np.copy(varweights[wafer][ang][band])
                    
    for ang in all_wafer_coadds:
        for band in all_wafer_coadds[ang]:
            all_wafer_coadds[ang][band] /= np.copy(all_wafer_coadd_weight[ang][band])

    
    nomangle_maps['all'] = all_wafer_coadds  
    
    return nomangle_maps, numanglebolos




def make_perband_coadds(bolo_maps,goodbolos = []):
    '''
    Coadds individual bolometer maps per-band.
    Uses inverse variance weighting of the single bolo maps.
    
    
    Arguments
    ---------
    bolo_maps : dict
        Dictionary of CenAMap class objects; these contain information about maps
        such as map noise, bolometer properties, and map sn.
             
    goodbolos : list
        List of bolometer names that pass cuts.
    
    Returns
    -------
    nomangle_maps : dict
        per-wafer, per-nominal angle, per-band,  coadds as 2D numpy arrays.
        dict struct --> nomangle_maps[wafer][nomang][band]
        key 'all' contains the wafer-averaged coadd
    
    numanglebolos : dict
        per-wafer, per-nominal angle, per-band counts of the number of bolos in each 
        coadd map
    
    all_wafer_coadds : dict
        wafer-averaged coadds per-angle, per-band
        dict struct --> all_wafer_coadds[nomang][band]
    
    '''
    
    BANDS = ['90','150','220']
    if len(goodbolos)==0:
        #print('No Good Bolos, returning')
        return {},{},{}
   
    varweights    = {}
    coaddmaps     = {}
    numperband    = {}
    for b in goodbolos:
        band   = bolo_maps[b].band
        if ~np.isnan(bolo_maps[b].noise) and bolo_maps[b].noise!=-1 and bolo_maps[b].noise !=0. and band!= -1:
            if band in coaddmaps:
                if bolo_maps[b].shape == np.shape(coaddmaps[band]):
                    coaddmaps[band]  += np.copy(bolo_maps[b].maps['T'])*(1./bolo_maps[b].noise) 
                    numperband[band] += 1.
                    varweights[band] += np.copy(1./bolo_maps[b].noise)

            else:
                for band in BANDS:
                    if band not in coaddmaps:
                        coaddmaps[band]  = np.zeros(bolo_maps[b].shape)
                        numperband[band] = 0.
                        varweights[band] = 0.
                coaddmaps[bolo_maps[b].band]  = np.copy(bolo_maps[b].maps['T'])*(1./bolo_maps[b].noise)
                numperband[bolo_maps[b].band] = 1.
                varweights[bolo_maps[b].band] += np.copy(1./bolo_maps[b].noise)
    for band in coaddmaps:
        coaddmaps[band]/=varweights[band]
                
    return coaddmaps, numperband

def load_perband_coaddmaps(coadd_file,xo=-1,yo=-1):
    '''
    Loads the per-band coadds saved by autoprocessing.
    Can perform reshaping via resampling if singlebolomapsize is specified.
    
    Arguments
    ---------
    coadd_file : str
        path to the coaddmaps saved by autoprocessing
            
    singlebolomapsize : int
        side length of the single-bolometer maps
    
    Returns
    -------
    summap : dictionary of 2D numpy arrays
        2D numpy arrays stored in a dictionary keyed by band.
        
    '''
    summap = {}
    for fr in core.G3File(coadd_file):
        if 'Wpol' in fr:
            if fr['Id']!='bsmap':
                band = fr['Id'].split('-')[1].split('GHz')[0]
                summap[band] = CenAMap(fr,has_pol=True)
                if xo == -1 and yo == -1:
                    summap[band].get_brightest_pixel()
                    summap[band].center_map()
                else:
                    summap[band].center_map(xo=xo,yo=yo)

    return summap




def calc_source_masked_variance(bolo_map,Noise_Mask=''):

    '''
    Calculate the variance of a masked map
    Checks the shapes of both maps and does the right thing.
    
    Arguments
    ---------
    bolo_map : 2D numpy array
        map to be masked and to of which to find the variance
        
    noise_mask : 2D numpy array
        mask map; multiplies bolo_map to mask source.
    
    Returns
    -------
    variance: float
        variance of the masked map (where bolo_map*noise_mask !=0)
        
    '''

    if Noise_Mask == '':
        Noise_Mask = np.zeros(np.shape(bolo_map))
        Noise_Mask[np.where(bolo_map>1.)]+=1.
        Noise_Mask = 1.- Noise_Mask
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
    try:
        variance = np.nanvar(bolo_map[np.nonzero(bolonoisemask)])
    except IndexError:
        variance = -1
    return variance



def interpolate_map_nans(initmap,interptype='linear'):
    '''
    Interpolates a map over NaNs
    
    Arguments
    ---------
    initmap : 2D numpy array
    
    Returns
    -------
    outmap: 2D-np array or []
        interpolated map if initmap has valid size; otherwise returns empty list.
        
    '''
    if initmap.shape[0] > 0 and initmap.shape[1] > 0:
        x = np.arange(0,initmap.shape[1])
        y = np.arange(0,initmap.shape[0])
        maskedmap = np.ma.masked_invalid(initmap)
        xx,yy = np.meshgrid(x,y)
        x1 = xx[~maskedmap.mask]
        y1=yy[~maskedmap.mask]
        newarr = initmap[~maskedmap.mask]
        GD1 = interpolate.griddata((x1,y1),newarr.ravel(),(xx,yy),method=interptype)
        return GD1
    else:
        return []

def map_is_square(Map):
    '''
    checks whether a map is square or not
    '''
    s = np.shape(Map)
    if s[0] == s[1]:
        return True
    else:
        return False


def downsample_map_nearest_neighboravg(a):
    '''
    Downsamples by finding average of nearest neighbor pixels for every-other 
    pixel of fullsampled map weights are = 1 for central pixel, and .5 for neighbors.
    
    Arguments
    ---------
    a : 2D numpy array
    
    Returns
    -------
    outmap: 2D-np array
        reshaped map 
    
    '''
    ashape = np.shape(a)
    outmap = np.zeros((ashape[0]//2,ashape[1]//2))
    w = np.asarray([[.5,.5,.5],[.5,1.,.5],[.5,.5,.5]])
    for i in range(1,ashape[0]-1,2):
        for j in range(1,ashape[1]-1,2):
            x = (i-1)//2
            y = (j-1)//2
            outmap[x,y] = np.sum((a[i-1:i+2,j-1:j+2]*w))/np.sum(w)
    return outmap


def downsample_map_by_averaging(a,shape):
    '''
    Downsamples a 2D array, a, by averaging.
    --uses nanmean to average
    
    Arguments
    ---------
    a : 2D numpy array
    
    shape: tuple
        shape to resample map into
    
    Returns
    -------
    2D numpy array, reshaped map
    
    '''
    
    if np.shape(a) == shape:
        print('Already the correct shape; returning.')
        return a
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

def check_and_resize_maps(map1,map2='',shape='',dstype='Avg'):
    '''
    Check that two maps are same size, otherwise, upsample/downsample to get
    them equivalent. Assert map is square.
    
    Arguments
    ---------
    map1 : 2D numpy array
        map to down/up-sample
        
    map2 : 2D numpy array or ''
        if 2D-numpy array: map to use for proper shape
        if '' : use shape specified in shape argument
        
    shape: tuple or ''
        if tuple: shape to up/down-sample map into
        if '' : will use map2 shape as shape for map1
       
    dstype: str
        method of downsampling for pixels, (Avg, NN) for average or nearest Neighbors.
        
    Returns
    -------
    map1, map2 
    
    '''
    s1 = np.shape(map1)
    try:
        s2 = np.shape(map2)
        
        if not Map_is_Square(map1) or not Map_is_Square(map2):
            print('One of the maps is not square... aborting')
            return [],[]
    
    except Exception:
        if shape!='':
            s2 = shape
        else:
            print('No shape to compare to! Returning')
            return [],[]
    
    if s1 == s2:
        return map1,map2
    
    if s1[0]> s2[0]:
        if s1[0] != 2*s2[0]:
            map1 = map1.repeat(2,axis=0).repeat(2,axis=1)
        if dstype == 'Avg':
            map1 = Downsample_Map_By_Averaging(map1,s2)
        elif dstype == 'NN':
            map1 = Downsample_Map_Nearest_NeighborAvg(map1)
    else:
        if s2[0] != 2*s1[0]:
            map2 = map2.repeat(2,axis=0).repeat(2,axis=1)
        if dstype == 'Avg':
            map1 = Downsample_Map_By_Averaging(map2,s1)
        elif dstype == 'NN':
            map1 = Downsample_Map_Nearest_NeighborAvg(map2)
    return map1,map2






def calc_cena_sn_using_mask(input_files,obsID,source,log=None,mask_location='./',output_dir='',plot_dir = './tmpplots/',PLOT=False):

    '''
    Calculates the signal to noise of CenA's nucleus, and it's lobes.
    Does nothing fancy but takes mean signal and mean noise in the signal or noise regions.
    
    Can also plot per-band FP maps of SN, as well as per-band histograms.
    
    Arguments
    ---------
    
    input_files : list
        calframe and individual bolomap paths, i.e. [calframe.g3, indbolomaps.g3]
    
    obsID : str
        observation ID, i.e. '64610968'
    
    source : str
        observed source, i.e. 'CenA'
        
    log : file object
        log file to write to, if default None, will open and write to ./log_out.txt
        writes a new log_out_n.txt if log_out_n-1.txt exists.
        
    mask_location : str
        path to mask location
        
    output_dir : str
        path to output directory

    PLOT : bool
        Boolean flag to make per-band S/N FP maps and histograms.
        To plot debugging plots (such as signal and noise regions), manually change params.
    ** TO DO : Add debugging arg for these plots **
    
    5/20/2019 Allen
    
    
    Saves 3 pickle files with bolometer signal and noise properties
    
    '''
    
    if log == None:
        lognum = 0
        if not os.path.exists('./log_out.txt'):
            log = open('log_out.txt','w')
        else:
            while os.path.exists('./log_out_%s.txt'%lognum):
                lognum+=1
            log = open('log_out_%s.txt'%lognum,'w')
                
    log.write('\n###################################################\n ')
    log.write('Welcome to the %s Lobe S/N Calculating Script.\n'%source)
    log.write(' ################################################### \n')   
     
    
    datafile = input_files[1]
    
    if os.path.exists(datafile):
        ind_map_data = core.G3File(datafile)
    else:
        log.write('%s does not exist; returning.\n'%datafile)
        return
    
    log.write('Reading in %s\n'%(datafile))
    
    num_frames = 0
    ## Get total number of bolometers per datafile
    while True:
        try:
            fr = ind_map_data.next()
            num_frames+=1
        except StopIteration:
            break
    
    log.write('Read in %s frames\n'%(num_frames))
    num_bolos = num_frames-3
    
    log.write('Loading Signal and Noise Masks\n')
    
    ## Load Masks
    upperlobe_signal = np.loadtxt(mask_location+'Upper_Lobe_Mask.txt')
    lowerlobe_signal = np.loadtxt(mask_location+'Lower_Lobe_Mask.txt')
    upperlobe_noise = np.loadtxt(mask_location+'Upper_Noise_Mask.txt')
    lowerlobe_noise = np.loadtxt(mask_location+'Lower_Noise_Mask.txt')
    ## Combine lobes and noise regions to reduce number of masks
    signal_mask = upperlobe_signal+lowerlobe_signal
    noise_mask  = lowerlobe_noise+upperlobe_noise
    ## Load the nuclear region mask
    nuclear_mask = np.loadtxt(mask_location+'Nuclear_Mask.txt')
    
    maskshape = np.shape(signal_mask)
    ms = maskshape[0]
    
    log.write('Starting SN Calculations\n')
    
    source_sn = {}
    source_signal = {}

    center_locs = {}
    center_locs['x'] = []
    center_locs['y'] = []

    bolo_noise= {}
    badbolos  = []
    numplots= 0
    bolonum = 0
    analnum = 0
    percent = 0
    ind_map_data = core.G3File(datafile)
    ## Loop over frames to analyze each bolometer.
    while True:
        try:
            fr = ind_map_data.next()
            if 'Wunpol' in fr:
                bolo = fr['Id']
                if bolo == 'map':
                    continue
                mapT = CenAMap(fr)
                mapT.get_brightest_pixel()
                mapT.center_map()
                mapT.interpolate_map_nans()
                 
                res  = mapT.res
                analnum+=1
                ## Print out the percent done (or log in the log file)
                if int(100*(analnum/num_bolos)) > percent:
                    log.write("Analyzing bolo %s out of %s : %s%% \n"%(analnum,num_bolos,percent))
                    log.write(' ------------------------------------------------------\n')
                    percent+=1

                ## Smooth the map with a .5 arcmin sigma ~= 1. arcmin fwhm
                dsx = mapT.shape[0]
                dsy = dsx
                
                ## Make sure masks and data are same shape; or make them so.
                if dsx > ms or dsy > ms:
                    data_im       = np.copy(mapT.maps['T'][dsx//2-ms//2:dsx//2+ms//2, dsy//2-ms//2:dsy//2+ms//2])
                    bolosigmask   = np.copy(signal_mask)
                    bolonoisemask = np.copy(noise_mask)
                    nuclearmask   = np.copy(nuclear_mask)
                elif dsx < ms or dsy < ms:
                    bolosigmask   = np.copy(signal_mask[ms//2-dsx//2:ms//2+dsx//2, ms//2-dsy//2:ms//2+dsy//2])
                    bolonoisemask = np.copy(noise_mask[ms//2-dsx//2:ms//2+dsx//2, ms//2-dsy//2:ms//2+dsy//2])
                    nuclearmask   = np.copy(nuclear_mask[ms//2-dsx//2:ms//2+dsx//2, ms//2-dsy//2:ms//2+dsy//2])
                    data_im       = np.copy(mapT.maps['T'])
                else:
                    data_im       = np.copy(mapT.maps['T'])
                    bolosigmask   = np.copy(signal_mask)
                    bolonoisemask = np.copy(noise_mask)
                    nuclearmask   = np.copy(nuclear_mask)
                    
                xlen,ylen = np.array(np.shape(data_im))
                ## Make sure you're square.
                if xlen == ylen:
                    lobesignal = np.array(data_im[np.nonzero(bolosigmask)],dtype=float)
                    numsignalpixels = len(lobesignal)
                    S   = np.nansum(lobesignal)/numsignalpixels
                    nuclearsignal = np.nansum(np.array(data_im[np.nonzero(nuclearmask)],dtype=float))/len(data_im[np.nonzero(nuclearmask)])
                    S = np.sqrt(S**2+(nuclearsignal)**2)
                    noise  = np.nanvar(data_im[np.nonzero(bolonoisemask)])
                    numnoisepixels = len(data_im[np.nonzero(bolonoisemask)])
                    N = np.sqrt(noise) 
                    
                    if ~np.isnan(S/N) and numplots < 0:

                        #plt.figure('noise hist')
                        #plt.hist(noise,bins=25)
                        #plt.title(bolo)
                        #plt.grid()
                        #plt.savefig('/big_scratch/axf295/tmpPlots/%s_noisehist.png'%bolo)
                        #plt.close()
                        plt.figure()
                        plt.imshow(mapT.maps['T'])
                        plt.colorbar()
                        plt.grid(color='w')
                        plt.savefig(output_dir+'%s_mapT.png'%bolo)

                        
                        plt.figure(bolo)
                        plt.subplot(1,2,1)
                        plt.imshow(data_im*(bolosigmask+nuclearmask))
                        plt.text(20,10, 'Sum: %.1f'%S, fontsize=12,color='w')
                        plt.text(500,20, 'sig-noise: %.1f'%(S-N), fontsize=10,color='w')
                        plt.colorbar()
                        plt.subplot(1,2,2)
                        plt.imshow(data_im*bolonoisemask)
                        plt.text(20,10, 'map noise (std): %.1f'%N, fontsize=12,color='w')
                        plt.xlabel('SNR = %.1f'%(S/N))
                        plt.colorbar()

                        plt.savefig(output_dir+'%s_SNPlot.png'%bolo)
                        plt.close('all')
                        numplots+=1
                    ## Apply some cuts on signal and noise; helps to remove glitches/nans
                    if np.any(data_im[np.nonzero(nuclearmask)] == np.nan) or np.any(noise == np.nan) or np.any(lobesignal==np.nan) or S>1e6 or N > 1e6:
                        source_sn[bolo] = -1.
                        source_signal[bolo] = [-1.,-1]
                        badbolos.append(bolo)
                    else:
                        if N!=0.0 and np.isfinite(S/N) and S/N < 1e6 and mapT.xo>0 and mapT.yo>0:
                            source_sn[bolo] = S/N
                            source_signal[bolo] = [nuclearsignal, N]
                            bolo_noise[bolo] = N
                            center_locs['x'].append(mapT.xo)
                            center_locs['y'].append(mapT.yo)
                        else:
                            source_sn[bolo] = -1.
                            source_signal[bolo] = [-1,-1.]
                            badbolos.append(bolo)
                            bolo_noise[bolo] = -1
        
        except StopIteration:
            break
      
    plt.figure()
    plt.hist(center_locs['x'],bins=201,histtype='step',label='x')
    plt.hist(center_locs['y'],bins=201,histtype='step',label='y')
    plt.grid()
    plt.legend()
    plt.savefig(output_dir+'x_and_y_center_hists.png')
    
    ## Write center locations to file
    with open(output_dir+'%s_center_location_Dictionary.pkl'%source, 'wb') as handle:
        pickle.dump(center_locs, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    log.write(str(len(badbolos))+ ' bolos dropped due to zeros or nans in map or SN\n')
    ## Write the Signal, noise, and SN to dictionaries
    with open(output_dir+'%s_SN_Dictionary.pkl'%source, 'wb') as handle:
        pickle.dump(source_sn, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    with open(output_dir+'Bolo_Noise_Dictionary.pkl', 'wb') as handle:
        pickle.dump(bolo_noise, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(output_dir+'Source_Signal_Dictionary.pkl', 'wb') as handle:
        pickle.dump(source_signal, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    if PLOT:
        log.write('Plotting the Lobe SN histograms\n')    

        offline_cal = input_files[0]

        if not os.path.exists(offline_cal):
            log.write('Defaulting to latest boloprops file\n')
            offline_cal = '/spt/user/production/calibration/boloproperties/62679603.g3'
        bolo_band = {}
        for frame in core.G3File(offline_cal):
            for bolo in frame['BolometerProperties'].keys():
                if not np.isnan(frame['BolometerProperties'][bolo].band):
                    bolo_band[bolo] = str(int(frame['BolometerProperties'][bolo].band/core.G3Units.GHz))
                else:
                    bolo_band[bolo] = -1


        bandcolors = {}
        bandcolors['90']  = 'b'
        bandcolors['150'] = 'g'
        bandcolors['220'] = 'r'


        SourceSN = {}
        SourceSN['90']  = []
        SourceSN['150'] = []
        SourceSN['220'] = []

        SourceS = {}
        SourceS['90']  = []
        SourceS['150'] = []
        SourceS['220'] = []
        for bolo in source_sn:
            if source_sn[bolo]>0. and bolo in bolo_band:
                if bolo_band[bolo]!=-1:
                    SourceSN[bolo_band[bolo]].append(source_sn[bolo])
                    if bolo in source_signal:
                        if source_signal[bolo][0]>0.:
                            SourceS[bolo_band[bolo]].append(source_signal[bolo][0])
        

        log.write('bolos SN added to dict\n')


        plt.figure('%s Source SN'%obsID)
        plt.title('%s Source SN Histogram'%obsID)
        SN_vmax = {}
        for band in SourceSN:
            medianSN = np.median(SourceSN[band])
            stdSN = np.std(SourceSN[band])
            SN_vmax[band] = medianSN+2*stdSN
            BINS = np.linspace(medianSN-3*stdSN,medianSN+3*stdSN,101)
            plt.hist(SourceSN[band],bins=BINS,color=bandcolors[band],alpha=.5,label=band+' GHz')

        plt.xlabel('%s SN'%source,fontsize=16)
        plt.legend()
        plt.savefig(output_dir+'%s_%s_SN_Hist_byband.png'%(obsID,source))

        plt.close()
        
        off = core.G3File(offline_cal).next()
        good_pixels= {}
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

                        if bolo_band[bolo] == band and bolo in source_sn:
                            bolomark = (2,0,NA)
                            plt.scatter(x,y,marker=bolomark,
                                        c=source_sn[bolo],vmin=0,
                                        vmax=VMAX,cmap='plasma',edgecolor='face')

                            numbolos+=1
            log.write('%s, %s GHz bolos plotted\n'%(numbolos,band))

            if numbolos != 0:
                plt.colorbar(label='Source SN')
            plt.title('Obs %s, %s GHz'%(obsID,band))
            plt.savefig(output_dir+'%s_%sLobe_SN_FP_map_%sGHz.png'%(obsID,source,band))
            plt.close()

        S_vmax = {}
        plt.figure('%s Source Signal'%obsID)
        plt.title('%s Source Signal Histogram'%obsID)
        for band in SourceS:
            medianS = np.median(SourceS[band])
            stdS = np.std(SourceS[band])
            S_vmax[band] = medianS+2*stdS
            BINS = np.linspace(medianS-3*stdS,medianS+3*stdS,101)
            plt.hist(SourceS[band],bins=BINS,color=bandcolors[band],alpha=.5,label=band+' GHz')
        plt.xlabel('%s Signal'%source,fontsize=16)
        plt.legend()
        plt.savefig(output_dir+'%s_%s_Signal_Hist_byband.png'%(obsID,source))
        plt.close()    
        for band in SourceS:
            plt.figure('Source Signal, %s GHz'%band,figsize=(10,10))
            numbolos = 0
            VMAX = S_vmax[band]

            for bolo in off['BolometerProperties'].keys():
                if np.isfinite(off['BolometerProperties'][bolo].pol_angle):
                    x = off['BolometerProperties'][bolo].x_offset
                    y = off['BolometerProperties'][bolo].y_offset
                    NA =int(np.rad2deg(off['BolometerProperties'][bolo].pol_angle%np.pi))

                    if bolo in bolo_band:
                        if bolo_band[bolo] == band and bolo in source_signal:
                            bolomark = (2,0,NA)
                            plt.scatter(x,y,marker=bolomark,c=source_signal[bolo][0],vmin=0,vmax=VMAX,cmap='plasma',edgecolor='face')

                            numbolos+=1
            log.write('%s, %s GHz bolos plotted\n'%(numbolos,band))
            if numbolos != 0:
                plt.colorbar(label='Source Signal [mK/pix]')
            plt.title('Obs %s, %s GHz'%(obsID,band))
            plt.savefig(output_dir+'%s_%s_Signal_FP_map_%sGHz.png'%(obsID,source,band))
            plt.close()

    log.write('\n ################################################### \n')
    log.write('Congrats! You are finished calculating %s masked SN.\n'%(source))
    log.write(' ################################################### \n')
    return log