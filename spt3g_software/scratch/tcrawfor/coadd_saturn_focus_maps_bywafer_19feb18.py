#from spt3g import core, dfmux, std_processing
##from spt3g.scratch.tcrawfor import tctools
#from spt3g.util import tctools
#import numpy as np
#import pickle, glob
#from scipy import ndimage
#
#obsid_min = 35296574
#obsid_max = 35315694
#obsid_fp = 35296574
#
#files = glob.glob('/spt/user/production/calibration/saturn-pixelraster/maps/*.g3')
#obsids = []
#for file1 in files:
#    obsids.append(np.int(file1.split('/')[-1].split('.')[0]))
#obsids = np.asarray(obsids)
#obsids.sort()
#obsids = obsids[np.where(np.logical_and(obsids >= obsid_min,obsids <= obsid_max))]
#
## get offsets first. would rather use independent fast_point, but none seem to have been processed
#xoffdict = {}
#yoffdict = {}
#mapdict = {}
#file1 = '/spt/user/production/calibration/saturn-pixelraster/maps/'+str(obsid_fp)+'.g3'
#f1 = core.G3File(file1)
#for frame in f1:
##    if frame.type is core.G3FrameType.Calibration:
##        bp = frame['BolometerProperties']
#    if frame.type is core.G3FrameType.Map:
#        name = frame['Id']
#        if name == 'bsmap':
#            continue
#        if 'Wunpol' in frame:
#            twt = np.asarray(frame['Wunpol'].TT)
#            twt_orig = twt.copy()
#            twt[np.where(twt == 0.)] = 1.
#        map_weighted = np.asarray(frame['T'])
#        reso_arcmin = (frame['T']).res/core.G3Units.arcmin
#        map_unw = map_weighted.copy()
#        map_unw /= twt
#        ny,nx = np.shape(map_unw)
##        amap_sm = ndimage.gaussian_filter(-map_unw,2./reso_arcmin)
#        amap_sm = -1.*map_unw.copy()
##        ycenter, xcenter = np.unravel_index(np.argmax(amap_sm),[ny,nx])
#        ycenter, xcenter = np.unravel_index(np.nanargmax(amap_sm),[ny,nx])
#        xoffdict[name] = xcenter - nx/2.
#        yoffdict[name] = ycenter - ny/2.
##        mapdict[name] = map_unw

# now make coadds
wafmaps = {}
uwn = np.unique(np.asarray([bp[name].wafer_id for name in xoffdict.keys()]))
#for obsid1 in obsids:
for obsid1 in obsids[3:]:
    file1 = '/spt/user/production/calibration/saturn-pixelraster/maps/'+str(obsid1)+'.g3'
    print(file1)
    dtemp = {}
    wtemp = {}
    for uwn1 in uwn:
        dtemp[uwn1] = np.zeros([ny,nx])
        wtemp[uwn1] = np.zeros([ny,nx])
    dtemp['all'] = np.zeros([ny,nx])
    wtemp['all'] = np.zeros([ny,nx])
    f1 = core.G3File(file1)
    for frame in f1:
#        if frame.type is core.G3FrameType.Calibration:
#            bp = frame['BolometerProperties']
        if frame.type is core.G3FrameType.Map:
            name = frame['Id']
            if name == 'bsmap':
                continue
            if bp[name].band != 2200.:
#            if bp[name].band != 1500.:
#            if bp[name].band != 900.:
                continue
            if 'Wunpol' in frame:
                twt = np.asarray(frame['Wunpol'].TT)
                twt_orig = twt.copy()
                twt[np.where(twt == 0.)] = 1.
            map_weighted = np.asarray(frame['T'])
            mwtemp = map_weighted.copy()
            twtemp = twt.copy()
            whbad = np.where(np.isfinite(map_weighted) == False)
            mwtemp[whbad] = 0.
            twtemp[whbad] = 0.
            whbad2 = np.where(np.abs(mwtemp) > 1e-10)
            mwtemp[whbad2] = 0.
            twtemp[whbad2] = 0.
            reso_arcmin = (frame['T']).res/core.G3Units.arcmin
            try:
                yoff = np.int(np.round(yoffdict[name]))
                xoff = np.int(np.round(xoffdict[name]))
                wafid = bp[name].wafer_id
                wtemp[wafid] += np.roll(np.roll(twtemp,-yoff,0),-xoff,1)
                dtemp[wafid] += np.roll(np.roll(mwtemp,-yoff,0),-xoff,1)
                wtemp['all'] += np.roll(np.roll(twtemp,-yoff,0),-xoff,1)
                dtemp['all'] += np.roll(np.roll(mwtemp,-yoff,0),-xoff,1)
            except:
                pass

    for wafid in dtemp.keys():
        wt2 = wtemp[wafid]
        wt2[np.where(wt2 == 0.)] = 1.
        dtemp[wafid] /= wt2
        dtemp[wafid] *= -1.

    wafmaps[obsid1] = dtemp

    print(notavariable)
