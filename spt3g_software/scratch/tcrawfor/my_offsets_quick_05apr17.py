from scipy import ndimage

xoffs_tc = {}
yoffs_tc = {}

xvec = np.arange(360)

for name in mapdict4.keys():
    amap = -mapdict4[name]
    amap[:,0:80] = 0.
    amap_sm = ndimage.gaussian_filter(amap,4)
    mapshape = np.shape(amap)
    ycenter, xcenter = np.unravel_index(np.argmax(np.abs(amap_sm)),[mapshape[0],mapshape[1]])
    yoffs_tc[name]= ycenter
    xoffs_tc[name]= xcenter
    if xcenter > 10 and xcenter < 350 and ycenter > 10 and ycenter < 350:
        cutout = amap[ycenter-10:ycenter+10,xcenter-10:xcenter+10]
        ycut = xvec[ycenter-10:ycenter+10]
        xcut = xvec[xcenter-10:xcenter+10]
        ycm = np.sum(np.sum(cutout,0)*ycut)/np.sum(cutout)
        xcm = np.sum(np.sum(cutout,0)*xcut)/np.sum(cutout)
        yoffs_tc[name]= ycm
        xoffs_tc[name]= xcm
    
