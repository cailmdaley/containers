#from spt3g import core, std_processing, dfmux, calibration
#from spt3g.mapmaker import mapmakerutils as mm
#from spt3g.mapmaker import summingmaps as sm
#from scipy import ndimage
#
#f1=core.G3File('/spt/user/ddutcher/coadds/weight1-3Hz_noW201wafCMpoly19mhpf300_51.g3')
#
#is_init = False
#for frame in f1:
#        if '150' in frame['Id']:
#                if not(is_init):
#                        tmap0 = frame['T']
#                        tmap150 = np.asarray(frame['T'])
#                        twts150 = np.asarray(frame['Wpol'].TT)
#                        is_init = True
#                else:
#                        tmap150 += np.asarray(frame['T'])
#                        twts150 += np.asarray(frame['Wpol'].TT)
#
#twts150 = twts150 + 0.0001
#tmap150 /= twts150
mapshape = np.shape(tmap150)
map = tmap150
nside = 20
thisstd = np.std(map[mapshape[0]/2-100:mapshape[0]/2+100,mapshape[1]/2-100:mapshape[1]/2+100])

ptsrc_list = np.loadtxt('/home/tcrawfor/code/spt3g_software/scratch/tcrawfor/at20g_list_kferguson.txt', usecols = [1,2,3], unpack = 1)
ras = np.asarray(ptsrc_list[0,:])
decs = np.asarray(ptsrc_list[1,:])

#create array of pixel locations of point sources
pixels = np.asarray(tmap0.angles_to_pixels(ras*core.G3Units.deg,decs*core.G3Units.deg))
xpixels = np.mod(pixels,mapshape[1])
ypixels = pixels/mapshape[1]
#whg = np.where(np.logical_and(np.abs(ypixels-mapshape[0]/2) < mapshape[0]/2,np.abs(xpixels-mapshape[1]/2) < mapshape[1]/2))
#ras_neg = ras.copy()
#ras_neg[np.where(ras > 180.)] -= 360.
whg = np.where(np.logical_and(np.abs(decs+56.) < 12.,np.abs(ras_neg) < 40.))
#whg = np.where(np.logical_and(np.abs(ypixels-890.) < 50.,np.abs(xpixels-1230.) < 50.))

##	bright_pixel = np.asarray(np.unravel_index(np.argmax(cutout),cutout.shape))
#ra_diffs = []
#ra_orig = []
#ra_unc = []
#dec_diffs = []
#dec_orig = []
#dec_unc = []
x_diffs = []
y_diffs = []
ras4diffs = []
decs4diffs = []

for j in np.arange(len(whg[0])):
        k = whg[0][j]
        y_cen = ypixels[k]
        x_cen = xpixels[k]
        y_min = int(y_cen - (nside/2))
        y_max = int(y_cen + (nside/2))
        x_min = int(x_cen - (nside/2))
        x_max = int(x_cen + (nside/2))
        #        cutout = ndimage.gaussian_filter(map[y_min:y_max, x_min:x_max],fwhm)
        cutout = map[y_min:y_max, x_min:x_max]
        #        if all([cutout.size != 0, np.isnan(map[y_cen][x_cen]) == False]):
        if np.max(cutout)/thisstd > 10.:
                bright_pixel_c = np.asarray(np.unravel_index(np.argmax(cutout),cutout.shape))
                bright_pixel = [bright_pixel_c[0] + y_min, bright_pixel_c[1] + x_min]
                x_diffs.append(bright_pixel[1] - x_cen)
                y_diffs.append(bright_pixel[0] - y_cen)
                ras4diffs.append(ras[k])
                decs4diffs.append(decs[k])
                
