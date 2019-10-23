#from spt3g import core, std_processing, dfmux, calibration
#from spt3g.mapmaker import mapmakerutils as mm
#from spt3g.mapmaker import summingmaps as sm
#from scipy import ndimage
#
#f1 = core.G3File('/spt/user/ddutcher/coadds/proj0_tonly_150GHz_0.25res_47.g3')
#frame = f1.next()
#mm.RemoveWeightModule(frame)
#map = np.asarray(frame['T'])
#mapshape = np.shape(map)
#nside = 60
#reso_arcmin = frame['T'].res/core.G3Units.arcmin
#fwhm = np.int(1.5/reso_arcmin)
#thisstd = np.std(map[mapshape[0]/2-100:mapshape[0]/2+100,mapshape[1]/2-100:mapshape[1]/2+100])
thisstd = np.std(ndimage.gaussian_filter(map[mapshape[0]/2-100:mapshape[0]/2+100,mapshape[1]/2-100:mapshape[1]/2+100],fwhm))
#
#ptsrc_list = np.loadtxt('/home/tcrawfor/code/spt3g_software/scratch/tcrawfor/at20g_list_kferguson.txt', usecols = [1,2,3], unpack = 1)
#ras = np.asarray(ptsrc_list[0,:])
#decs = np.asarray(ptsrc_list[1,:])

#create array of pixel locations of point sources
pixels = np.asarray(frame['T'].angles_to_pixels(ras*core.G3Units.deg,decs*core.G3Units.deg))
xpixels = np.mod(pixels,mapshape[1])
ypixels = pixels/mapshape[1]
#whg = np.where(np.logical_and(np.abs(ypixels-mapshape[0]/2) < mapshape[0]/2,np.abs(xpixels-mapshape[1]/2) < mapshape[1]/2))
ras_neg = ras.copy()
ras_neg[np.where(ras > 180.)] -= 360.
#whg = np.where(np.logical_and(np.abs(decs+56.) < 12.,np.abs(ras_neg) < 40.))
whg = np.where(np.logical_and(np.abs(decs+56.) < 14.,np.abs(ras_neg) < 50.))
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
ra_diffs = []
dec_diffs = []
ras4diffs = []
decs4diffs = []

for j in np.arange(len(whg[0])):
#for j in np.arange(10):
        k = whg[0][j]
        y_cen = ypixels[k]
        x_cen = xpixels[k]
        y_min = int(y_cen - (nside/2))
        y_max = int(y_cen + (nside/2))
        x_min = int(x_cen - (nside/2))
        x_max = int(x_cen + (nside/2))
        cutout = ndimage.gaussian_filter(map[y_min:y_max, x_min:x_max],fwhm)
#        cutout = map[y_min:y_max, x_min:x_max]
        #        if all([cutout.size != 0, np.isnan(map[y_cen][x_cen]) == False]):
#        if np.max(cutout)/thisstd > 10.:
        if np.max(cutout)/thisstd > 3.:
                bright_pixel_c = np.asarray(np.unravel_index(np.argmax(cutout),cutout.shape))
                bright_pixel = [bright_pixel_c[0] + y_min, bright_pixel_c[1] + x_min]
                x_diffs.append(x_cen - bright_pixel[1])
                y_diffs.append(y_cen - bright_pixel[0])
                radec_temp = frame['T'].pixel_to_angle(bright_pixel[1],bright_pixel[0])/core.G3Units.deg
                ra_diffs.append(radec_temp[0] - ras_neg[k])
                dec_diffs.append(radec_temp[1] - decs[k])
                ras4diffs.append(ras_neg[k])
                decs4diffs.append(decs[k])
#                if np.abs(ras[k]-325.4) < 0.2 and np.abs(decs[k]+64.187) < 0.2:
#                        print(notavariable)
                
x_diffs = np.asarray(x_diffs)
y_diffs = np.asarray(y_diffs)
ra_diffs = np.asarray(ra_diffs)
dec_diffs = np.asarray(dec_diffs)
ras4diffs = np.asarray(ras4diffs)
decs4diffs = np.asarray(decs4diffs)
figure()
plot(decs4diffs,ra_diffs*3600.,'o')
ylim(0,800)
xlabel('dec [deg]')
ylabel('measured RA minus AT20G RA [arcsec]')
savefig('dra_vs_dec_with_3g_software.png')
clf()
plot(decs4diffs,dec_diffs*3600.,'o')
ylim(0,120)
xlabel('dec [deg]')
ylabel('measured dec minus AT20G dec [arcsec]')
savefig('ddec_vs_dec_with_3g_software.png')
