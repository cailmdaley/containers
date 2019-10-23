import argparse
import sky_local as sky
import numpy as np
import matplotlib.pyplot as plt
import glob

map0 = np.load('/scratch/sri/3g_maps/KFerguson_201811xx/dec'+str(args.source)+'_150GHz_Tmaps_0.npz')['MAP']

ptsrc_list = np.loadtxt('/home/tcrawfor/code/spt3g_software/scratch/tcrawfor/at20g_list_kferguson.txt', usecols = [1,2,3], unpack = 1)

#create array of pixel locations of point sources
flipped_pixels = sky.ang2Pix(ptsrc_list,[0,dec0],args.resol_arcmin,map0.shape,proj=args.proj)

pixels = [map0.shape[0] - flipped_pixels[0],map0.shape[1] - flipped_pixels[1]]

obsid = j[38:46]
map = np.load(j)['MAP']
#	bright_pixel = np.asarray(np.unravel_index(np.argmax(cutout),cutout.shape))
ra_diffs = []
ra_orig = []
ra_unc = []
dec_diffs = []
dec_orig = []
dec_unc = []

for k in bright:
        y_cen,x_cen,lum = pixels_in_field[k]
        y_min = int(y_cen - (args.width_in_pixels/2))
        y_max = int(y_cen + (args.width_in_pixels/2))
        x_min = int(x_cen - (args.width_in_pixels/2))
        x_max = int(x_cen + (args.width_in_pixels/2))
        cutout = map[y_min:y_max, x_min:x_max]
        if all([cutout.size != 0, np.isnan(map[y_cen][x_cen]) == False]):
                bright_pixel_c = np.asarray(np.unravel_index(np.argmax(cutout),cutout.shape))
                bright_pixel = [bright_pixel_c[0] + y_min, bright_pixel_c[1] + x_min]
                reflipped_cen = [map0.shape[0] - y_cen, map0.shape[1] - x_cen]
                flipped_bright = [map0.shape[0] - bright_pixel[0], map0.shape[1] - bright_pixel[1]]

                radec_cen = sky.pix2Ang(reflipped_cen,np.asarray([0,dec0]),args.resol_arcmin,map0.shape,proj=args.proj)
                radec_bright = sky.pix2Ang(flipped_bright,np.asarray([0,dec0]),args.resol_arcmin,map0.shape,proj=args.proj)
                diff_ra = radec_bright[0] - radec_cen[0]
                diff_dec = radec_bright[1] - radec_cen[1]

                ra_diffs.append(diff_ra)
                ra_orig.append(radec_cen[0])
                dec_diffs.append(diff_dec)
                dec_orig.append(radec_cen[1])
                
                var_ra_diffs = []
                var_dec_diffs = []
                N = 10
                for this in range (0,N):
                        temp = np.copy(cutout)[:]
                        #				print(cutout[bright_pixel_c[0]-1:bright_pixel_c[0]+2,bright_pixel_c[1]-1:bright_pixel_c[1]+2])
                        temp[bright_pixel_c[0],bright_pixel_c[1]] = np.nan
                        rms = np.sqrt(np.nanmean(np.square(temp)))
                        cutout_variance = np.asarray(np.random.normal(0,rms,temp.shape))
                        var_cutout = np.asarray(cutout) + cutout_variance
                        var_bright_pixel_c = np.asarray(np.unravel_index(np.argmax(var_cutout),var_cutout.shape))
                        var_bright_pixel = [var_bright_pixel_c[0] + y_min, var_bright_pixel_c[1] + x_min]
                        var_flipped_bright = [map0.shape[0] - var_bright_pixel[0], map0.shape[1] - var_bright_pixel[1]]
                        var_radec_bright = sky.pix2Ang(var_flipped_bright,np.asarray([0,dec0]),args.resol_arcmin,map0.shape,proj=args.proj)
                        var_diff_ra = var_radec_bright[0] - radec_cen[0]
                        var_diff_dec = var_radec_bright[1] - radec_cen[1]
                        var_ra_diffs.append(var_diff_ra)
                        var_dec_diffs.append(var_diff_dec)
                        #			print(np.asarray(var_ra_diffs) - diff_ra); print(np.asarray(var_dec_diffs) - diff_dec);
                        #			print(np.std(var_ra_diffs)); print(np.std(var_dec_diffs))
                        ra_unc.append(np.std(var_ra_diffs)/np.sqrt(N))
                        dec_unc.append(np.std(var_dec_diffs)/np.sqrt(N))

                        #		print(ra_unc); print(dec_unc)
                        plt.subplot(221); plt.errorbar(x=ra_orig,y=ra_diffs,yerr=ra_unc,fmt='bo'); plt.title('delta-RA vs. RA (deg)')#; plt.ylim(y_lim_RA)
                        plt.subplot(222); plt.errorbar(x=ra_orig,y=dec_diffs,yerr=dec_unc,fmt='bo'); plt.title('delta-dec vs. RA (deg)')#; plt.ylim(y_lim_dec)
                        plt.subplot(223); plt.errorbar(x=dec_orig,y=ra_diffs,yerr=ra_unc,fmt='bo'); plt.title('delta-RA vs. dec (deg)')#; plt.ylim(y_lim_RA)
                        plt.subplot(224); plt.errorbar(x=dec_orig,y=dec_diffs,yerr=dec_unc,fmt='bo'); plt.title('delta-dec vs. dec (deg)')#; plt.ylim(y_lim_dec)
                        mng = plt.get_current_fig_manager()
                        mng.resize(*mng.window.maxsize()) 
                        plt.savefig('/home/kferguson/pointing_offsets/dec'+str(args.source)+'_obs_'+str(obsid)+'.png')
                        plt.close()
