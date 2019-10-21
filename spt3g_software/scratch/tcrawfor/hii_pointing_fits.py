import numpy as np
import cPickle as pkl
import glob
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math as m
from astropy.modeling import models, fitting
from spt3g import core, coordinateutils
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Ellipse

# stolen from Matt Young

######################
# Settings...
######################
verbose = True
output_dir = '/spt/user/tcrawfor/public/'
pkl_savefile = output_dir+'hii_pointing_fits.pkl'
pdf_savefile = 'hii_pointing_fits.pdf'
obs_dir = '/spt/data/bolodata/downsampled/' # Where the observation data lives
map_dir = '/spt/user/production/calibration/' # Where the autoprocessed maps live
ob_start = 50000000#36115066
ob_stop  = 60000000#36464817
fit_area = 40 # arcmins squared
cmap = cm.bone_r
cmap_res = cm.coolwarm
vmin, vmax = -0.1, 1
vmin_res, vmax_res = -0.1, 0.1 

######################
# Setting up Variables...
######################
offset_dict = {}
benchpos_dict = {}
sources = ['RCW38','MAT5A']
freqs = ['90GHz','150GHz','220GHz']
fig = plt.figure(123,figsize=(15,10))

#######################
# Finding autoprocessed maps and fitting offsets...
#######################
for _source in sources:
    offset_dict[_source] = {}
    map_source_dir = map_dir+_source+'/maps/' 
    map_files = np.sort([_file for _file in glob.glob(map_source_dir+'*.g3') if ob_start <= int(_file.split('/')[-1].split('.')[0]) <= ob_stop])
    pdf = PdfPages(output_dir+_source+'_'+pdf_savefile)

    for _file in map_files:
        plt.clf()
        maps = list(core.G3File(_file))
        obsID = _file.split('/')[-1].split('.')[0]
        if obsID == '34563908': continue # Glitched map
        offset_dict[_source][obsID] = {}
        
        # Load and save the bench information (useful for early observations with changing focus)
        ob_source_dir = obs_dir+_source+'/'+obsID+'/'
        ob_file = core.G3File(ob_source_dir+'0000.g3').next()
        benchpos_dict[obsID] = ob_file['BenchCommandedPosition']
        ob_time = ob_file['ObservationStart']

        # Skip the pdf page if no valid fits
        skip_page = 0

        for _i, _freq in enumerate(freqs):
            # Load the map and properties
            current_map = maps[_i]['T']/maps[_i]['Wunpol'].TT
            resolution = current_map.res*(180./np.pi)*60. # in arcminutes
            cutout_mapsize = fit_area / resolution
            y_size, x_size = np.shape(current_map)
            x_map_centre = (x_size-1)/2.
            y_map_centre = (y_size-1)/2.
            xmin = int(m.floor(x_map_centre - cutout_mapsize/2.))
            xmax = int(m.ceil(x_map_centre + cutout_mapsize/2.))
            ymin = int(m.floor(y_map_centre - cutout_mapsize/2.))
            ymax = int(m.ceil(y_map_centre + cutout_mapsize/2.))

            # Set up the model for fitting
            y, x = np.mgrid[0:y_size,0:x_size]
            map_cutout = np.asarray(current_map)[ymin:ymax,xmin:xmax]
            if np.isnan(map_cutout).any():
                skip_page += 1
                continue # If the map cutout contains nans, lets skip it for now
            y_guess, x_guess = np.unravel_index(map_cutout.argmax(), map_cutout.shape)
            g_init = models.Gaussian2D(amplitude=0.1, x_mean=x_guess+xmin, y_mean=y_guess+ymin, x_stddev=2, y_stddev=2)
            fit_g = fitting.LevMarLSQFitter()
            
            # Fit and calculate source displacement
            g = fit_g(g_init, x[ymin:ymax,xmin:xmax], y[ymin:ymax,xmin:xmax], map_cutout)
            fit_amplitude = g.amplitude.value
            x_centroid = g.x_mean.value
            y_centroid = g.y_mean.value
            x_stddev =g.x_stddev.value
            y_stddev = g.y_stddev.value
            delta_RA = resolution*(x_centroid - (x_size/2. - 0.5)) # displacement in arcmins. 0.5 pixel edge offset
            delta_dec = resolution*(y_centroid - (y_size/2. - 0.5)) # displacement in arcmins. 0.5 pixel edge offset
            pointing_error = np.sqrt(delta_RA**2+delta_dec**2)
            
            if verbose:
                print '\nFitting %s obsID %s %s' %(_source,obsID,_freq)
                print 'X centroid:', x_centroid, '\tdelta RA: %.2f arcmins' %delta_RA
                print 'Y centroid:', y_centroid, '\tdelta dec: %.2f arcmins' %delta_dec
                print '- Total pointing error: %.2f arcmins' %pointing_error
            
            # Need a way of identifying  bad fits and potentially correcting them
            #if fit_amplitude == 0.1 and x_centroid == 180 and y_centroid == 180:
            #    offset_dict[_source].pop(obsID, None)
            #    continue
            
            # Store results
            offset_dict[_source][obsID][_freq] = {'amplitude'  : fit_amplitude,
                                                  'x_centroid' : x_centroid, #pixel
                                                  'y_centroid' : y_centroid, #pixel
                                                  'delta_RA'   : delta_RA,   #arcmins
                                                  'delta_dec'  : delta_dec,  #arcmins
                                                  'total_error': pointing_error} #arcmins

            # Plot the result
            ax = plt.subplot(2,3,_i+1)
            im = plt.imshow(np.asarray(current_map),interpolation='none',cmap=cmap,vmin=vmin,vmax=vmax)
            plt.xlim([xmin,xmax])
            plt.ylim([ymin,ymax])
            plt.title(_freq)
            # Annotate the plot with fit parameters and overplot the fit
            text_ref = [xmin+fit_area*0.1, ymax-fit_area*0.65]
            ax.text(text_ref[0],text_ref[1],r"""
$A: %.2f$
$\sigma_{RA}: %.2f '$
$\sigma_{dec}: %.2f '$
$\Delta{RA}: %.2f '$
$\Delta{dec}: %.2f '$""" %(fit_amplitude,x_stddev*resolution,y_stddev*resolution,delta_RA,delta_dec))
            e = Ellipse(xy=(x_centroid,y_centroid), width=2*x_stddev, height=2*y_stddev, edgecolor='r', fc='None', lw=1)
            ax.add_artist(e)
            # Plot centroid and crosshairs
            plt.plot(x_centroid,y_centroid,'r.',markersize=2)
            plt.plot([xmin,xmax],[y_map_centre,y_map_centre],'k--')
            plt.plot([x_map_centre,x_map_centre],[ymin,ymax],'k--')
            # Manually set the ticks
            x_ticks = [x_map_centre-(fit_area/2.)/resolution,x_map_centre-(fit_area/4.)/resolution,x_map_centre,x_map_centre+(fit_area/4.)/resolution,x_map_centre+(fit_area/2.)/resolution]
            y_ticks = [y_map_centre-(fit_area/2.)/resolution,y_map_centre-(fit_area/4.)/resolution,y_map_centre,y_map_centre+(fit_area/4.)/resolution,y_map_centre+(fit_area/2.)/resolution]
            plt.xticks(x_ticks)
            plt.yticks(y_ticks)
            # Modify the tick values from pixels to arcmin
            ax.set_xticklabels([str(int((_tick-x_map_centre)*resolution))+"'" for _tick in x_ticks])
            ax.set_yticklabels([str(int((_tick-y_map_centre)*resolution))+"'" for _tick in y_ticks])

            # Plot the residual
            ax_res = plt.subplot(2,3,_i+4)
            residual = np.asarray(current_map) - g(x,y)
            im_res = plt.imshow(residual,interpolation='none',cmap=cmap_res,vmin=vmin_res,vmax=vmax_res)
            plt.xlim([xmin,xmax])
            plt.ylim([ymin,ymax])
            plt.title(_freq+' residuals')
            plt.xticks(x_ticks)
            plt.yticks(y_ticks)
            ax_res.set_xticklabels([str(int((_tick-x_map_centre)*resolution))+"'" for _tick in x_ticks])
            ax_res.set_yticklabels([str(int((_tick-y_map_centre)*resolution))+"'" for _tick in y_ticks])

        if skip_page == 3: 
            offset_dict[_source].pop(obsID, None)
            continue

        # Add page title and colourbars
        plt.suptitle(ob_time.GetFileFormatString()+' '+_source+' obsID '+obsID,fontsize=25)
        fig.subplots_adjust(right=0.85)
        top_ax_pos = ax.get_position()
        bot_ax_pos = ax_res.get_position()
        cbar_ax = fig.add_axes([0.9, top_ax_pos.y0, 0.03, top_ax_pos.height])
        cbar_ax_res = fig.add_axes([0.9, bot_ax_pos.y0, 0.03, top_ax_pos.height])
        cbar = fig.colorbar(im, label='Tcmb',cax=cbar_ax)
        cbar_res = fig.colorbar(im_res, label='Tcmb',cax=cbar_ax_res)
        pdf.savefig()

    pdf.close() 

pkl.dump(offset_dict,open(pkl_savefile,'wr'))
pkl.dump(benchpos_dict,open(output_dir+'benchpos.pkl','wr'))
