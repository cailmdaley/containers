from sptpol_software.analysis.pointing import pointing_tools_3g as pt
import sptpol_software.analysis.pointing.consolidate_source_scan as css
#import sptpol_software.analysis.pointing.fit_hii_for_pointing as fh2
import sptpol_software.util.time as time
import numpy as np
import pylab as py
import os
import socket # we need this to get the hostname

focus_position = None #Let fit_hii_for_pointing figure out whether
                      #a source observation is at a good focus_position
                      #that matches the CMB field focus.

config_dir = os.environ.data['SPTPOL_SOFTWARE']+'/config_files/'
pointing_dir = os.environ.data['SPTPOL_SOFTWARE']+'/analysis/pointing/'
pointing_dir_source_scan = '/data/jhenning/pointing/'
old_tilt_config = config_dir+'sptpol_offline_tilts_20141202_000000.config'
old_hii_config = config_dir+'sptpol_offline_hii_params_20150810_201500.config'
source_summary_file = pointing_dir_source_scan+'source_scan_structure_20150218_201500.pkl'

############################Az Tilt Pointing Config Generation#########################
#What was the last date to have az tilt fits for (redo the last few to catch
#dates where arc data was missing for intermittent days)?
mjd,a2,a3,tilt_ang = pt.read_tilt_config(old_tilt_config, file_type='config')
start_az_date = str(time.SptDatetime(mjd[-1]-1))

#Generate new lists of cleaned az tilt start and stop times.
print 'Obtaining az tilt fits for observations after:', start_az_date
startfile, stopfile = pt.get_az_tilt_start_files(start_time=time.SptDatetime(start_az_date), runlog_location='/buffer/gcplogs/')
tilt_fit_file = pt.az_tilt(startfile,stopfile, plotting=1, skip_to=start_az_date, az_tilt_fit_dir='az_tilt_fits_corrected', arcdir='/buffer/arc/', skip_day_fraction=0.0)

#Write a new config file, appending to the data from an old config file.
new_tilt_config = pt.write_tilt_config(tilt_fit_file=tilt_fit_file,config_file=old_tilt_config)

#Write the new tilt file to the config_dir.
os.system('\cp -f '+new_tilt_config+' '+old_tilt_config)
####################################################################################


#######################HII Region Pointing Config Generation#######################
#Append new source observations to source file.

# arg.  We can't do this in the same way on anal and cloud
if socket.gethostname() == 'ana':
    t = css.consolidate_source_scan(start_date='01-Jan-2015:00:00:00', 
                                    stop_date=None, 
                                    fitsdir='/data/sptdat/auxdata/', 
                                    sources_to_use=['rcw38','mat5a','cena'],
                                    filename=source_summary_file, 
                                    savename=None, 
                                    dosave=True, 
                                    doall=False)
else:
    t = css.consolidate_source_scan(start_date='01-Jan-2015:00:00:00', 
                                    stop_date=None, 
                                    fitsdir='/data/sptdat/auxdata/', 
                                    sources_to_use=['rcw38','mat5a','cena'],
                                    filename=source_summary_file, 
                                    savename=None, 
                                    dosave=True, 
                                    doall=False)


#What was the last date to have hii fits for?
mjd,a0,a1,a4,a5,a6,az0 = pt.read_hii_config(old_hii_config)
hii_fit_start = str(time.SptDatetime(mjd[-1]))

#Do the new HII region fits.
print 'Obtaining HII fits for observations after:', hii_fit_start
hii_savename = fh2.fit_hii_for_pointing(source_summary_file, dodebug=False, \
                             test_mode=False,skip_cal=True,skip_elnod=True,\
                             skip_ampVsNoise=True, use_this_source='rcw38', \
                             use_this_a5=None, use_this_a6=None, \
                             fit_a5=False, \
                             datestart=hii_fit_start,\
                             focus_position_input=focus_position)

#Write a new HII pointing config file, appending to data in the most current config file.
new_hii_config = pt.write_hii_config(hii_savename, config_file=old_hii_config,focus=focus_position)

#Write the new HII file to the config_dir.
os.system('\cp -f '+new_hii_config+' '+old_hii_config)
########################################################################################


