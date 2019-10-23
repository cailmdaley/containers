from spt3g.pointing import pointing_tools as pt
from spt3g import core, gcp, std_processing
import numpy as np
import os
import socket # we need this to get the hostname


#What was the last date to have az tilt fits for (redo the last few to catch
#dates where arc data was missing for intermittent days)?
#start_az_date = core.G3Time('20-Jan-2017:00:00:00')

start_az_date = pt.find_newest_tilts(tilts_dir='/home/jhenning/data/pointing/tilt_fits/')

#Generate new lists of cleaned az tilt start and stop times.
print 'Obtaining az tilt fits for observations after:', start_az_date.Summary()
startfile, stopfile = pt.get_az_tilt_start_files(start_time=start_az_date.GetFileFormatString(), 
                                                 runlog_location='/spt_data/gcplogs/')

mjd_avg, tilt_lat, tilt_ha, tilt_mag, tilt_angle = pt.az_tilt(startfile,stopfile, plotting=1, 
                                                    skip_to=start_az_date, 
                                                    az_tilt_fit_dir=os.path.join(os.environ['HOME'],'data/pointing','tilt_fits'), 
                                                    arcdir='/spt_data/arc/', 
                                                    skip_day_fraction=0.0)

