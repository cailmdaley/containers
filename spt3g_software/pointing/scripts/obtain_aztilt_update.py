from spt3g.pointing import pointing_tools as pt
from spt3g import core
import os, glob
import platform

if platform.node() == 'anal.spt':
    tilts_dir = '/poleanalysis/sptdaq/azTilts/'
    gcplogs_dir = '/spt_data/gcplogs/'
    arcdir = '/spt_data/arc/'
else:
    tilts_dir = os.path.join(os.getenv('SPT3G_SOFTWARE_PATH'), 'pointing/aztilts/')
    gcplogs_dir = '/spt/data/gcplogs/'
    arcdir = '/spt/data/arc/'

#What was the last date to have az tilt fits for (redo the last few to catch
#dates where arc data was missing for intermittent days)?
start_az_date, obs_id = pt.find_newest_tilts(tilts_dir=tilts_dir)

#Generate new lists of cleaned az tilt start and stop times.
print('Obtaining az tilt fits for observations after:', start_az_date.Summary())
tilt_times = pt.find_tilt_times(start_time=start_az_date, runlog_location=gcplogs_dir)

mjd_avg, tilt_lat, tilt_ha, tilt_mag, tilt_angle = pt.az_tilt(tilt_times,
                                                    plotting=1, 
                                                    skip_to=start_az_date, 
                                                    az_tilt_fit_dir=tilts_dir, 
                                                    arcdir=arcdir, 
                                                    skip_day_fraction=0.0)

#Find the newest tilt now that we potentially analyzed a few more.
new_az_date, obs_id = pt.find_newest_tilts(tilts_dir=tilts_dir)

new_command = pt.get_online_tilt_correction(obs_id, tilts_dir=tilts_dir)
