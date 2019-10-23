import numpy as np
import os
import re
from spt3g import core

"""
This script takes the (RA, Dec) of the N brightest sources and makes 
a text file for each listing the observations frames that contain that
source.
"""

npoints = 10

coords = np.zeros([npoints, 2])

ps_file = os.path.join(os.getenv('SPT3G_SOFTWARE_PATH'), 'beams',
                        'mapmaking_scripts',
                        '1500d_3band_10sigma_ptsrc.txt')

ps_list = np.loadtxt(ps_file)
file_list = {}
for p in range(npoints):
    coords[p] = [ps_list[p][1], ps_list[p][2]]
    file_list[str(p)] = []


fields = ['ra0hdec-44.75', 'ra0hdec-52.25', 'ra0hdec-59.75', 'ra0hdec-67.25']
data_path = os.path.join('/spt', 'data', 'bolodata', '{samplerate}',
                         '{field}')

dec_pad = 1.5 # degrees extra from boresight on either side to include
# focal plane is ~+/1 deg in elevation so this is conservative

def get_dec(fl, positions=['first', 'last']):
    if positions == ['first']:
        for fr in fl:
            if fr.type == core.G3FrameType.Scan:
                dec0 = fr['OnlineBoresightDec'][0]/core.G3Units.degrees
                return dec0
    elif positions == ['last']:
        for fr in fl:
            if fr.type == core.G3FrameType.Scan:
                dec0 = fr['OnlineBoresightDec'][-1]/core.G3Units.degrees
        return dec0
    else:
        dec0 = None
        for fr in fl:
            if fr.type == core.G3FrameType.Scan:
                if dec0 is None:
                    # first dec
                    dec0 = fr['OnlineBoresightDec'][0]/core.G3Units.degrees
                    dec1 = fr['OnlineBoresightDec'][-1]/core.G3Units.degrees
                else:
                    dec1 = fr['OnlineBoresightDec'][-1]/core.G3Units.degrees
                    # eventually last dec
        return dec0, dec1

def is_dec_in_obs(dec_start, dec_end, file_name, file_list):
    # include padding since dec is boresight only
    dec_start += dec_pad
    dec_end -= dec_pad
    #print('searching in decs',dec_start, dec_end)
    for pi, p in enumerate(coords):
        if np.logical_and(dec_start > p[1], dec_end < p[1]):
            file_list[str(pi)].append(file_name)
            print('Adding {} to {}'.format(file_name, pi))

# Now, for each observation frame, check at what dec it and the next starts.
# If a coord is within that dec range, add the obs to that file list
for field in sorted(fields):
    print(field)
    map_dir = data_path.format(samplerate='fullrate', field=field)
    for obs in sorted(os.listdir(map_dir)):
        if int(obs[:2]) < 47 or int(obs[:2]) > 57:
            # Daniel's subset is obsids 47* to 57*
            continue
        print(obs)
        g3_files = os.listdir(os.path.join(map_dir, obs))
        # throw out non-observation files
        g3_files = list(filter(lambda x: re.match("[0-9][0-9][0-9][0-9].g3", x),
                               g3_files))
        # sort just in case
        g3_files = sorted(g3_files)
        first = g3_files[0]
        last = g3_files[-1]

        for i, f in enumerate(g3_files):
            print(f)
            fl = core.G3File(os.path.join(map_dir, obs, f))
            if f == last or first == last:
                dec_start, dec_end = get_dec(fl, positions=['first', 'last'])
                if first != last:
                    # check previous file, as usual
                    is_dec_in_obs(dec0, dec_start, file_name, file_list)
                    
                # now, search within the file
                file_name = os.path.join(data_path.format(samplerate='fullrate',
                                                          field=field),
                                         obs, f)
                is_dec_in_obs(dec_start, dec_end, file_name, file_list)
            elif f == first:
                # Set up starting dec and file name
                file_name = os.path.join(data_path.format(samplerate='fullrate',
                                                          field=field),
                                         obs, f)
                dec0 = get_dec(fl, positions=['first'])
            else:
                # Get starting dec. Search between this and previous file's
                # starting dec. If a source is found, previous file is added
                # to that source's list.
                dec1 = get_dec(fl, positions=['first'])
                is_dec_in_obs(dec0, dec1, file_name, file_list)
                # Now this file's starting dec/file are set to be first.
                dec0 = dec1
                file_name = os.path.join(data_path.format(samplerate='fullrate',
                                                          field=field),
                                         obs, f)
        
        
    np.savez('/spt/user/agambrel/configs/point_source_obs_{}.npz'.format(field),
             **file_list)
