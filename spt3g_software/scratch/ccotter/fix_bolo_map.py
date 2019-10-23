'''
This script is used to identify mismatched bolometers given a set of
observations. First it finds a set of bolometers that appear to be well-matched.
Then it calculates the average position of the bolometer from the observations
and subtracts this from the expected bolometer position. Making a distribution
of these deltas allows us to have a confidence score for how likely a
bolometer is mismatched. If the delta of a bolometer is greater than three
standard deviations away from its expected position, we say it's mismatched.


TODO:
Need to change comparison so that it doesn't depend on
number of observations

Will need a better way to define well-matched bolos

Could add pixel consistency check to make sure all bolos
in a pixel are next to each other

Should check if wafers are causing bumps in delx/dely plot
'''


from spt3g import core, calibration
import matplotlib.pyplot as plt
import numpy as np
import copy
import csv
from optparse import OptionParser
from os import listdir

# setup command line options
usage = 'Usage: fix_bolo.py -t TARGET [options]'

parser = OptionParser(usage=usage)
parser.add_option('-t', '--target', action='store', type='string', dest='target', metavar='TARGET', help='The target of the observations; e.g. "RCW38-pixelraster"')
parser.add_option('-b', '--observations', action='store', type='string', dest='obs', metavar='FILE', help='File containing a list of newline separated observation times to use. Defaults to all observations in nominal target directory.')
parser.add_option('-m', '--map', action='store', type='string', dest='mapping', metavar='FILE', help='File containing the mapping of bolometers to LC board (can be downloaded from bitbucket). Doesn\'t output LC or mapping changes if not provided.')
parser.add_option('-o', '--output', action='store', type='string', dest='out', metavar='FILE', default='bad_bolos.txt', help='Output for bad bolometers')
parser.add_option('-p', '--plots', action='store_true', dest='plot', default=False, help='Generate plots of delta x and delta y for reference bolometers')

(options, args) = parser.parse_args()

if (options.target == None):
    print(usage)
    exit()

# paths to nominal and offline data
nom_path = '/spt/data/bolodata/downsampled/' + options.target + '/'
off_path = '/spt/user/production/calibration/' + options.target + '/'

# chosen observations with decent data
good_times = []
if (options.obs != None):
    with open(options.obs, 'r') as f:
        for line in f:
            good_times.append(line.strip())
else:
    good_times = listdir(nom_path)


# dictionary containing well-matched bolometers and their associated offsets
ref_bolos = dict()
# use the first observation nominal data to populate the dictionary with all
# possible bolometers and trim it after going through all data
temp = [frame for frame in core.G3File(nom_path + good_times[0] + '/nominal_online_cal.g3')]
for key, item in temp[0]['NominalBolometerProperties'].items():
    ref_bolos[key] = dict(nom_x=item.x_offset, nom_y=item.y_offset)

# fill in all bolometers that have offsets for each observation
for direc in good_times:
    off_data = [frame for frame in core.G3File(off_path + direc + '.g3')]
    off_bolo_x = off_data[0]['PointingOffsetX']
    off_bolo_y = off_data[0]['PointingOffsetY']

    for key, item in off_data[0]['PointingOffsetX'].items():
        ref_bolos[key][direc] = dict(x=item)
    for key, item in off_data[0]['PointingOffsetY'].items():
        ref_bolos[key][direc]['y'] = item

# count number of times bolos appear in an observation
# NOTE: 12 is the max length for 10 obs because of nominal offsets
'''
count10 = 0
count9 = 0
count8 = 0
count7 = 0
count6 = 0
ref_copy = dict(ref_bolos)
for key, item in ref_copy.items():
    if (len(item) == 10):
        count10 += 1
    elif (len(item) == 9):
        count9+= 1
    elif (len(item) == 8):
        count8 += 1
    elif (len(item) == 7):
        count7 += 1
    elif (len(item) == 6):
        count6 += 1
print(count10)
print(count9)
print(count8)
print(count7)
print(count6)
'''

# make a copy of all bolometer data to be used later
all_bolos = copy.deepcopy(ref_bolos)
# add a field to keep track of number of observations
for key, item in all_bolos.items():
    item['num_obs'] = len(item) - 2 # - 2 for nominal fields

# read mapping file and store information in all_bolos
if (options.mapping != None):
    with open(options.mapping, newline='') as csvfile:
        map_reader = csv.reader(csvfile, delimiter='\t')
        for row in map_reader:
            if row[1] in all_bolos.keys():
                all_bolos[row[1]]['lc'] = row[0]

# clean up the reference bolometers
ref_copy = copy.deepcopy(ref_bolos)
for key, item in ref_copy.items():
    # throw out any bolometers not in all observations
    if (len(item) != 12):
        ref_bolos.pop(key)
        continue
    # remove any bolometers that have offsets separated by more than 0.001
    x = [obs['x'] for key, obs in item.items() if isinstance(obs, dict)]
    if max(x) - min(x) > 0.001:
        ref_bolos.pop(key)
        continue
    y = [obs['y'] for key, obs in item.items() if isinstance(obs, dict)]
    if max(y) - min(y) > 0.001:
        ref_bolos.pop(key)

# make histograms of the offsets for each bolometer and overplot the nominal offset
'''
for key, obs in ref_bolos.items():
    plt.axvline(x=obs['nom_x'], c='r')
    x = [item['x'] for key, item in obs.items() if isinstance(item, dict)]
    n, bins, patches = plt.hist(x, bins=10)
    plt.axis([min(x), max(x), 0, 10])
    plt.show()
    plt.close()
'''

# calculate mean offset of each reference bolometer and subtract from nominal offset
delx = []
dely = []
for key, item in ref_bolos.items():
    x = [obs['x'] for key, obs in item.items() if isinstance(obs, dict)]
    xhat = np.mean(x)
    delx.append(item['nom_x'] - xhat)
    y = [obs['y'] for key, obs in item.items() if isinstance(obs, dict)]
    yhat = np.mean(y)
    #NOTE: y offsets are flipped relative to nominal so multiply by -1
    dely.append(item['nom_y'] + yhat)


# calculate mean and std of the deltas to get a confidence of a bolometer
# being mapped correctly
mean_delx = np.mean(delx)
std_delx = np.std(delx)
mean_dely = np.mean(dely)
std_dely = np.std(dely)


if (options.plot):
    # plot the reference bolometer deltas along with the mean and std
    n, bins, patches = plt.hist(delx, bins=100)
    plt.axvline(x=mean_delx, c='r', label='mean')
    plt.axvline(x=mean_delx + 3 * std_delx, c='g', label='3 sigma')
    plt.axvline(x=mean_delx - 3 * std_delx, c='g')
    plt.axis([min(delx), max(delx), 0, max(n)])
    plt.xlabel('Separation of mean offset from nominal')
    plt.ylabel('Number of reference bolometers')
    plt.title('Distribution of x offsets')
    plt.legend()
    plt.show()
    plt.close()

    n, bins, patches = plt.hist(dely, bins=100)
    plt.axvline(x=mean_dely, c='r', label='mean')
    plt.axvline(x=mean_dely + 3 * std_dely, c='g', label='3 sigma')
    plt.axvline(x=mean_dely - 3 * std_dely, c='g')
    plt.axis([min(dely), max(dely), 0, max(n)])
    plt.xlabel('Separation of mean offset from nominal')
    plt.ylabel('Number of reference bolometers')
    plt.title('Distribution of y offsets')
    plt.legend()
    plt.show()
    plt.close()

# calculate deltas for all bolometers
for key, item in all_bolos.items():
    if (item['num_obs'] > 0):
        x = [obs['x'] for key, obs in item.items() if isinstance(obs, dict)]
        xhat = np.mean(x)
        item['delx'] = item['nom_x'] - xhat
        y = [obs['y'] for key, obs in item.items() if isinstance(obs, dict)]
        yhat = np.mean(y)
        #NOTE: y offsets are flipped relative to nominal so multiply by -1
        item['dely'] = item['nom_y'] + yhat



# write all bad bolos with their offset
with open(options.out, 'w') as output:
    for key, item in all_bolos.items():
        if (item['num_obs'] >= 8):
            if (item['delx'] > mean_delx + 3 * std_delx
                    or item['delx'] < mean_delx - 3 * std_delx
                    or item['dely'] > mean_dely + 3 * std_dely
                    or item['dely'] < mean_dely - 3 * std_dely):
                sep_x = item['delx'] - mean_delx
                sep_y = item['dely'] - mean_dely
                if 'lc' in item.keys():
                    output.write(key + '\t' + item['lc'] + '\t' + str(sep_x) + '\t' + str(sep_y) + '\n')
                else:
                    output.write('Bad bolometer without mapping: ' + key + '\n')












