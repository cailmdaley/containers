from spt3g import core, calibration
import numpy as np
import scipy as sp
import sys
import os, glob
import argparse as ap

P = ap.ArgumentParser(description='Construct set of relative pointing offsets from RCW38 fastpoints',
              formatter_class=ap.ArgumentDefaultsHelpFormatter)

P.add_argument('-o', '--output', action='store', default='offsets.g3',
           help='output filename')
P.add_argument('-n', '--min_number_good_obs', action='store', default = '5',
           help='minimum number of succesful (non-nan or inf) observations per bolometer')
P.add_argument('-a', '--after_obs', action='store', default='40000000', 
           help='only use observations on or after this id')
P.add_argument('-e', '--end_obs', action='store', default='1e8',
           help='only use observations before this id')
P.add_argument('-v', '--verbose', action='store_true', default=False)
args = P.parse_args()

def get_obsid(path):
    return os.path.splitext(os.path.basename(path))[0]

obs_root = '/spt/user/production/calibration/RCW38-pixelraster/'
#obs_lst = [get_obsid(p) for p in glob.glob(obs_root + '*.g3') if (int(get_obsid(p)) >= int(args.after_obs) and int(get_obsid(p)) <= int(args.end_obs))]

good_obs = [
'32611890',
'33387301',
'33480202',
'33677605',
'33924973',
'34028129',
'34413272',
'34882327',
'35421843',
'36047851',
'36169722',
'36282797',
'36464817',
'36593069',
'36724500',
'36856066',
'36888994',
'36985606']

obs_lst = [get_obsid(p) for p in glob.glob(obs_root + '*.g3') if get_obsid(p) in good_obs]
obs_lst = [p for p in obs_lst if (int(p) >= int(args.after_obs) and int(p) <= int(args.end_obs))]

if args.verbose:
    print "found %d observations" %(len(obs_lst))

#get a list of all bolos
example_boloprops = core.G3File('/spt/data/bolodata/downsampled/RCW38-pixelraster/47017623/offline_calibration.g3').next()['BolometerProperties']
all_bolos = example_boloprops.keys()

x_pt = {b: np.array([np.nan]*len(obs_lst)) for b in all_bolos}
y_pt = {b: np.array([np.nan]*len(obs_lst)) for b in all_bolos}

for i,obs in enumerate(obs_lst):
    f = core.G3File(obs_root + obs + '.g3')
    fr = f.next()
    assert(set(fr['PointingOffsetX'].keys()) == set(fr['PointingOffsetY'].keys()))
    bolos = fr['PointingOffsetX'].keys()
    for b in bolos:
        x_pt[b][i] = fr['PointingOffsetX'][b] if not np.isinf(fr['PointingOffsetX'][b]) else np.nan
        y_pt[b][i] = fr['PointingOffsetY'][b] if not np.isinf(fr['PointingOffsetY'][b]) else np.nan

#find and subtract the mean focal plane position of the common pixels in each obs to get rid of boresight jitter between obs
common_pixels = []
for b in all_bolos:
    nanx = np.isnan(x_pt[b])
    nany = np.isnan(y_pt[b])
    if np.logical_or(nanx, nany).sum() == 0:
        common_pixels.append(b)

if args.verbose:
    print "%d common pixels found" %(len(common_pixels))
fiducial = ['015.4.1.3.2028', '005.6.2.3.1860', '005.6.1.3.3104', '005.9.2.1.1722', '005.13.1.4.4774']
common_pixels = fiducial

x_avg = np.zeros(len(obs_lst))
y_avg = np.zeros(len(obs_lst))
for cp in common_pixels:
    x_avg += x_pt[cp]
    y_avg += y_pt[cp]
x_avg /= len(common_pixels)
y_avg /= len(common_pixels)

def pass_cut(bolo, cut):
    cut = cut*core.G3Units.arcmin
    ma_x = np.nanmax(x_pt[bolo])
    mi_x = np.nanmin(x_pt[bolo])
    ma_y = np.nanmax(y_pt[bolo])
    mi_y = np.nanmin(y_pt[bolo])
    mdx = np.abs(ma_x - mi_x)
    mdy = np.abs(ma_y - mi_y)
    return ((mdx < cut) and (mdy < cut))


# take median offset for bolo if passes required # of observations cut, otherwise nan
x_out = core.G3MapDouble()
y_out = core.G3MapDouble()
for b in all_bolos:
    i_obs_x = np.invert(np.isnan(x_pt[b]))
    i_obs_y = np.invert(np.isnan(y_pt[b]))
    
    if (i_obs_x.sum() >= int(args.min_number_good_obs) and i_obs_y.sum() >= int(args.min_number_good_obs) and pass_cut(b, 100.)):  
        x_out[b] = np.nanmedian(x_pt[b] - x_avg) 
        y_out[b] = np.nanmedian(y_pt[b] - y_avg)
    else:
        x_out[b] = np.nan
        y_out[b] = np.nan

# subtract average of new offsets and add average of old offsets 
old_pixel_x_avg = []
old_pixel_y_avg = []
new_pixel_x_avg = []
new_pixel_y_avg = []
for b in all_bolos:
    if not (np.isnan(x_out[b]) or np.isnan(y_out[b]) or 
            np.isnan(example_boloprops[b].x_offset) or np.isnan(example_boloprops[b].y_offset)):
        new_pixel_x_avg.append(x_out[b])
        new_pixel_y_avg.append(y_out[b])
        old_pixel_x_avg.append(example_boloprops[b].x_offset)
        old_pixel_y_avg.append(example_boloprops[b].y_offset)

old_pixel_x_avg = np.mean(old_pixel_x_avg)
old_pixel_y_avg = np.mean(old_pixel_y_avg)
new_pixel_x_avg = np.mean(new_pixel_x_avg)
new_pixel_y_avg = np.mean(new_pixel_y_avg) 

n_success = 0
n_nans = 0
wafer_lst = []

for b in all_bolos:
    x_out[b] += old_pixel_x_avg - new_pixel_x_avg
    y_out[b] += old_pixel_y_avg - new_pixel_y_avg
    if np.isnan(x_out[b]) or np.isnan(y_out[b]):
        n_nans += 1
    else:
        n_success += 1
        wafer_lst.append(example_boloprops[b].wafer_id)
if args.verbose:
    print "total %d bolos\ngood pointing %d bolos\nfailed pointing %d bolos"%(len(all_bolos), n_success, n_nans) 
    wafers = set(wafer_lst)
    for w in wafers:
        print "%s: %d good bolos"%(w, wafer_lst.count(w))

outfr = core.G3Frame(core.G3FrameType.Calibration)

outfr['PointingOffsetX'] = x_out
outfr['PointingOffsetY'] = y_out

wr = core.G3Writer(filename=args.output)
wr.Process(outfr)
