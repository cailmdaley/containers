# Sam Guns, last edit 21 December 2018
#
# This script fits a 2D Gaussian beam profile (amplitude, center, width, rotation, offset) to the brightest pixel in every map frame in a g3 file, and writes the resulting fit parameters and covariance matrix to an output g3 file.
#
# Usage:
#
# fit_beamprofile.py inputs (-r) (-x) (-o)
#
# Inputs can be one or more g3 files or folders containing g3 files. Fit information is written into new keys. With the -r option, subfolders will also be opened and searched for g3 files. The original file will be overwritten (all the original data is saved, but new keys are added) only when using the -x option. 
# Options:
#
# -r : search through subfolders for g3 files
# -x : write fit parameter keys to original g3 file (keeps all data)
# -o : overwrite fit parameters if already present
# Standard output format is {original file}_fitdata.g3 unless -x is used. 
#
# Important error: script will crash if brightest pixel is less than 4 arcmin from map edge. 

import sys, os, glob, copy, re, datetime, csv
import numpy as np
import cPickle as pickle
from spt3g import core, mapmaker
import argparse as ap
import scipy as sc

P = ap.ArgumentParser(description='Single bolometer maps with boresight pointing',
              formatter_class=ap.ArgumentDefaultsHelpFormatter)

P.add_argument('inputs', action='store', nargs='+', default=[],
           help='Input files')
P.add_argument('-r', '--recursive', action='store_true', default=False, help='search through all subfolders')
P.add_argument('-x', '--replace', action='store_true', default=False, help='replace original file with analysis output (keeps all data)')
P.add_argument('-o', '--overwrite', action='store_true', default=False, help='overwrite previous fit information, if present (i.e., retry fitting).')
args = P.parse_args()

input_files_and_folders = []
input_files = []

for f in args.inputs:
    if '*' in f:
        input_files_and_folders.extend(glob.glob(f))
    else:
        input_files_and_folders.append(f)
for f in input_files_and_folders:
    if os.path.isfile(f):
        input_files.append(f)
    else:
        f = f if f[-1] == '/' else f + '/'
        if not args.recursive:
            input_files.extend(glob.glob(f+'*.g3'))
        else:
            for folder,sub,files in os.walk(f):
                for filename in files:
                    input_files.append(os.path.join(folder, filename))
input_files = list(set([f for f in input_files if f[-3:] == ".g3"]))

def limit_angle(theta):
    return theta % (2.*np.pi)

def window(T, wsize):
    xmin = (T.shape[1] - wsize) // 2
    ymin = (T.shape[0] - wsize) // 2
    xmax, ymax = xmin + wsize, ymin + wsize
    return (xmin,ymin), np.array(T, copy=True)[ymin:ymax,xmin:xmax]

def px_to_arcmin(x, res):
    return x*res/core.G3Units.arcmin

def px_to_deg(x, res):
    return x*res/core.G3Units.degree

def noise_annulus(Tmap, rmin, rmax):
    #todo: implement custom center
    cx,cy = Tmap.shape[1]//2, Tmap.shape[0]//2
    narr = []
    for (y,x),val in np.ndenumerate(Tmap):
        d = np.sqrt((y-cy)**2 + (x-cx)**2)
        if d >= rmin and d <= rmax:
            narr.append(val)
    return np.nanmean(narr), np.nanstd(narr)

def twoD_Gaussian(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    x, y = xy
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp(-(a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
    return g.ravel()

def gaussfit(Tmap, guess):
    x,y = np.arange(Tmap.shape[1]), np.arange(Tmap.shape[0])
    xv, yv = np.meshgrid(x,y)
    data = np.ravel(Tmap)
    pred_params, uncert_cov = sc.optimize.curve_fit(twoD_Gaussian, (xv, yv), data, p0=guess, maxfev=10000)
    #data_fitted = twoD_Gaussian((xv, yv), *pred_params)
    errors = np.sqrt(np.diag(uncert_cov))
    return pred_params, errors, uncert_cov

def Analyze(frame):
    bands = ['90GHz','150GHz','220GHz']
    if 'Fit attempted' in frame and not args.overwrite:
        return frame

    if 'Fit attempted' in frame and args.overwrite:
        keys = ['Fit attempted', 'Fit success', 'Parameter Description', 'Parameters', 'Errors', 'Covariance matrix', 'Wafer', 'Band', 'Centroid RA', 'Centroid Dec', 'Noise mean', 'Noise std'] 
        for k in keys:
            frame.pop(k, 0)

    fitted = False        
    if 'Id' in frame and 'T' in frame and 'Wunpol' in frame:
        map_id = frame['Id']
        wafer = map_id[:4].upper() if map_id[0].upper() == 'W' else 'All'
        band = 'All'
        for b in bands:
            if b in map_id:
                band = b
        T = mapmaker.mapmakerutils.remove_weight_t(frame['T'], frame['Wunpol'])
        px_size_x, px_size_y = T.shape[1],T.shape[0]
        res = frame['T'].res
        window_size = int(8.*core.G3Units.arcmin/res)
        (x0,y0), T = window(T, window_size)
        cy, cx = np.unravel_index(np.nanargmax(T), T.shape)
        amp0 = np.nanmax(T)
        noise_mean, noise_err = noise_annulus(T,4.*core.G3Units.arcmin/res,10.*core.G3Units.arcmin/res)
        guess = [amp0, cx, cy, 1.*core.G3Units.arcmin/res, 1.*core.G3Units.arcmin/res, 0., noise_mean]
        frame['Fit attempted'] = True
        fitted = True
        if np.all(T==0):
            fitted = False
        try:
            params, errors, covariance_matrix = gaussfit(T, guess)
        except (RuntimeError, ValueError):
            fitted = False
            print "fit failed in map %s" %(map_id)

    if fitted:
        params[1] = params[1] + x0
        params[2] = params[2] + y0
        params[5], errors[5] = limit_angle(params[5]), limit_angle(errors[5])
        (ra,dec) = frame['T'].pixel_to_angle(int(params[1]),int(params[2]))

        frame['Parameter Description'] = "[amplitude (map units), centroid x (px), centroid y (px), sigma x (px), sigma y (px), rotation angle (rad), offset (map units)]"
        frame['Parameters'] = core.G3VectorDouble(params)
        frame['Errors'] = core.G3VectorDouble(errors)
        frame['Covariance matrix'] = core.G3VectorDouble(covariance_matrix)
        frame['Wafer'] = wafer
        frame['Band'] = band
        frame['Centroid RA'] = ra
        frame['Centroid Dec'] = dec
        frame['Noise mean'] = noise_mean
        frame['Noise std'] = noise_err

    frame['Fit success'] = fitted
    return frame

for obs_path in input_files:
    p = core.G3Pipeline()
    p.Add(core.G3Reader, filename=obs_path)
    output = obs_path[:-3] + '_fitdata.g3'
    p.Add(Analyze)
    p.Add(core.G3Writer, filename=output)
    p.Run()

    if args.replace:
        cmd = "mv %s %s" %(output, obs_path)
        os.system(cmd)



                     


    

    






