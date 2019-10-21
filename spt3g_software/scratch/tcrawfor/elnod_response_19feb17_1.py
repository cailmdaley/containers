import numpy as np
import scipy
from scipy import ndimage
import pickle
import argparse as ap
from spt3g import core, std_processing, gcp
from spt3g.util import tctools
import os

file1 = '/spt/data/bolodata/downsampled/elnod/4284593/0000.g3'
cfile1 = '/spt/data/bolodata/downsampled/elnod/4284593/nominal_online_cal.g3'

f2 = core.G3File(cfile1)
bp = f2.next()['NominalBolometerProperties']

bdict0 = {}
el_raw_0 = []
az_raw_0 = []
f1 = core.G3File(file1)
for frame in f1:
    if frame.type is core.G3FrameType.Scan:
        if 'Turnaround' not in frame:
            names = frame['RawTimestreams_I'].keys()
            if len(bdict0) == 0:
                for name in names:
                    bdict0[name] = []
                    bdict0[name].append(frame['RawTimestreams_I'][name])
            else:
                for name in names:
                    bdict0[name].append(frame['RawTimestreams_I'][name])
            el_raw_0.append(frame['RawBoresightEl'])
            az_raw_0.append(frame['RawBoresightAz'])
        
bdict = {}
for name in names:
    bdict[name] = tctools.list_to_array(bdict0[name])
el_raw = tctools.list_to_array(el_raw_0)
az_raw = tctools.list_to_array(az_raw_0)

# kludge
nmax = 2980
el_raw = el_raw[0:nmax]
az_raw = az_raw[0:nmax]
for name in names:
    bdict[name] = bdict[name][0:nmax]

template = 1./np.sin(el_raw)
template -= np.mean(template)

elnod_dict = {}
for name in bdict.keys():
    dirtemp = {}
    bdtemp = np.asarray(bdict[name])
    if np.sum(bdtemp) == 0.:
        continue
    
    slope, intercept, r_value, p_value, std_err = \
        scipy.stats.linregress(template,bdtemp)
    dirtemp['elnod_response'] = slope
    dirtemp['physical_name'] = bp[name].physical_name
    model = template*slope + intercept
    resids = bdtemp - model
    dirtemp['elnod_sn'] = 0.
    resids = bdtemp - model
    resid_rms = np.std(resids)
    if resid_rms > 0.:
        dirtemp['elnod_sn'] = np.sqrt(np.sum((template*slope/resid_rms)**2))
#    if np.abs(std_err) > 0.:
#        dirtemp['elnod_sn'] = np.abs(slope/std_err)
    elnod_dict[name] = dirtemp

elnod_response = np.asarray([elnod_dict[name]['elnod_response'] for name in elnod_dict.keys()])
elnod_sn = np.asarray([elnod_dict[name]['elnod_sn'] for name in elnod_dict.keys()])
