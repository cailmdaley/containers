import numpy as np
import scipy
from scipy import ndimage
import pickle
from spt3g import core, std_processing, gcp
import os
import glob

# get files
fields = ['ra0hdec-44.75','ra0hdec-52.25','ra0hdec-59.75','ra0hdec-67.25']
bands = ['90','150','220']
file_dict = {}
for band in bands:
    file_dict[band] = {}
    for field in fields:
        files = glob.glob('/spt/data/onlinemaps/'+field+'/*_'+band+'GHz*.g3.gz')
        files.sort()
        file_dict[band][field] = files

#ny = 12000
#dec_center = -57.5
#reso_arcmin = 0.25
ny = 1500
dec_center = -57.5
reso_arcmin = 2.0

# declination vector
dec = np.arange(ny)*reso_arcmin/60.
dec -= np.mean(dec)
dec += dec_center
wvec0 = np.cos(dec_center*np.pi/180.)/np.cos(dec*np.pi/180.)

# define field weights
wvec_dict = {}
nhann = np.int(np.round(2.*60./reso_arcmin))
htemp = np.hanning(2*nhann)
hann_left = htemp[0:nhann]
hann_right = htemp[nhann:]
nobs_dict = {}
for field in fields:
    wvec_dict[field] = wvec0.copy()
    dec_mid = -np.float(field.split('-')[-1])
    dec_min = dec_mid - 4.75
    dec_max = dec_mid + 4.75
    ind1 = np.max(np.where(dec <= dec_min)[0])
    ind2 = np.min(np.where(dec >= dec_max)[0])
    wvec_dict[field][0:ind1] = 0.
    wvec_dict[field][ind1:ind1+nhann] *= hann_left
    wvec_dict[field][ind2:] = 0.
    wvec_dict[field][ind2-nhann:ind2] *= hann_right
    nobs_dict[field] = {}

# now make slice by adding one copy of field weights for every obs
slice_dict = {}
for band in bands:
    slice_dict[band] = {}
    for field in fields:
        nobs_dict[field][band] = {}
month_names = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct']
m2do = [4,5,6]
for mind in m2do:
    month = month_names[mind-1]
    next_month = month_names[mind]
    start_time = core.G3Time('01-'+month+'-2019:00:00:00')
    start_obsid = std_processing.time_to_obsid(start_time)
    end_time = core.G3Time('01-'+next_month+'-2019:00:00:00')
    end_obsid = std_processing.time_to_obsid(end_time)
    for band in bands:
        slice_dict[band][month] = np.zeros(len(dec))
        for field in fields:
            nobs_dict[field][band][month] = 0
            for file1 in file_dict[band][field]:
                obsid1 = np.int(file1.split('/')[-1].split('_')[0])
                if obsid1 > start_obsid and obsid1 < end_obsid:
                    slice_dict[band][month] += wvec_dict[field]
                    nobs_dict[field][band][month] += 1

file_dict_2 = {}
nobs_dict_2 = {}
for field in fields:
    files_2 = glob.glob('/spt/data/bolodata/downsampled/'+field+'/*/0000.g3')
    files_2.sort()
    file_dict_2[field] = files_2
    nobs_dict_2[field] = {}
for mind in m2do:
    month = month_names[mind-1]
    next_month = month_names[mind]
    start_time = core.G3Time('01-'+month+'-2019:00:00:00')
    start_obsid = std_processing.time_to_obsid(start_time)
    end_time = core.G3Time('01-'+next_month+'-2019:00:00:00')
    end_obsid = std_processing.time_to_obsid(end_time)
    for field in fields:
        nobs_dict_2[field][month] = 0
        for file1 in file_dict_2[field]:
            obsid1 = np.int(file1.split('/')[-2])
            if obsid1 > start_obsid and obsid1 < end_obsid:
                nobs_dict_2[field][month] += 1
