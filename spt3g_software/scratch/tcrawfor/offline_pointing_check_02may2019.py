import numpy as np
import scipy
from scipy import ndimage
import pickle
import argparse as ap
from spt3g import core, std_processing, gcp
from spt3g.mapmaker import mapmakerutils as mm 
import os
import glob
from spt3g.util import gauss_fit

mapfiles = glob.glob('/spt/user/weiquan/notch_filter_test_04-check_high_res_maps/*.g3')
mapfiles.sort()
obsids = [obsid.split('/')[-1] for obsid in glob.glob('/spt/user/ddutcher/ra0hdec-52.25/test_2019_maps/*')]
obsids_52 = []
mapfiles_52 = []
for obsid in obsids:
    for mf in mapfiles:
        if obsid in mf:
            mapfiles_52.append(mf)
            obsids_52.append(obsid)
    
ycenter = 6861
xcenter = 4854

xcen_dict = {}
ycen_dict = {}
map_small_dict = {}

npts_model = 81
xmodel = np.zeros([npts_model,npts_model])
ymodel = np.zeros([npts_model,npts_model])
xtemp = np.arange(npts_model,dtype=int)
for i in np.arange(npts_model):
    xmodel[i,:] = xtemp
    ymodel[i,:] = i
    
for obsid,mfile in zip(obsids_52,mapfiles_52):
    f1 = core.G3File(mfile)
    print(mfile)
    print(obsid)
    for frame in f1:
        if frame.type == core.G3FrameType.Map:
            mframe = frame
        mm.RemoveWeightModule(frame)
        map_unw = np.asarray(mframe['T'])
        map_small = map_unw[ycenter-npts_model/2:ycenter+npts_model/2,xcenter-npts_model/2:xcenter+npts_model/2]
        map_small_dict[obsid] = map_small
        ptemp = np.asarray(gauss_fit.fit2Dgaussian(map_small))
        xcen_dict[obsid] = ptemp[2]
        ycen_dict[obsid] = ptemp[3]


