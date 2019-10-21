import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as pl
import pickle as pk
import argparse as ap
from spt3g import core, std_processing
import sys
sys.path.append('/home/axf295/2019/code/spt3g_software/polcal/scratch/axf295/python/')
import CenA_Map_Utils as CMU
import General_Utils  as GU
import Plotting_Utils as PU



# Usage: Plot_FitParams.py -obs [ObsIDs] -o /path/to/outputdir/ -s CenA -P bool -op /path/to/plot_directory/
P = ap.ArgumentParser(description='Calc Signal to Noise of Single Bolo Source Maps',
              formatter_class=ap.ArgumentDefaultsHelpFormatter)

P.add_argument('-obs','--ObsIDs', action='store',nargs='+', default=[], help='CenA Observation IDs')

P.add_argument('-s', '--source', action='store', default='CenA', help='name of source')


## Static data locations
dataloc = '/big_scratch/axf295/2019/CenA/Analyzed_Data/'


## Set wafers to analyze by adding/removing from this list
wafer_colors_18 = {'W188':'c', 'W174':'b', 'W177':'purple', 'W176':'y', 'W172':'k', 'W180':'orange', 'W181':'brown','W203':'r'}## Removed w2xx #'w203':'r', 'w201':'g',  'w187':'darkblue'

wafer_colors_19 = {'W172':'r', 'W174':'g', 'W176':'c', 'W177':'b', 'W180':'purple', 'W181':'y', 'W188':'k', 'W203':'orange', 'W204':'brown', 'W206':'darkblue'}

## Load inputs
args     = P.parse_args()
source   = args.source
obsids   = args.ObsIDs

## Function to get all obsIDs from 2018-2019
if obsids == []:
    obsids = GU.get_cena_obsids()

failedobs = {}
    
for obs in obsids:
    print(obs)
    nominal_angles = {}
    bolo2band = {}
    bolo_loc  = {}
    pname     = {}
    xy = {}
    AB = {}
    unique_nominal_angles = []
    fitparamloc = dataloc+'%s/'%obs
    output_dir  = fitparamloc+'Polarization_Fit_Plots/'
    try:
        coadd_params   = pk.load(open(fitparamloc+'%s_NomAngCoadd_BestFitParams.pkl'%obs,'rb'))
        singlebolo_params  = pk.load(open(fitparamloc+'%s_Individual_Bolometer_FitParams.pkl'%obs,'rb'))

    except Exception:
        failedobs[obs]='Data File(s) Not Found'
        continue
        
    
    GU.make_directory(output_dir)

    if obs != 'All':
        offline_cal = '/spt/user/production/calibration/calframe/CenA-pixelraster/%s.g3'%obs #'/spt/user/production/calibration/boloproperties/60000000.g3'#
    else:
        ## use most recent boloprops file
        offline_cal = '/spt/user/production/calibration/boloproperties/60000000.g3'
    
    boloprops = {}    
    for frame in core.G3File(offline_cal):
        for bolo in frame['BolometerProperties'].keys():
            boloprops[bolo] =  frame['BolometerProperties'][bolo]
            
    PU.hist_x_minus_y(singlebolo_params,boloprops,savedir=output_dir)
    PU.hist_singlebolo_polangle_residuals(singlebolo_params,boloprops,savedir=output_dir)
    PU.plot_coadd_fit_params(coadd_params,perwafer=True,savedir = output_dir)
    
    
print('%s Failed Obs:'%len(failedobs.keys()))
for obs in failedobs:
    print(obs,failedobs[obs])


