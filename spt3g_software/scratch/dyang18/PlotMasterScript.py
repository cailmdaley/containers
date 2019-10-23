'''
Overview: Use Condor_submit function to generate and submit a script that make plots for a specific
observation, eventually I want this script to generate plots for all observations of a specific type
'''
from astropy.time import Time
import argparse as ap
from spt3g import core, std_processing, mapmaker, dfmux, calibration, xtalk, todfilter, coordinateutils,cluster
import scipy.stats
import os, numpy
from spt3g.std_processing import obsid_to_g3time
from glob import glob

wafer_lst=['W136' 'W139' 'W142' 'W147' 'W148' 'W152' 'W153' 'W157' 'W158' 'W162']

P = ap.ArgumentParser(description='Making Plots for a observation',
              formatter_class=ap.ArgumentDefaultsHelpFormatter)
#P.add_argument('input_files', action='store', nargs='+', default=[],
#          help='Input files')

P.add_argument('source', action='store',
           default=None, help='name of source')

P.add_argument('-o','--observation', action='store', default=None, help='observation ID')
P.add_argument('-d','--daterange', action='store', default=None, help='date range for timeseries plot like 2017-08')


arg = P.parse_args()


bands=['90.0','150.0','220.0']
script=arg.source+'PlotCondor.py'
args=[ ]
log_root='/scratch/dyang18'
output_root='/spt/user/dyang18/output'
input_files=['/spt/user/dyang18/'+arg.source+'.tar.gz', '/spt/user/dyang18/'+arg.source+ 'BolometerProperties.tar.gz']
if arg.daterange is None or arg.observation is None:
    output_files=[arg.source+".tar.gz"]
else:
    output_files=[arg.source+".tar.gz"]
#bundle these things together in order to transfer one giant thing
grid_proxy='/home/dyang18/user.cert'

cluster.condor_submit(script=script, args=args, caller='python', jobname=None,
                  log_root=log_root, output_root=output_root, input_files=input_files,
                  output_files=output_files, grid_proxy=grid_proxy, aux_input_files=[],
                  aux_output_files=[], spt3g_env=True, python3=False,
                  user_code='', requirements=None, request_cpus=1,
                  request_memory=6*core.G3Units.GB,
                  request_disk=10*core.G3Units.GB, queue=1,
                  when_to_transfer_output='ON_EXIT', create_only=False,
                  globus_uri='gsiftp://gridftp.grid.uchicago.edu:2811/cephfs')
