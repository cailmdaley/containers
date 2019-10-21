'''
Script to make individual bolometer maps with boresight pointing for use in
pointing calibration (e.g. for RCW38).
'''
import argparse as ap
from spt3g import core, std_processing, mapmaker, dfmux, calibration, xtalk, coordinateutils
from spt3g.mapmaker.mapmakerutils import MakeMap
import os
import numpy as np

#python map_rcw38.py /spt/data/bolodata/downsampled/RCW38-pixelraster/16508393/offline_calibration.g3 /spt/data/bolodata/downsampled/RCW38-pixelraster/16508393/0000.g3 /spt/data/bolodata/downsampled/RCW38-pixelraster/16508393/0001.g3 /spt/data/bolodata/downsampled/RCW38-pixelraster/16508393/0002.g3 -o /spt/user/benderan/noise_test/rcw38_maps/rcw38_16508393.g3 -s rcw38

# Usage: makebolomaps.py <input files.g3> -o outputmaps.g3 -s rcw38
P = ap.ArgumentParser(description='Single bolometer maps with boresight pointing',
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input_files', action='store', nargs='+', default=[],
               help='Input files')
P.add_argument('-o', '--output', action='store', default='output.g3',
               help='Output filename')
P.add_argument('--input-ts', action='store',
               choices=['RawTimestreams_I', 'RawTimestreams_Q'],
               default='RawTimestreams_I', help='Input timestream type')
P.add_argument('-s', '--source', action='store', 
               default='rcw38', help='name of source')
#P.add_argument('-c', '--calpath', action='store', default=None,
#               help='path to directory of calibrator results')
#P.add_argument('-x', '--xtalkpath', action='store', default=None,
#               help='path to directory of crosstalk results')
P.add_argument('-r', '--res', action='store', 
               default=0.5, help='resolution [arcmin]')
args = P.parse_args()

res = float(args.res) * core.G3Units.arcmin


#Get the list of bolos to use from the calibrator scan
cal_scan='/spt/user/production/calibration/calibrator/16507005.g3'
cal_data=list(core.G3File(cal_scan))
all_bolos=cal_data[0]['CalibratorResponseSN'].keys()
cal_sn=np.zeros(len(all_bolos))
for ii,kk in enumerate(all_bolos):
    cal_sn[ii]=cal_data[0]['CalibratorResponseSN'][kk]
bolos=np.ndarray.tolist(np.array(all_bolos)[np.where(cal_sn>5)])
#bolos=bolos[0:2]
print 'using '+np.str(len(bolos))+' bolos'



for fname in args.input_files:
    for frame in core.G3File(fname):
        if args.input_ts in frame:
            starttime=frame[args.input_ts].start
            break


'''
# Guess the list of bolos to use
bolos = None
if not bolos:
    for fname in args.input_files:
        for frame in core.G3File(fname):
            if args.input_ts in frame:
                bolos = frame[args.input_ts].keys()
                starttime = frame[args.input_ts].start
                break
        if bolos is not None:
            break
'''
# Generate map stub
smstub = std_processing.CreateSourceMapStub(
    args.source, x_len = 3.*core.G3Units.deg/res,
    y_len = 3.*core.G3Units.deg/res, res = res,
    proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
    at_time = starttime)

pipe = core.G3Pipeline()

input = []

input += args.input_files
pipe.Add(core.G3Reader, filename=input)
pipe.Add(core.DeduplicateMetadata)

# Combine our various in-progress cal data
pipe.Add(calibration.build_cal_frames.MergeCalibrationFrames)

# convert to watts
pipe.Add(dfmux.ConvertTimestreamUnits, Input=args.input_ts, Output='BoloMapTimestreams', Units=core.G3TimestreamUnits.Power)

pipe.Add(std_processing.MakeBoresightBolometerProperties)

# Cut turnarounds
pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])

pipe.Add(core.Dump)

# Kick off maps

pipe.Add(MakeMap,
         # Timestream filtering -- these really shouldn't be part of MakeMap
         poly_order = 4,
         use_dynamic_source_filter = True, # Enable dynamic PS filtering

         # Actual map-making options
         map_in = smstub,
         map_id = 'bsmap',
         ts_in_key = 'BoloMapTimestreams',

         make_polarized = False,
         do_weight = True,
         fill_in_unity_weights = True,
         use_boresight_pointing = True,
         individual_bolos_to_map = bolos,
         boresight_ra_key = 'OnlineBoresightRa',
         boresight_dec_key = 'OnlineBoresightDec')

pipe.Add(lambda fr: fr.type != core.G3FrameType.Scan) # Drop TOD


pipe.Add(core.G3Writer, filename=args.output)
pipe.Run()


'''
savedir='/spt/user/benderan/noise_test/rcw38_maps/gif_plots/'
f=list(core.G3File('/spt/user/benderan/noise_test/rcw38_maps/rcw38_16268551.g3'))
for ff in range(5,len(f)):                                                          
    plt.imshow(f[ff]['T'])
    plt.title(f[ff]['Id'])
    plt.show()
    plt.savefig(os.path.join(savedir,'rcw38_20170710_'+str(ff)+'.png'))
    plt.clf()


'''

