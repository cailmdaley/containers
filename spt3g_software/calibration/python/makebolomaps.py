from spt3g import core, mapmaker, dfmux, coordinateutils, std_processing
from spt3g.mapmaker.mapmakerutils import MakeMap
from spt3g.calibration.apply_t_cal import ApplyTCalibration
import sys

import argparse as ap

def makebolomaps(input_files, output='output.g3', source='rcw38', input_ts='RawTimestreams_I', bolos=[],
                 units = core.G3TimestreamUnits.Power, res=0.5 * core.G3Units.arcmin, map_width=3.3 * core.G3Units.deg):
    '''
    Purpose: Make a set of individual bolometer maps

    Inputs:
      input_files [list]:          list of input intermediate data files in g3 format
      output [str]:                name and full path of output file
      source [str]:                source name for this map
      input_ts [str]:              input timestream key
      bolos [list]:                list of bolo names to use -- uses all bolos in file by default

    Outputs: (none)
      Writes a g3 file to [output]
    '''
    
    smstub = std_processing.CreateSourceMapStub(
        source, x_len = map_width/res, y_len = map_width/res, res = res,
        proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea)

    # Grub through the input files to get the set of bolometers
    bolos = bolos
    if not bolos:
        for fname in input_files:
            for frame in core.G3File(fname):
                if 'WiringMap' in frame:
                    bolos = frame['WiringMap'].keys()
                break
            if input_ts in frame:
                bolos = frame[input_ts].keys()
                break
            if bolos is not None and len(bolos) > 0:
                break

    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=input_files)

    # Cut turnarounds
    pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])

    # Convert to target units
    if units == core.G3TimestreamUnits.Tcmb:
        pipe.Add(dfmux.ConvertTimestreamUnits, Input=input_ts, Output='BoloMapTimestreamsWatts', Units=core.G3TimestreamUnits.Power)
        pipe.Add(ApplyTCalibration, Input='BoloMapTimestreamsWatts', Output='BoloMapTimestreams') # , SkyCal='RCW38FluxCalibration') <- For RCW38 rather Kcmb
    else:
        pipe.Add(dfmux.ConvertTimestreamUnits, Input=input_ts, Output='BoloMapTimestreams', Units=units)

    # Kick off maps
    pipe.Add(MakeMap,
         # Timestream filtering -- these really shouldn't be part of MakeMap
         poly_order = 4,
         use_dynamic_source_filter = True, # Enable dynamic PS filtering
         dynamic_masking_debug = False,

         # Actual map-making options
         map_in = smstub,
         map_id = 'bsmap',
         ts_in_key = 'BoloMapTimestreams',

         make_polarized = False,
         do_weight = True,
         fill_in_unity_weights = True,
         individual_bolos_to_map = bolos,
         boresight_ra_key = 'BoresightRa',
         boresight_dec_key = 'BoresightDec')

    pipe.Add(core.Dump)
    pipe.Add(lambda fr: fr.type != core.G3FrameType.Scan)
    pipe.Add(core.G3Writer, filename=output)
    pipe.Run(profile = True)

