#!/usr/bin/env python

from spt3g import core, dfmux, gcp
import socket, argparse, os

parser = argparse.ArgumentParser(description='Record dfmux data to disk')
parser.add_argument('hardware_map', metavar='/path/to/hwm', help='Path to hardware map directory')
parser.add_argument('output', metavar='/path/to/files/file-%03u.g3', help='Path for output files')
parser.add_argument('port', metavar=8734, help='Port to listen on')

parser.add_argument('-v', dest='verbose', action='store_true', help='Verbose mode (print all frames)')
parser.add_argument('--max_file_size', dest='max_file_size', default=1024, help='Maximum file size in MB (default 1024)')
args = parser.parse_args()

import pywtl.common.DFML.HWM.HardwareMap as HWM

print('Initializing hardware map')
hwm = HWM.HardwareMap(args.hardware_map)

print('Beginning data acquisition')
# Set up DfMux consumer
pipe = core.G3Pipeline()
builder = dfmux.DfMuxBuilder([int(hwm(b, 'dfmux_id')) for b in hwm.boards])

# Set up listener and point at the event builder
collector = dfmux.LegacyDfMuxCollector(port=args.port, builder=builder)
pipe.Add(builder)

# Insert current hardware map into data stream. This is critical to get the
# board ID -> IP mapping needed to do anything useful with the data
pipe.Add(dfmux.DfmlHardwareMapInjector, dfml_hwm=hwm)

# Collect housekeeping when GCP asks for it. Housekeeping consumer runs in a
# different process to reduce pipeline stalls waiting for boards to respond.
pipe.Add(gcp.GCPSignalledHousekeeping)
# Rather than setting subprocess=True, instantiate subproc by hand with a
# long queue length (8 seconds @ 152 Hz) to avoid standard 50 frame
# (= 0.3 second) limit. Typical response time from the boards to housekeeping
# collection requests is 2 seconds.
pipe.Add(core.Subproc(dfmux.LegacyHousekeepingConsumer(), name='Housekeeping Consumer', maxqueuelen=1000))
pipe.Add(gcp.GCPHousekeepingTee)

if args.verbose:
	pipe.Add(core.Dump)

pipe.Add(core.G3MultiFileWriter, filename=args.output, size_limit=args.max_file_size*1024*1024)

collector.Start()
pipe.Run()

