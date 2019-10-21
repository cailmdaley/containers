#!/usr/bin/env python

from spt3g import core, dfmux, gcp
import socket, argparse

parser = argparse.ArgumentParser(description='Record dfmux data to a NetCDF file', prog='ledgerman')
parser.add_argument('hardware_map', metavar='/path/to/hwm.yaml', help='Path to hardware map YAML file')

parser.add_argument('-v', dest='verbose', action='store_true', help='Verbose mode (print all frames)')
parser.add_argument('-a', dest='align', action='store_true', help='Align sampling. This has to happen once, but will break any existing DAN loops when run.')
parser.add_argument('-s', dest='system_time', action='store_true', help='Replace board time with system time when data received. Useful if your board is in IRIG_TEST mode.')
parser.add_argument('--no_bprops', action='store_true')
args = parser.parse_args()

# Import pydfmux later since it can take a while
import pydfmux

print('Initializing hardware map and boards')
hwm = pydfmux.load_session(open(args.hardware_map, 'r'))['hardware_map']
if hwm.query(pydfmux.IceCrate).count() > 0:
        hwm.query(pydfmux.IceCrate).resolve()

if args.align:
	print('Aligning board sampling, this will break any existing DAN loops!')
	hwm.query(pydfmux.IceBoard).set_fir_stage(6)
	hwm.query(pydfmux.IceBoard).align_sampling()

print('Beginning data acquisition')

# Set up DfMux consumer
pipe = core.G3Pipeline()
builder = dfmux.DfMuxBuilder([int(board.serial) for board in hwm.query(pydfmux.IceBoard)])

# Get the local IP to connect to the boards by opening a test connection.
testsock = socket.create_connection(('iceboard' + str(hwm.query(pydfmux.IceBoard).first().serial) + '.local', 80))
local_ip = testsock.getsockname()[0]
testsock.close()

# Build mapping dictionary for old (64x) firmware
v2_mapping = {'iceboard' + str(b.serial) + '.local': int(b.serial) for b in hwm.query(pydfmux.IceBoard)}

collector = dfmux.DfMuxCollector(local_ip, builder, v2_mapping)
pipe.Add(builder)

# Insert current hardware map into data stream. This is critical to get the
# board ID -> IP mapping needed to do anything useful with the data
pipe.Add(dfmux.PyDfMuxHardwareMapInjector, pydfmux_hwm=hwm)

if not args.no_bprops:
        pipe.Add(dfmux.PyDfMuxBolometerPropertiesInjector, pydfmux_hwm=hwm)

if args.system_time:
	def sub_system_time(frame):
		if frame.type != core.G3FrameType.Timepoint:
			return
		del frame['EventHeader']
		frame['EventHeader'] = core.G3Time.Now()
	pipe.Add(sub_system_time)

if args.verbose:
	pipe.Add(core.Dump)

pipe.Add(gcp.GCPSignalledHousekeeping)
pipe.Add(dfmux.HousekeepingConsumer)

pipe.Add(core.G3ThrottledNetworkSender,
         port = 8675,
         max_queue_size = 10,
         frame_decimation = {core.G3FrameType.Timepoint: 10})

collector.Start()
pipe.Run()

