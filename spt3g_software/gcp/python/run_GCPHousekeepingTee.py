#!/usr/bin/env python

from spt3g import core, dfmux, gcp, auxdaq
import socket, argparse, os

parser = argparse.ArgumentParser(description='Send Housekeeping data to the GCP')
parser.add_argument('--hardware_map', dest='hardware_map', default=os.path.join(os.environ.get('HWM_DIR'),'2018','hwm_pole_run2','hwm_squids.yaml'), 
		    metavar='/path/to/hwm.yaml', help='Path to hardware map YAML file')
parser.add_argument('-v', dest='verbose', action='store_true', help='Verbose mode (print all frames)')
args = parser.parse_args()

# Import pydfmux later since it can take a while
import pydfmux

#Sanity-check hardware_map
if not os.path.exists(args.hardware_map):
    raise RuntimeError("Given path to hardware_map does not exist: %s" % args.hardware_map)

print('Initializing hardware map and boards')
hwm = pydfmux.load_session(open(args.hardware_map, 'r'))['hardware_map']
if hwm.query(pydfmux.core.dfmux.IceCrate).count() > 0:
    hwm.query(pydfmux.core.dfmux.IceCrate).resolve()

print('Beginning data acquisition')
# Set up DfMux consumer
pipe = core.G3Pipeline()
builder = dfmux.DfMuxBuilder([int(board.serial) for board in hwm.query(pydfmux.IceBoard)])

# Get the local IP(s) to use to connect to the boards by opening test
# connections. Using a set rather than a list deduplicates the results.
local_ips = set()
for board in hwm.query(pydfmux.core.dfmux.IceBoard):
    print('Found iceboard{}. Attempting to connect:'.format(board.serial))
    try:
        testsock = socket.create_connection(('iceboard{}.local'.format(board.serial), 80))
        sockname = testsock.getsockname()[0]
        local_ips.add(sockname)
        testsock.close()
        print('Successfully connected to iceboard{} at {}'.format(board.serial, sockname))
    except Exception as e:
        print('Failed to open socket to iceboard{}'.format(board.serial))
        print(e)
        raise(e)
print('Creating listeners for %d boards on interfaces: %s' % (hwm.query(pydfmux.core.dfmux.IceBoard).count(), ', '.join(local_ips)))

# Set up listeners per network segment and point them at the event builder
collectors = [dfmux.DfMuxCollector(ip, builder) for ip in local_ips]
pipe.Add(builder)

# Insert current hardware map into data stream. This is critical to get the
# board ID -> IP mapping needed to do anything useful with the data
print('Add PyDfMuxHardwareMapInjector')
pipe.Add(dfmux.PyDfMuxHardwareMapInjector, pydfmux_hwm=hwm)

# Collect housekeeping when GCP asks for it. Housekeeping consumer runs in a
# different process to reduce pipeline stalls waiting for boards to respond.
print('Add GCPHousekeepingTee')
pipe.Add(gcp.GCPSignalledHousekeeping)
# Rather than setting subprocess=True, instantiate subproc by hand with a
# long queue length (8 seconds @ 152 Hz) to avoid standard 50 frame
# (= 0.3 second) limit. Typical response time from the boards to housekeeping
# collection requests is 2 seconds.
pipe.Add(core.Subproc(dfmux.HousekeepingConsumer(), name='Housekeeping Consumer', maxqueuelen=10000))
pipe.Add(gcp.GCPHousekeepingTee(verbose=True))

if args.verbose:
    pipe.Add(core.Dump)

for collector in collectors:
    collector.Start()
pipe.Run()

