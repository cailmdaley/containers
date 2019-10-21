#!/usr/bin/env python

from spt3g import core, dfmux, gcp, auxdaq
import socket, argparse, os

pipe = core.G3Pipeline()
# boards = [123, 54, 151, 96, 64]
boards = [123]
builder = dfmux.DfMuxBuilder(boards)

# Get the local IP(s) to use to connect to the boards by opening test
# connections. Using a set rather than a list deduplicates the results.
local_ips = set()
for board in boards:
    testsock = socket.create_connection(('iceboard' + '%04d'%board + '.local', 80))
    local_ips.add(testsock.getsockname()[0])
    testsock.close()
# print('Creating listeners for boards on interfaces: %s' % (', '.join(local_ips)))
# local_ips = set([
print(local_ips)

# Set up listeners per network segment and point them at the event builder
collectors = [dfmux.DfMuxCollector(ip, builder) for ip in local_ips]
pipe.Add(builder)

# pipe.Add(core.Dump)

for collector in collectors:
    collector.Start()
pipe.Run(profile=True)

# Shut everything down
for collector in collectors:
    collector.Stop()

# C++ global destructor runs after python has closed, sometimes
# Setting the logger to None like this explicitly avoids segfaults
core.G3Logger.global_logger = None
