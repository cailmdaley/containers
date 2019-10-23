# rawdump_daq.py
#
# Simple python script to control acquisition of rawdump data. We take data
# in 10 minute file blocks, repeating indefinitely the 'iceboard_da.py' script
# issues a kill command.
# 
# Adam Anderson
# adama@fnal.gov

import os, sys, signal
import datetime as dt

if len(sys.argv) != 3:
    print 'usage: python rawdump_daq.py [output filename]'

record_time = 10*60    # [sec]

rawdump_binary_path=sys.argv[2]
rawdump_binary='rdrecord'

while True:
    dtnow = dt.datetime.utcnow()
    out_filename = 'rawdump_%d'%dtnow.year + '%d'%dtnow.month + '%d_'%dtnow.day + \
                   '%d'%dtnow.hour + '%d'%dtnow.minute + '%d'%dtnow.second
    
    daq_pid = os.spawnv(os.NOWAIT, rawdump_binary_path + rawdump_binary,
                        [rawdump_binary_path + rawdump_binary, out_filename])
    
    # wait for specified length, then kill daq and start a new file
    time.sleep(record_time)
    os.kill(daq_pid, signal.SIGKILL)
    
