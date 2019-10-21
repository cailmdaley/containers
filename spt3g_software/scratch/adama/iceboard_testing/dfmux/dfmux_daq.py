# dfmux_daq.py
#
# A control script for data acquisition from the dfmux board.
# This script takes no command line arguments, but see the
# "OPTIONS" section for available options. They should be self-
# explanatory.
#
# Adam Anderson
# adama@fnal.gov

from spt3g import core, dfmux
import pydfmux
import os
import datetime as dt

# OPTIONS
hwm_file = '/home/sptdaq/iceboard_testing/hwm_03_hf_bias_skyarcs/dongle.yaml'
out_path = '/data/sptdaq/iceboard_testing/dfmux/raw/'
set_0_bias = False   # true disables DAN and runs with 0 carrier amplitude, otherwise overbias_and_null
use_irig = False     # true uses the IRIG for timekeeping

print('Initializing hardware map and boards')
hwm = pydfmux.load_session(open(hwm_file, 'rb'))['hardware_map']

for board in hwm.query(pydfmux.Dfmux):
	print 'Initializing board'
	board.clear_all()
	board.set_mezzanine_power(True, 1)
	board.set_mezzanine_power(True, 2)

	print 'Setting time source'
        if use_irig == True:
                board.set_timestamp_port(board.TIMESTAMP_PORT.SMA_B)
        else:
                board.set_timestamp_port(board.TIMESTAMP_PORT.TEST)
	
print('Aligning board sampling')
hwm.query(pydfmux.core.dfmux.IceBoard).set_fir_stage(6)
hwm.query(pydfmux.core.dfmux.IceBoard).align_sampling()

print('Configuring bolometers')
if set_0_bias == True:
        for bolo in hwm.query(pydfmux.core.dfmux.Bolometer):
                frequency = bolo.readout_channel.channel_map.lc_channel.frequency
                bolo.readout_channel.set_frequency(frequency, 'Hz', 'CARRIER')
                bolo.readout_channel.set_frequency(frequency, 'Hz', 'NULLER')
                bolo.readout_channel.set_frequency(frequency, 'Hz', 'DEMOD')
else:
     hwm.query(pydfmux.core.dfmux.Bolometer).overbias_and_null(serialize=True)   

# set up the pipeline
datenow = dt.datetime.utcnow()
pipe = core.G3Pipeline()
builder = dfmux.DfMuxBuilder(len(hwm.query(pydfmux.core.dfmux.IceBoard).all()))
collector = dfmux.DfMuxCollector("192.168.5.200", builder)
pipe.Add(builder)

# Insert current hardware map into data stream. This is critical to get the
# board ID -> IP mapping needed to do anything useful with the data
pipe.Add(dfmux.PyDfMuxHardwareMapInjector, pydfmux_hwm=hwm)
#pipe.Add(core.Dump)
pipe.Add(core.G3MultiFileWriter, filename=out_path+'/dfmux_%4d%02d%02d_%02d%02d%02d_%%02u.g3'%(datenow.year, datenow.month, datenow.day, datenow.hour, datenow.minute, datenow.second), size_limit=(1024**3))

print('starting DfMux collector')
collector.Start()
print('acquiring data')
pipe.Run()
