import pydfmux
from pydfmux.core import dfmux
from pydfmux.algorithms.bolo import overbias_and_null
from pydfmux.algorithms.initialization.initialize_hardware import initialize_iceboard, initialize_squidcontroller

hwm = pydfmux.load_session(open('/home/sptdaq/pydfmux/spt3g/hardware_maps/dongle.yaml', 'r'))['hardware_map']

for board in hwm.query(pydfmux.Dfmux):
	print 'Initializing board'
	board.clear_all()
	board.set_mezzanine_power(True, 1)
	board.set_mezzanine_power(True, 2)
	print 'Setting time source'
	board.set_timestamp_port(board.TIMESTAMP_PORT.TEST)
        #board.set_timestamp_port(board.TIMESTAMP_PORT.SMA_B)

#hwm.query(pydfmux.Bolometer).call_with(overbias_and_null.overbias_and_null) # crashes for unknown reason

