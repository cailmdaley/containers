import pydfmux, sys

#Usage: rdinit.py squidtodump

modtodump = int(sys.argv[1])

d = pydfmux.Dfmux(serial='0112')
d.set_mezzanine_power(True, 1)
d.set_mezzanine_power(True, 2)
d._set_fpga_tx_ip('255.255.255.255')
d._set_fpga_tx_mac('ff:ff:ff:ff:ff:ff')
d._fpga_spi_poke(0x110000, 0x80000000L + modtodump - 1)
raw_input("Rawdump begun from module %s. Press enter to stop. Don't kill!" % modtodump)

d._set_fpga_tx_ip('192.168.254.1')
d._set_fpga_tx_mac('84:7e:40:6d:37:e9')
d._fpga_spi_poke(0x110000, 0)

