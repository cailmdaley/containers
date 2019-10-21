from spt3g import core, calibration
import numpy, sys

# Dumps the focal plane positions from a single-calibration-frame .g3 file
# to text. Can handle calibration files with a bolometer properties map
# or raw PointingOffsetX/Y data. In the former case, dumps the bands to a third
# column.

fr = core.G3File(sys.argv[1]).next()
print(fr)

if 'BolometerProperties' in fr:
	numpy.savetxt(sys.argv[2], [(x.x_offset, x.y_offset, x.band/core.G3Units.GHz) for x in fr['BolometerProperties'].values()])
elif 'NominalBolometerProperties' in fr:
	numpy.savetxt(sys.argv[2], [(x.x_offset, x.y_offset, x.band/core.G3Units.GHz) for x in fr['BolometerProperties'].values()])
else:
	numpy.savetxt(sys.argv[2], numpy.asarray(list(zip(fr['PointingOffsetX'].values(), fr['PointingOffsetY'].values()))))

