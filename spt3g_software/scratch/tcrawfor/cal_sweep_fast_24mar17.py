import numpy as np
from spt3g import core, std_processing, dfmux
from spt3g.util import tctools
import glob

dir1 = '/spt/user/production/calibration/calibrator/'
files = glob.glob(dir1 + '71*.g3')
files.sort()
files = files[0:6]

freqs = [6.,30.]
frame6 = (core.G3File(files[0])).next()
frame30 = (core.G3File(files[3])).next()

