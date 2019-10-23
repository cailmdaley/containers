import numpy as np
import os
from spt3g import core, gcp, std_processing

#start_time='22-Jan-2017:20:06:05'
#stop_time='22-Jan-2017:20:06:10'
start_time='28-Jan-2017:07:37:38'
stop_time='28-Jan-2017:07:37:39'
#start_time='24-Jan-2017:04:46:12'
#stop_time='24-Jan-2017:04:46:13'
arcdir='/buffer/arc/'

pipe = core.G3Pipeline()
pipe.Add(std_processing.ARCTimerangeReader, start_time=start_time, stop_time=stop_time, basedir=arcdir)
fntemp = '/home/tcrawfor/temp_eraseme_dump.g3'
pipe.Add(core.G3Writer, filename=fntemp)
pipe.Run()

frames = [frame for frame in core.G3File(fntemp)]
