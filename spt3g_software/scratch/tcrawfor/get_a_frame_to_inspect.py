import numpy as np
import os
from spt3g import core, gcp, std_processing

#start_time='22-Jan-2017:20:06:05'
#stop_time='22-Jan-2017:20:06:10'
#start_time='28-Jan-2017:17:45:50'
#stop_time='28-Jan-2017:17:45:50'
start_time='28-Jan-2017:19:38:01'
stop_time='28-Jan-2017:19:38:01'
#start_time='28-Jan-2017:23:06:11'
#stop_time='28-Jan-2017:23:06:11'
#start_time='24-Jan-2017:04:46:12'
#stop_time='24-Jan-2017:04:46:13'
#start_time='29-Jan-2017:00:00:01'
#stop_time='29-Jan-2017:00:00:01'
arcdir='/buffer/arc/'

pipe = core.G3Pipeline()
pipe.Add(std_processing.ARCTimerangeReader, start_time=start_time, stop_time=stop_time, basedir=arcdir)
fntemp = '/home/tcrawfor/temp_eraseme_dump.g3'
pipe.Add(core.G3Writer, filename=fntemp)
pipe.Run()

for frame in core.G3File(fntemp):
    print frame['array']['frame']['utc']
    print np.min(frame['antenna0']['tracker']['actual'][0])
    print np.max(frame['antenna0']['tracker']['actual'][0])
    print frame
