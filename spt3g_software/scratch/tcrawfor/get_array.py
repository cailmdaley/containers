import numpy
from spt3g import core, std_processing, gcp

def grab_boards(f, boardlist = []):
#    keylist = ['TrackerPointing',
#               'ACUStatus',
#               'SourceName',
#               'array',
#               'TrackerStatus',
#               'GCPFeatureBits',
    key = 'array'
#    key = 'TrackerPointing'
    try:
        boardlist.append(f[key])
    except:
        pass

#start_time='29-Jan-2017:11:00:15'
#stop_time='29-Jan-2017:11:01:15'
start_time='08-Feb-2017:11:00:00'
stop_time='09-Feb-2017:09:00:00'
#start_time='29-Jan-2017:11:05:41'
#stop_time='29-Jan-2017:11:06:41'
#    start_time='28-Jan-2017:00:20:00'
#    stop_time='28-Jan-2017:00:20:00'
arcdir='/spt_data/arc/'

data1 = []

pipe1 = core.G3Pipeline()
pipe1.Add(std_processing.ARCTimerangeReader, start_time=start_time, stop_time=stop_time, basedir=arcdir)
pipe1.Add(gcp.ARCExtract)
pipe1.Add(grab_boards, boardlist=data1)
#    pipe1.Add(core.G3Writer, filename=fntemp1)
pipe1.Run()

