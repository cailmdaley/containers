import numpy
from spt3g import core, std_processing, gcp
#from spt3g.scratch.tcrawfor import tctools
from spt3g.util import tctools

def grab_boards(f, boardlist = []):
#    keylist = ['TrackerPointing',
#               'ACUStatus',
#               'SourceName',
#               'array',
#               'TrackerStatus',
#               'GCPFeatureBits',
    key = 'antenna0'
#    key = 'TrackerPointing'
    try:
        boardlist.append(f[key])
    except:
        pass

#start_time='29-Jan-2017:11:00:15'
#stop_time='29-Jan-2017:11:01:15'
#start_time='29-Jan-2017:11:05:41'
#stop_time='29-Jan-2017:11:06:40'
#start_time='28-Jan-2017:00:20:00'
#stop_time='28-Jan-2017:00:20:00'
#start_time='29-Jan-2017:10:26:00'
#stop_time='29-Jan-2017:10:26:59'
#start_time = '31-Jan-2017:19:07:00'
#stop_time = '31-Jan-2017:19:13:00'
#start_time = '31-Jan-2017:17:30:00'
#stop_time = '31-Jan-2017:17:31:00'
#start_time = '30-Jan-2017:10:00:00'
#stop_time = '30-Jan-2017:10:01:00'
#start_time = '02-Feb-2017:16:00:00'
#stop_time = '02-Feb-2017:16:06:00'
#start_time = '02-Feb-2017:15:51:43'
#stop_time = '02-Feb-2017:15:52:34'
#start_time = '02-Feb-2017:23:30:00'
#stop_time = '02-Feb-2017:23:30:10'
#start_time = '04-Feb-2017:16:18:50'
#stop_time = '04-Feb-2017:18:11:29'
#start_time = '02-Feb-2017:17:38:39'
#stop_time = '02-Feb-2017:17:39:28'
#start_time = '05-Feb-2017:18:42:17.000000320'
#stop_time = '05-Feb-2017:18:43:06.000000320'
#start_time = '20-Feb-2018:04:30:00'
#stop_time = '20-Feb-2018:04:32:00'
#start_time = '03-Feb-2018:10:00:00'
#stop_time = '03-Feb-2018:10:03:00'
#start_time = '08-Mar-2018:09:30:00'
#stop_time = '08-Mar-2018:10:30:00'
start_time = '11-Jan-2018:16:00:00'
stop_time = '11-Jan-2018:19:00:00'
#start_time = '11-Jan-2018:15:00:00'
#stop_time = '11-Jan-2018:16:00:00'


#arcdir='/spt_data/arc/'
arcdir='/spt/data/arc/'

data1 = []

pipe1 = core.G3Pipeline()
pipe1.Add(std_processing.ARCTimerangeReader, start_time=start_time, stop_time=stop_time, basedir=arcdir)
pipe1.Add(gcp.ARCExtract)
pipe1.Add(grab_boards, boardlist=data1)
#    pipe1.Add(core.G3Writer, filename=fntemp1)
pipe1.Run()

aztemp = []
eltemp = []
azrtemp = []
elrtemp = []
eaztemp = []
eeltemp = []
for frame in data1:
    aztemp.append(frame['tracker']['actual'][0])
    eltemp.append(frame['tracker']['actual'][1])
    azrtemp.append(frame['tracker']['actual_rates'][0])
    elrtemp.append(frame['tracker']['actual_rates'][1])
    eaztemp.append(frame['tracker']['errors'][0])
    eeltemp.append(frame['tracker']['errors'][1])

naztemp = tctools.list_to_array(aztemp)/core.G3Units.deg
neltemp = tctools.list_to_array(eltemp)/core.G3Units.deg
nazrtemp = tctools.list_to_array(azrtemp)/(core.G3Units.deg/core.G3Units.s)
nelrtemp = tctools.list_to_array(elrtemp)/(core.G3Units.deg/core.G3Units.s)
neaztemp = tctools.list_to_array(eaztemp)/core.G3Units.deg
neeltemp = tctools.list_to_array(eeltemp)/core.G3Units.deg

