import numpy as np
from spt3g import core, std_processing, gcp
#from spt3g.scratch.tcrawfor import tctools

def grab_boards(f, boardlist = []):
    keylist = ['TrackerPointing',
               'ACUStatus',
               'SourceName',
               'array',
               'TrackerStatus',
               'GCPFeatureBits',
               'antenna0']
    try:
        for key in keylist:
            boardlist.append(f[key])
    except:
        pass

#start_time = '31-Jan-2017:17:57:09.930000000'
#stop_time = '31-Jan-2017:17:58:09.930000000'
start_time = '04-Feb-2017:18:27:17'
stop_time = '04-Feb-2017:18:28:07'

arcdir='/spt_data/arc/'

data1 = []

pipe1 = core.G3Pipeline()
pipe1.Add(std_processing.ARCTimerangeReader, start_time=start_time, stop_time=stop_time, basedir=arcdir)
pipe1.Add(gcp.ARCExtract)
pipe1.Add(grab_boards, boardlist=data1)
#    pipe1.Add(core.G3Writer, filename=fntemp1)
pipe1.Run()

nseconds = len(data1)/7

tstemp = []
anttemp = []
for i in np.arange(nseconds):
    tstemp.append(data1[7*i+4])
    anttemp.append(data1[7*i+6])

azpos = []
azpos2 = []
elpos = []
elpos2 = []
for i in np.arange(nseconds): 
    azpos.append(tstemp[i].az_pos)
    azpos2.append(anttemp[i]['tracker']['actual'][0])
    elpos.append(tstemp[i].el_pos)
    elpos2.append(anttemp[i]['tracker']['actual'][1])

nazpos = tctools.list_to_array(azpos)
nazpos2 = tctools.list_to_array(azpos2)
nelpos = tctools.list_to_array(elpos)
nelpos2 = tctools.list_to_array(elpos2)
