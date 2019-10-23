import numpy
from spt3g import core, std_processing, gcp
#from spt3g.scratch.tcrawfor import tctools
#from spt3g.util import tctools

def list_to_array(list_in):

    npts_tot = 0
    for item in list_in:
        npts_tot += len(item)
    array_out = np.zeros(npts_tot)
    npts = 0
    for item in list_in:
        array_out[npts:npts+len(item)] = item
        npts += len(item)

    return array_out

def grab_reg(f, reg='actual', data1=[]):
    key = 'antenna0'
    key2 = 'tracker'
    try:
        data1.append(f[key][key2][reg])
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
#start_time = '31-Jan-2017:17:02:00'
#stop_time = '31-Jan-2017:18:45:00'
#start_time = '05-Aug-2017:14:06:29'
#stop_time = '05-Aug-2017:16:09:19'
#start_time = '05-Aug-2017:16:09:19'
#stop_time = '05-Aug-2017:16:19:19'
#start_time = '23-Apr-2018:05:00:00'
#stop_time = '23-Apr-2018:05:30:00'
#start_time = '01-Jul-2018:12:50:00'
#stop_time = '01-Jul-2018:12:59:00'
#start_time = '09-Jul-2018:09:44:00'
#stop_time = '09-Jul-2018:09:53:00'
#start_time = '12-Jul-2018:16:52:30'
#stop_time = '12-Jul-2018:16:55:00'
#start_time = '25-Feb-2019:02:34:24'
#stop_time = '25-Feb-2019:05:14:46'
#start_time = '25-Feb-2019:18:02:09'
#stop_time = '25-Feb-2019:20:41:54'
#start_time = '15-Jun-2019:10:35:00'
#stop_time = '15-Jun-2019:10:53:00'
start_time = '14-Jul-2019:00:37:46'
stop_time = '14-Jul-2019:00:39:30'

arcdir='/spt/data/arc/'

data1 = []
data2 = []

pipe1 = core.G3Pipeline()
pipe1.Add(std_processing.ARCTimerangeReader, start_time=start_time, stop_time=stop_time, basedir=arcdir)
pipe1.Add(gcp.ARCExtract)
pipe1.Add(grab_reg, reg='actual', data1=data1)
pipe1.Add(grab_reg, reg='lst', data1=data2)
#    pipe1.Add(core.G3Writer, filename=fntemp1)
pipe1.Run()

aztemp = []
eltemp = []
for data in data1:
    aztemp.append(data[0])
    eltemp.append(data[1])
lsttemp = []
for data in data2:
    lsttemp.append(data[0])

naztemp = list_to_array(aztemp)/core.G3Units.deg
neltemp = list_to_array(eltemp)/core.G3Units.deg
#ndectemp = -neltemp
#nlsttemp = list_to_array(lsttemp)/core.G3Units.h
#nratemp = 
