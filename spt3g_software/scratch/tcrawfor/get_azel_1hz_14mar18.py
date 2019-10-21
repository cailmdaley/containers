import numpy
from spt3g import core, std_processing, gcp
#from spt3g.scratch.tcrawfor import tctools
from spt3g.util import tctools

def grab_azel(f, azlist = [], ellist = []):
    try:
        azlist.append(f['antenna0']['tracker']['actual'][0][0])
        ellist.append(f['antenna0']['tracker']['actual'][1][0])
    except:
        pass

#start_time = '11-Jan-2018:16:00:00'
#stop_time = '11-Jan-2018:19:00:00'
#stop_time = '11-Jan-2018:16:10:00'
#start_time = '10-Jan-2018:16:00:00'
start_time = '11-Jan-2018:09:00:00'
stop_time = '12-Jan-2018:19:00:00'

#arcdir='/spt_data/arc/'
arcdir='/spt/data/arc/'

azdata = []
eldata = []

pipe1 = core.G3Pipeline()
pipe1.Add(std_processing.ARCTimerangeReader, start_time=start_time, stop_time=stop_time, basedir=arcdir)
pipe1.Add(gcp.ARCExtract)
pipe1.Add(grab_azel, azlist=azdata, ellist=eldata)
pipe1.Run()

naztemp = np.asarray(azdata)/core.G3Units.deg
neltemp = np.asarray(eldata)/core.G3Units.deg

