import numpy
from spt3g import core, std_processing, gcp

def grab_data(f, data1 = []):
    try:
        data1.append(f['array']['cryo']['temperature'])
    except:
        pass

start_time=core.G3Time('07-Feb-2017:11:00:00')
#start_time=core.G3Time('08-Feb-2017:11:00:00')
nhr = 2
arcdir='/spt_data/arc/'

data1 = []

for i in np.arange(nhr):
    this_start_time = start_time + 1.*core.G3Units.hour
    this_stop_time = this_start_time + 9.*core.G3Units.s
    pipe1 = core.G3Pipeline()
    pipe1.Add(std_processing.ARCTimerangeReader, start_time=this_start_time, stop_time=this_stop_time, basedir=arcdir)
    pipe1.Add(gcp.ARCExtract)
    pipe1.Add(grab_data, data1=data1)
#    pipe1.Add(core.G3Writer, filename=fntemp1)
    pipe1.Run()

