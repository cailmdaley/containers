import glob
import numpy
from spt3g import core, std_processing, gcp

def grab_data(f, data1 = []):
    try:
        data1.append(f['array']['weather']['windSpeed'])
    except:
        pass

#start_time=core.G3Time('01-Aug-2018:00:00:00')
start_time=core.G3Time('10-Aug-2018:00:00:00')
#nhr = 120
nhr = 20

arcdir='/spt/data/arc/'

windspeed = np.zeros(nhr)

for i in np.arange(nhr):
    data1 = []
#    this_start_time = start_time + 2.*float(i)*core.G3Units.hour
    this_start_time = start_time + 12.*float(i)*core.G3Units.hour
    this_stop_time = this_start_time + 10.*core.G3Units.s
    pipe1 = core.G3Pipeline()
    pipe1.Add(std_processing.ARCTimerangeReader, start_time=this_start_time, stop_time=this_stop_time, basedir=arcdir)
    pipe1.Add(gcp.ARCExtract)
    pipe1.Add(grab_data, data1=data1)
#    pipe1.Add(core.G3Writer, filename=fntemp1)
    pipe1.Run()
    print this_start_time
    print data1[0]
    data2 = np.asarray([d1.value for d1 in data1])
    windspeed[i] = np.mean(data2)



