import numpy
from spt3g import core, std_processing, gcp
from spt3g.util import tctools

def grab_features(f, data1 = []):
    try:
        data1.append(f['array']['frame']['features'])
    except:
        pass

arcdir='/spt_data/arc/'
arcdir='/spt/data/arc/'

# change this to your favorite chunk of time
start_time=core.G3Time('20171027_034800')
stop_time=core.G3Time('20171027_035900')

data1 = []

pipe1 = core.G3Pipeline()
pipe1.Add(std_processing.ARCTimerangeReader, start_time=start_time, stop_time=stop_time, basedir=arcdir)
pipe1.Add(grab_features, data1 = data1)
pipe1.Run()

features = np.zeros(len(data1),dtype=int)
for i in np.arange(len(data1)):
    features[i] = data1[i].value


