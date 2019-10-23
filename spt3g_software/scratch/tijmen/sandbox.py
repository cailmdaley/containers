import numpy as np
import scipy.signal
from spt3g import core
pipe = core.G3Pipeline()

ts = []
def add_my_data(data):
  if data.type == core.G3FrameType.Scan:
    ts.append(data)

def detrend_frame(frame):
  if frame.type != core.G3FrameType.Scan:
    return
  tsmap = core.G3TimestreamMap()
  for bolo in frame['CalTimestreams'].keys():
    tsmap[bolo] = core.G3Timestream(scipy.signal.detrend(frame['CalTimestreams'][bolo]),
                                    units=frame['CalTimestreams'][bolo].units)
  frame['FilteredTimestreams'] = tsmap

pipe.Add(core.G3Reader, filename='/data52/nwhitehorn/sptpol.g3')
pipe.Add(detrend_frame)
pipe.Add(add_my_data)
pipe.Run()

calts = ts[1]['FilteredTimestreams']
from matplotlib import pyplot as plt
plt.ion()
for key in calts.keys():
  value = calts[key]
  #this_timestream = value - np.median(value) 
  plt.plot(value)
