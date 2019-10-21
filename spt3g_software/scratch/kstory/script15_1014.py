from spt3g import core, dfmux, calibration
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal


# f = core.G3File('/data52/nwhitehorn/sptpol.g3')
# c1 = f.next()
# f1 = f.next()

# # Get timestream for one bolometer
# pipe = core.G3Pipeline()
# pipe.Add(core.G3Reader, filename='/data52/nwhitehorn/sptpol.g3')
# pipe.Add(append_ts)
# pipe.Run(profile=True)

# ts = np.concatenate(bolo_data)
# plt.plot(ts)


# collect and Filter timestreams

def filter_ts(frame):
    if frame.type != core.G3FrameType.Scan:
        return
    print frame['CalTimestreams'].start
    tsmap = core.G3TimestreamMap()
    for bolo in frame['CalTimestreams'].keys():
        tsmap[bolo] = core.G3Timestream( signal.detrend(frame['CalTimestreams'][bolo]), units=frame['CalTimestreams'][bolo].units)

    frame['Filtered'] = tsmap

def collect_ts_one_bolo(frame, boloID, outputarray):
    if frame.type == core.G3FrameType.Scan:
        #bolo_data.append(frame['CalTimestreams']['D1.A2.15.Y'])
        outputarray.append(frame['Filtered'][boloID])

def collect_ts(frame, data):
    if frame.type != core.G3FrameType.Scan:
        return
    data.append(frame)


# Get timestream for one bolometer
bolo_data = [] # this will be a list for frames

pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename='/data52/nwhitehorn/sptpol.g3')
pipe.Add(filter_ts)
#pipe.Add(collect_ts, boloID='D1.A2.15.Y', outputarray=bolo_data)
pipe.Add(collect_ts, data=bolo_data)
pipe.Run(profile=True)

calts = bolo_data[1]['Filtered']
for key in calts.keys():
    value = calts[key]
    plt.plot(value)
