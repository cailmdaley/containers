from spt3g import core, dfmux
import time
import numpy as np

fr = core.G3Frame(core.G3FrameType.Timepoint)
tm = core.G3TimestreamMap()
for i in range(100):
    tm['%d'%i] = core.G3Timestream(np.random.normal(size = 1000))
fr['TimestreamMap'] = tm

def VomitModule(frame):
    return [fr]
pipe = core.G3Pipeline()
pipe.Add(VomitModule)
pipe.Add(core.G3NetworkSender, port=8445)
pipe.Run()

