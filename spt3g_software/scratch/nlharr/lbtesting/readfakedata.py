from spt3g import core, dfmux
import time
import numpy as np

pipe = core.G3Pipeline()
pipe.Add(core.G3NetworkReceiver, port = 8446, hostname = 'localhost')
pipe.Add(core.Dump)
pipe.Run()
