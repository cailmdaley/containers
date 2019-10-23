from spt3g import core, auxdaq
import time

#core.set_log_level(core.G3LogLevel.LOG_DEBUG)

cd = auxdaq.CalibratorDAQ()
i= 0
while 1:
    time.sleep(0.1)
    i+=2000
    frame = core.G3Frame(core.G3FrameType.Timepoint)
    frame['EventHeader'] = core.G3Time.Now()
    out = cd(frame)
    #for i1,j1 in frame.iteritems(): print i1,j1
    for f in out:
        print f

