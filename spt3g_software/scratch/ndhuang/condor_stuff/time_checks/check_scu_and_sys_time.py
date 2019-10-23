import os
import sys
import time
import numpy as np
from spt3g import core, gcp
from frame_grabber import FrameGrabber

class scu_check(object):
    def __init__(self):
        self.outages = 0
        self.outage_time = 0.
        self.last = 1
    
    def __call__(self, fr):
        if fr.type != core.G3FrameType.GcpSlow:
            return
        timeLocked = np.zeros(101)
        timeLocked[0] = self.last
        timeLocked[1:] = fr['antenna0']['scu']['timeLocked'][0]
        inds = np.where(timeLocked == 0)[0]
        self.last = timeLocked[-1]
        if len(inds) == 0:
            return
        self.outage_time += len(inds) * .01
        self.outages += len(np.where(np.diff(inds) > 1)[0]) + 1

def drop_receiver(fr):
    fr.pop('receiver', None)

class array_check(object):
    def __init__(self):
        self.last_time = None
        
    def __call__(self, fr):
        if fr.type != core.G3FrameType.GcpSlow:
            return
        if self.last_time is None:
            self.last_time = fr['array']['frame']['utc']
            return
        this_time = fr['array']['frame']['utc']
        bad = (this_time.time - self.last_time.time) / core.G3Units.s != 1.0
        self.last_time = fr['array']['frame']['utc']
        return bad

if __name__ == '__main__':
    time.sleep(int(np.random.rand() * 100))
    arcfile = sys.argv[1]
    outdir = sys.argv[2]
    bn = os.path.basename(arcfile).split('.')[0]
    
    scu = scu_check()
    pipe = core.G3Pipeline()
    pipe.Add(gcp.ARCFileReader, filename = arcfile)
    pipe.Add(drop_receiver)
    pipe.Add(scu)
    pipe.Add(FrameGrabber, condition = array_check(), extra_frames = 1)
    pipe.Add(core.G3Writer, filename = os.path.join(outdir, bn + '.g3'))
    pipe.Run()
    
    f = open(os.path.join(outdir, bn + '_scu.txt'), 'w')
    if scu.outages > 0:
        f.write('number:\t {:d}\n'.format(scu.outages))
        f.write('time:\t{:.02f}\n'.format(scu.outage_time))
    f.close()
    
