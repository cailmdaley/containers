import os
import argparse
from collections import deque
import numpy as np
from spt3g import core, gcp
from frame_grabber import FrameGrabber, FrameGrabberWriter
from whathaveyoudonetotime import TimeTest

class ACUStatusTest(object):
    def __init__(self):
        self.last_controls = deque([], 2)

    def _check_control(self, fr):
        if len(self.last_controls) < 2:
            self.last_controls.append(fr['antenna0']['tracker']['inControl'])
            return False
        mean_cont = np.mean([self.last_controls[0], self.last_controls[1]])
        return mean_cont > .75
        
    def __call__(self, fr):
        if fr.type == core.G3FrameType.EndProcessing:
            return
        acu_status = fr['antenna0']['acu']['acu_status'].value
        return (self._check_control(fr) and acu_status != 224 and acu_status != 192)

def grab_receiver_stuff(fr):
    if fr.type == core.G3FrameType.EndProcessing:
        return
    fr['receiver_frame'] = fr['receiver']['frame']
    fr['receiver_utc'] = fr['receiver']['bolometers']['utc'][0]
    del fr['receiver']

def getfilename(fr, dir):
    try:
        return os.path.join(dir, 
                            fr['antenna0']['tracker']['utc'][0][0].GetFileFormatString() + 
                            '_stripped.g3.gz')
    except:
        return

def time_update(fr, hours = 2):
    if fr.type == core.G3FrameType.EndProcessing:
        return
    if fr['array']['frame']['utc'].time / core.G3Units.s % (60 * 60 * hours) == 0:
        print fr['array']['frame']['utc']

def do(arcfiles, time_dir, acu_dir):
    time_ = TimeTest()
    acu_ = ACUStatusTest()
    acu_grabber = FrameGrabber(acu_, 300)
    time_grabberwriter = FrameGrabberWriter(time_, lambda fr: getfilename(fr, time_dir)
                                            , 30)
    pipe = core.G3Pipeline()
    pipe.Add(gcp.ARCFileReader(arcfiles))
    pipe.Add(grab_receiver_stuff)
    pipe.Add(time_update, hours = .5)
    pipe.Add(time_grabberwriter)
    pipe.Add(acu_grabber)
    pipe.Add(core.G3MultiFileWriter, filename = lambda fr, seq: getfilename(fr, acu_dir),
             size_limit = 1024**6, divide_on = acu_grabber.issplit)
    pipe.Run()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', type = str, default = '/data31/ndhuang/acu_timing_badness')
    parser.add_argument('arcfiles', nargs = '+')
    args = parser.parse_args()
    do(args.arcfiles, os.path.join(args.directory, 'time'), os.path.join(args.directory, 'acu'))

