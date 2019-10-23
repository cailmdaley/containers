import os
import argparse
from spt3g import core, gcp
from frame_grabber import FrameGrabber

class TimeTest(object):
    def __init__(self):
        self.last_time = None

    def __call__(self, fr):
        if fr.type == core.G3FrameType.EndProcessing:
            return 
        time = fr['antenna0']['tracker']['utc'][0][0].time / core.G3Units.s
        if self.last_time is None:
            self.last_time = time
            return False
        dt = time - self.last_time
        self.last_time = time
        if abs(dt - 1.0) > 0:
            return True

def filename(fr):
    try:
        return fr['antenna0']['tracker']['utc'][0][0].GetFileFormatString() + '_stripped.g3.gz'
    except:
        return

def strip_bolodata(fr):
    try:
        del fr['receiver']
        return fr
    except:
        return

def time_update(fr, hours = 2):
    if fr.type == core.G3FrameType.EndProcessing:
        return
    if fr['array']['frame']['utc'].time / core.G3Units.s % (60 * 60 * hours) == 0:
        print fr['array']['frame']['utc']

def file_update(fr):
    try:
        if "FRAMEGRABBER_SPLIT" in fr:
            print filename(fr)
    except:
        return

def do(arcfiles, directory, extra_frames):
    condition = TimeTest()
    grabber = FrameGrabber(condition, extra_frames)
    pipe = core.G3Pipeline()
    pipe.Add(gcp.ARCFileReader(arcfiles))
    pipe.Add(strip_bolodata)
    pipe.Add(time_update, hours = 2)
    pipe.Add(grabber)
    pipe.Add(file_update)
    pipe.Add(core.G3MultiFileWriter, filename = lambda fr, seq: os.path.join(directory, filename(fr)), 
             size_limit = 1024**6, divide_on = grabber.issplit)
    pipe.Run()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', type = str, default = '/data31/ndhuang/acu_timing_badness')
    parser.add_argument('arcfiles', nargs = '+')
    args = parser.parse_args()
    do(args.arcfiles, args.directory, 30)
