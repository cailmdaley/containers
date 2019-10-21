import argparse
import numpy as np
from spt3g import core, gcp
from frame_grabber import FrameGrabber

def tracker_check(fr):
    dt = np.diff(np.array(fr['antenna0']['tracker']['utc']).astype(float))
    return (dt != .01 * core.G3Units.s).any()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('arcfiles', nargs = '+')
    parser.add_argument('-o', '--output', default = 'output-%04u.g3')
    args = parser.parse_args()
    grabber = FrameGrabber(tracker_check, extra_frames = 8)
    pipe = core.G3Pipeline()
    pipe.Add(gcp.ARCFileReader, filename = args.arcfiles)
    pipe.Add(gcp.ARCExtract)
    pipe.Add(grabber)
    pipe.Add(core.G3MultiFileWriter, filename = args.output, 
             size_limit = 1024 * 1024 * 1024,
             divide_on = grabber.issplit)
    pipe.Run()
