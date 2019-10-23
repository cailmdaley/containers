import argparse
import numpy as np
from spt3g import core, gcp
from frame_grabber import FrameGrabber

def acu_check(fr):
    dt = np.diff(fr['antenna0']['tracker']['utc'].astype(float))
    return (dt != .01 * core.G3Unit.s).any()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('arcfiles', nargs = '+')
    parser.add_argument('-o', '--output', default = 'output-%04u.g3')
    pipe = core.G3Pipeline()
    pipe.Add(gcp.ARCFileReader, filename = arcfile)
