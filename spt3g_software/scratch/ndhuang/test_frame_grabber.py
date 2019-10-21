'''
A happy little testing script.
'''

from collections import deque
from spt3g import core
from frame_grabber import FrameGrabber
class TestFrames(object):
    def __init__(self, extras = 1):
        self.count = 0
        self.extras = extras

    def __call__(self, fr,):
        fr = core.G3Frame(core.G3FrameType.Timepoint)
        if self.count < 3:
            out = 1
        elif self.count < 6:
            out = 0
        elif self.count < 8:
            out = 1
        elif self.count < 12:
            out = 0
        elif self.count < 12 + self.extras:
            out = 1
        else:
            return []
        fr['test'] = core.G3Int(out)
        self.count += 1
        return fr

class Accumulator():
    def __init__(self):
        self.d = deque()
        self.breaks = -1
    def __call__(self, fr):
        if fr.type == core.G3FrameType.EndProcessing:
            return
        if FrameGrabber.issplit(fr):
            self.d.append(deque())
            self.breaks += 1
        self.d[self.breaks].append(fr['test'])
        return

def print_(fr):
    if fr.type == core.G3FrameType.EndProcessing:
        return
    print fr['test']

def print_fr(fr):
    print fr

def test():
    grabber = FrameGrabber(lambda fr: fr['test'] == 0, 1)
    frames = TestFrames(1)
    acc = Accumulator()
    pipe = core.G3Pipeline()
    pipe.Add(frames)
    pipe.Add(grabber)
    pipe.Add(acc)
    pipe.Run()
    assert(len(acc.d) == 2)
    assert(list(acc.d[0]) == [1, 0, 0, 0, 1])
    assert(list(acc.d[1]) == [1, 0, 0, 0, 0, 1])

    grabber = FrameGrabber(lambda fr: fr['test'] == 0, 4)
    frames = TestFrames(1)
    acc = Accumulator()
    pipe = core.G3Pipeline()
    pipe.Add(frames)
    pipe.Add(grabber)
    pipe.Add(acc)
    pipe.Run()
    assert(len(acc.d) == 1)
    assert(list(acc.d[0]) == [1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1])

    grabber = FrameGrabber(lambda fr: fr['test'] == 0, 4)
    frames = TestFrames(5)
    acc = Accumulator()
    pipe = core.G3Pipeline()
    pipe.Add(frames)
    pipe.Add(grabber)
    pipe.Add(acc)
    pipe.Run()
    assert(len(acc.d) == 1)
    assert(list(acc.d[0]) == [1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1])

test()
