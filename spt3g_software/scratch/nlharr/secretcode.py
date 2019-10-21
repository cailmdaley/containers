from spt3g import core

class FrameGrabber(object):
    def __init__(self):
        self.frames = []
    def __call__(self, frame):
        self.frames.append(frame)


#example code
fg = FrameGrabber()

pipe = core.G3Pipeline()
pipe.Add(AllYourProcessing)
pipe.Add(fg)
pipe.Run()

#your frames now live in fg.frames
