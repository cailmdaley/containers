from spt3g import core
import cPickle as pickle



class GenFrame(object):
    def __init__(self):
        self.sent = False
    def __call__(self, frame):
        if self.sent:
            return []
        test_frame = core.G3Frame(core.G3FrameType.Wiring)
        test_frame['three'] = 3
        test_frame['four'] = 3
        test_frame['five'] = 'five'
        pickle.dump(test_frame, open("test_frame.pkl", 'w'), protocol = 2)
        self.sent = True
        return test_frame

pipe = core.G3Pipeline()

pipe.Add(GenFrame)
pipe.Add(core.G3Writer, filename = 'test_frame.g3')
pipe.Run()
