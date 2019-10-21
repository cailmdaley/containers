from spt3g import core

def CopyKey(frame, key_src, key_dest, frame_type):
    if frame.type == frame_type:
        frame[key_dest] = frame[key_src]

def FilterFramesOfTypes(frame, frame_type, types_to_remove, keys_to_save = None):
    if keys_to_save == None:
        keys_to_save = []
    if frame.type == frame_type:
        for k in frame.keys():
            if type(frame[k]) in types_to_remove and (not k in keys_to_save):
                del frame[k]
@core.indexmod
class FrameExtractor(object):
    def __init__(self, frame_type = None):
        self.frames = []
        self._frame_type = frame_type
    def __call__(self, frame):
        if (frame.type == self._frame_type or frame.type == None):
            self.frames.append(frame)


def PrintScanNumber(frame):
    if 'ScanNumber' in frame:
        print('scan_number: %d' %int(frame['ScanNumber']))



def JustSaveSomeKeys(frame, ks, frame_type):
    if frame.type == frame_type:
        new_frame = core.G3Frame(frame.type)
        for k in ks:
            new_frame[k] = frame[k]
        return [new_frame]
