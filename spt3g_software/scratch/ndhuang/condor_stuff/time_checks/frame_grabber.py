from collections import deque
from spt3g import core

@core.indexmod
class FrameGrabber(object):
    '''
    Grab contiguous sets of frames that meet `condition` (and possibly frames before and after).
    This module adds a key 'FRAMEGRABBER_SPLIT' to the first frame of each contiguous section.  
    If you use FrameGrabber.issplit to detect the splits, the key is deleted.

    INPUTS
    ------
    condition: callable
        The condition under which we should return this frame.  
        Returns True if the frame should be returned.
        Uses the following signature:
        condtion(frame, **kwargs)
    extra_frames: int, optional [0]
        The number of frames before and after frames meeting `condition` to return.
        For example, if `extra_frames` is set to 3, and frames 109 and 110 meet `condition`,
        frames 106-113 will be returned.
    condition_args: dict, optional [{}]
        Extra keyword arguments to `condition`
    '''
    def __init__(self, condition = None, extra_frames = 0, **condition_args):
        assert(condition is not None)
        self.condition = condition
        self.extra_frames = extra_frames
        self.last_frames = deque([], extra_frames + 1)
        self.since_bad = extra_frames + 1
        self.extra_args = condition_args

    def __call__(self, fr):
        if fr.type == core.G3FrameType.EndProcessing:
            return

        self.last_frames.append(fr)
        if self.condition(fr, **self.extra_args):
            if self.since_bad <= self.extra_frames:
                # grab all the frames since the last time we had a bad one
                frames_to_save = list(self.last_frames)[-(self.since_bad + 1):]
            else:
                frames_to_save = list(self.last_frames)
                frames_to_save[0]['FRAMEGRABBER_SPLIT'] = True
            self.since_bad = 0
        else:
            if self.since_bad < self.extra_frames:
                frames_to_save = fr
            else:
                frames_to_save = False
            self.since_bad += 1
        return frames_to_save

    @staticmethod
    def issplit(fr):
        '''
        Test this frame to see if we should split on it.  
        Useful for writing contiguous frames to disk.
        '''
        split = 'FRAMEGRABBER_SPLIT' in fr
        if split:
            del fr['FRAMEGRABBER_SPLIT']
        return split

class FrameGrabberWriter(object):
    '''
    Grab contiguous sets of frames that meet `condition` (and possibly frames 
    before and after) and write them to disk.  All frames are passed to the 
    next module.

    INPUTS
    ------
    condition: callable
        The condition under which we should return this frame.  
        Returns True if the frame should be returned.
        Uses the following signature:
        condtion(frame, **kwargs)
    extra_frames: int, optional [0]
        The number of frames before and after frames meeting `condition` to return.
        For example, if `extra_frames` is set to 3, and frames 109 and 110 
        meet `condition`, frames 106-113 will be returned.
    filenamer: callable
        A callable that takes a frame as its argument a returns a string, which will be 
        the new filename.
    condition_args: dict, optional [{}]
        Extra keyword arguments to `condition`
    '''
    def __init__(self, condition, filenamer, extra_frames = 0, **condition_args):
        self.condition = condition
        self.extra_frames = extra_frames
        self.last_frames = deque([], extra_frames + 1)
        self.since_bad = extra_frames + 1
        self.extra_args = condition_args
        self.filenamer = filenamer
        self.writer = None

    def _newfile(self, fr):
        if self.writer is not None:
            # close the last writer
            self.writer(core.G3Frame(core.G3FrameType.EndProcessing))
        self.writer = core.G3Writer(self.filenamer(fr))

    def __call__(self, fr):
        if fr.type == core.G3FrameType.EndProcessing:
            if self.writer is not None:
                self.writer(fr)
            return

        self.last_frames.append(fr)
        if self.condition(fr, **self.extra_args):
            if self.since_bad <= self.extra_frames:
                # grab all the frames since the last time we had a bad one
                frames_to_save = list(self.last_frames)[-(self.since_bad + 1):]
            else:
                frames_to_save = list(self.last_frames)
                self._newfile(frames_to_save[0])
            self.since_bad = 0
        else:
            if self.since_bad < self.extra_frames:
                frames_to_save = fr
            else:
                frames_to_save = False
            self.since_bad += 1
        if frames_to_save:
            if not isinstance(frames_to_save, core.G3Frame):
                for fr in frames_to_save:
                    self.writer(fr)
            else:
                self.writer(frames_to_save)
