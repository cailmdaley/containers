from spt3g import core, gcp
import os, math
import numpy as np

def ARCForTimerange(start, stop, basedir = '/spt/data/arc',
                    no_more_data_error=True):
    '''
    Find archive files for a given time range in `dir`.
    The time range of the returned arcfiles is guaranteed to span both `start`
    `stop`, but may not be contiguous iff `dir` does not contain contiguous 
    arcfiles.
    '''
    if not isinstance(start, core.G3Time):
        start = core.G3Time(start)
    if not isinstance(stop, core.G3Time):
        stop = core.G3Time(stop)
    
    dirs = os.listdir(basedir)
    if 'this_is_an_arcfile_summary_dir' in dirs:
        summary_dir = basedir
        t_dayonly = t.split('_')
        t_dayonly[1] = '000000'
        t_dayonly = '_'.join(t_dayonly)
        if t_dayonly in dirs:
            dirs.sort()
            dir = os.path.join(basedir, dirs[dirs.index(t_dayonly) - 2])
        else:
            dirs.append(t)
            dirs.remove('this_is_an_arcfile_summary_dir')
            dirs.sort()
            dir = os.path.join(basedir, dirs[dirs.index(t) - 1])
        files = []
        out = []
        i_dir = 0
        while out is not None:
            out = ARCForTimeRange(start, stop, dir, False)
            i_dir += 1
            dir = os.path.join(basedir, dirs[dirs.index(t) - 1 + i_dir])
            files += out
        return out
    else:
        summary_dir = None
        dir = basedir
        files = dirs
        files = [f for f in files if f.endswith('.dat') or f.endswith('.dat.gz')]
        if len(files) == 0:
            if no_more_data_error:
                raise FileNotFoundError('No ARC files in %s' % dir)
            else:
                return None
        files = sorted(files)
        i_start = np.searchsorted(files, start.GetFileFormatString() + '.dat',
                                  side = 'right') - 1
        i_stop = np.searchsorted(files, stop.GetFileFormatString() + '.dat',
                                 side = 'left') + 1
        files = [os.path.join(dir, f) for f in files[i_start:i_stop]]
        return files
    

def ARCForTime(t, basedir):
    if isinstance(t, core.G3Time):
        t = t.GetFileFormatString()
    elif isinstance(t, str):
        t = core.G3Time(t)
        t = t.GetFileFormatString()
    
    dirs = os.listdir(basedir)
    if 'this_is_an_arcfile_summary_dir' in dirs:
        t_dayonly = t.split('_')
        t_dayonly[1] = '000000'
        t_dayonly = '_'.join(t_dayonly)
        if t_dayonly in dirs:
            dirs.sort()
            dir = os.path.join(basedir, dirs[dirs.index(t_dayonly) - 2])
        else:
            dirs.append(t)
            dirs.remove('this_is_an_arcfile_summary_dir')
            dirs.sort()
            dir = os.path.join(basedir, dirs[dirs.index(t) - 1])
        files = os.listdir(dir)
    else:
        dir = basedir
        files = dirs

    files = [f for f in files if f.endswith('.dat') or f.endswith('.dat.gz')]
    if len(files) == 0:
        raise FileNotFoundError('No ARC files in %s' % dir)
    if t + '.dat' in files:
        return os.path.join(dir, t + '.dat')
    elif t + '.dat.gz' in files:
        return os.path.join(dir, t + '.dat.gz')
    else:
        files.append(t + '.dat')
        files.sort()
        return os.path.join(dir, files[files.index(t + '.dat') - 1])

def NextARCFile(path, no_more_data_error=True):
    '''
    Path corresponding to the ARC file following this one in time.
    '''
    arcdir = os.path.dirname(path)
    arcfiles = os.listdir(arcdir)
    arcfiles.sort()
    arcfiles = [f for f in arcfiles if f.endswith('.dat') or f.endswith('.dat.gz')]
    i = arcfiles.index(os.path.basename(path))
    if i+1 < len(arcfiles):
        return os.path.join(arcdir, arcfiles[i+1])
    else:
        arcparentdir = os.path.dirname(arcdir)
        dirs = os.listdir(arcparentdir)
        if 'this_is_an_arcfile_summary_dir' in dirs:
            dirs.remove('this_is_an_arcfile_summary_dir')
        else:
            if no_more_data_error:
                core.log_fatal('No more data in ARC directory.', unit='ARCTimerangeReader')
            else:
                return None
        dirs.sort()
        i = dirs.index(os.path.basename(arcdir))
        if i+1 >= len(arcfiles):
            if no_more_data_error:
                core.log_fatal('No more data in ARC directory.', unit='ARCTimerangeReader')
            else:
                return None
        arcdir = os.path.join(arcparentdir, dirs[i+1])
        arcfiles = os.listdir(arcdir)
        arcfiles.sort()
        try:
            # Some dirs have duplicates. Who knows why?
            i = arcfiles.index(os.path.basename(path)) + 1
        except (IndexError, ValueError):
            i = 0
        return os.path.join(arcdir, arcfiles[i])

@core.indexmod
class ARCTimerangeReader(object):
    '''
    Read data from an ARC directory between a given start and stop time.
    '''
    def __init__(self, start_time, stop_time, basedir='/spt/data/arc',
                 no_more_data_error=True):
        self.arcfiles = ARCForTimerange(start_time, stop_time, basedir = basedir)
        self.cur_file = self.arcfiles.pop(0)
        core.log_info('Opening %s' % self.cur_file, unit='ARCTimerangeReader')
        self.reader = gcp.ARCFileReader(filename=self.cur_file)
        if isinstance(stop_time, str):
            self.stop_time = core.G3Time(stop_time)
        else:
            self.stop_time = stop_time
        if isinstance(start_time, str):
            self.start_time = core.G3Time(start_time)
        else:
            self.start_time = start_time
        self.no_more_data_error = no_more_data_error
    def __call__(self, frame):
        t = core.G3Time(0) # 1970 predates SPT
        assert(frame is None)
        while t < self.start_time:
            nextframes = self.reader(None)
            if len(nextframes) == 0:
                if len(self.arcfiles) == 0:
                    return []
                self.cur_file = self.arcfiles.pop(0)
                core.log_info('Opening %s' % self.cur_file, unit='ARCTimerangeReader')
                self.reader = gcp.ARCFileReader(filename=self.cur_file)
                continue
    
            assert(len(nextframes) == 1)
            t = nextframes[0]['array']['frame']['utc']
        if self.stop_time is not None and t > self.stop_time:
            return []
        else:
            return nextframes

@core.indexmod
class ARCInterposer(object):
    '''
    Read data from an ARC directory to provide metadata for an existing set
    of time point frames. GcpSlow frames from the ARC files will be
    interleaved with the timepoint frames in the original stream. Any
    non-timepoint frames (wiring frames, for example) will be passed through
    silently to the next module.
    '''
    def __init__(self, basedir='/data/sptdat/arcfile_directory'):
        self.basedir = basedir
        self.cur_file = None
        self.reader = None
        self.gcp_time = core.G3Time(0) # 1970 predates SPT
    def __call__(self, frame):
        if frame.type != core.G3FrameType.Timepoint:
            return [frame]

        t = frame['EventHeader']
        if math.floor(self.gcp_time.time / core.G3Units.s) == math.floor(t.time / core.G3Units.s): # GCP slow frames are 1 second long at second boundaries
            return [frame]

        # Check that the while loop below will actually execute at least once
        if math.floor(self.gcp_time.time / core.G3Units.s) > math.floor(t.time / core.G3Units.s):
            core.log_fatal('Time has gone backwards: %s > %s' % (self.gcp_time, t), unit='ARCInterposer')

        # Loop through the file and find the GCP frame for the second including this sample
        while math.floor(self.gcp_time.time / core.G3Units.s) < math.floor(t.time / core.G3Units.s):
            if self.reader is None:
                new_file = ARCForTime(t, self.basedir)
                if self.cur_file == new_file:
                    # Correct for roundoff; avoid loops
                    new_file = NextARCFile(new_file)
                self.cur_file = new_file
                core.log_info('Opening %s for data at %s' % (self.cur_file, t), unit='ARCTimerangeReader')
                self.reader = gcp.ARCFileReader(filename=self.cur_file)
            nextframes = self.reader(None)
            if len(nextframes) == 0:
                self.reader = None
                continue # This will trigger the code above
    
            assert(len(nextframes) == 1)
            # The "array time" in GCP frames is the time at which
            # the frame *ends*, not when it starts, so subtract
            # one second
            self.gcp_time = core.G3Time(nextframes[0]['array']['frame']['utc'].time - core.G3Units.s)
        return [nextframes[0], frame]


