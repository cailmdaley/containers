import numpy, sys
if sys.version_info[0] >= 3:
    from functools import reduce
from spt3g import core, gcp, dfmux

@core.indexmod
class InsertGCPScanBoundaries(object):
    '''
    Inserts new scan frames when the GCP scan status, observation ID, or
    field change. This will create scan frames for both the main part of
    the scan and the turnaround. Meant to be used before DfMuxCollator.
    The turnaround scan will be marked with a "Turnaround" key in the
    scan frame set to True.
    '''
    def __init__(self):
        self.last_flag = None
        self.last_obs = None
        self.last_feature_bits = None
        self.dropout = False
        self.scan_starts = []
    def __call__(self, frame):
        if frame.type == core.G3FrameType.GcpSlow:
            # When we get GCP slow frames, watch the tracker
            # scan_flag field. When it flips, start a scan.
            # False means a turnaround. True means the telescope
            # is going.

            # Also make sure to break if the ObservationID or
            # SourceName change.

            # The slow frame will be emitted immediately following
            # the scan start frame. For scans that start part way
            # through the slow frame, the slow frame will be
            # duplicated and occur in both the old and new scan
            # to ease pointing interpolation.

            flags = frame['TrackerStatus'].scan_flag
            times = frame['TrackerStatus'].time

            # Check if obs has changed, start a new scan if so
            source = frame['SourceName']
            if 'obs_id' in frame['antenna0']['tracker']:
                obsid = frame['antenna0']['tracker']['obs_id'].value
            else:
                obsid = -1
            # If the antenna process hiccups, don't break the scan, because
            # we might be able to use the data anyway.
            if obsid == 0:
                obsid = self.last_obs[1]
                self.dropout = True
            else:
                self.dropout = False
            obs = (source, obsid)
            if obs != self.last_obs:
                self.last_flag = None
                self.last_obs = obs

            # The same if feature bits have changed
            feature_bits = list(frame['GCPFeatureBits'])
            if feature_bits != self.last_feature_bits:
                self.last_flag = None
                self.last_feature_bits = feature_bits
            
            # Find indices (if any) where the scan flag changed
            changes = []
            if self.last_flag != flags[0]:
                changes = [-1]
            changes += list(numpy.where(
              numpy.diff(numpy.asarray(flags, dtype=int)) != 0)[0])
            self.last_flag = flags[-1]

            for change in changes:
                self.scan_starts.append((times[int(change+1)],
                  flags[int(change+1)], core.G3Frame(frame)))
            if len(changes) == 0 and len(self.scan_starts) > 0:
                # hokay, so...
                # This only triggers if this frame had no changes, but there are
                # still cached changes.  This only happens if we have dropped Timepoint
                # frames because they are unneeded, and have old data hanging around.
                # We should, therefore, drop all the cached data, replace it with
                # this frame, and move on.
                self.scan_starts = [(times[0], flags[0], core.G3Frame(frame))]
                return []
            elif len(changes) == 0 or changes[0] != -1:
                return [frame]
            else:
                # If the first sample is already in the next
                # scan, the old scan does not need a copy.
                return []

        if frame.type == core.G3FrameType.Timepoint:
            t = frame['EventHeader']
            out = []
            while len(self.scan_starts) > 0 and t >= self.scan_starts[0][0]:
                scan = core.G3Frame(core.G3FrameType.Scan)
                if not self.scan_starts[0][1]:
                    scan['Turnaround'] = True
                if self.dropout:
                    scan['GCPDropout'] = True
                out = [scan, self.scan_starts[0][2]]
                del self.scan_starts[0]
            out.append(frame)
            return out

        return [frame]

@core.indexmod
class LimitScanLengths(object):
    '''
    Insert synthetic scan frames every N timepoint frames in the absence
    of scan frames coming from another source. This has the effect of
    preventing the synthesis of scan frames by DfMuxCollator with more than
    N timepoints. Scan frames created by this module with have a boolean
    value named "Synthetic" set to True.
    '''
    def __init__(self, N=100000):
        self.N = N
        self.ticker = 0
    def __call__(self, frame):
        if frame.type == core.G3FrameType.Scan:
            self.ticker = 0
        if frame.type == core.G3FrameType.Timepoint:
            self.ticker += 1
            if self.ticker >= self.N:
                out = core.G3Frame(core.G3FrameType.Scan)
                out['Synthetic'] = True
                self.ticker = 0
                return [out, frame]

@core.indexmod
class MoveGCPScanDataToScanFrame(object):
    '''
    Moves a standard set of data from GCP slow frames into Scan frames.
    '''
    def __init__(self, drop_scans_with_no_gcp_data=True):
        self.buffer = []
        self.drop_no_gcp = drop_scans_with_no_gcp_data
    def __call__(self, frame):
        out = []

        if len(self.buffer) > 0 and (frame.type == core.G3FrameType.Scan or frame.type == core.G3FrameType.EndProcessing):
            gcp_frames = [f for f in self.buffer if f.type == core.G3FrameType.GcpSlow]
            scan = self.buffer[0]
            assert(scan.type == core.G3FrameType.Scan)
            if len(gcp_frames) == 0 and self.drop_no_gcp:
                # Nonsense scans with no real data can happen at script termination sometimes
                self.buffer = []
                if frame.type == core.G3FrameType.EndProcessing:
                    return
                else:
                    return []
            scan['SourceName'] = gcp_frames[0]['SourceName']
            scan['GCPFeatureBits'] = gcp_frames[0]['GCPFeatureBits']
            scan['TrackerStatus'] = reduce(lambda a,b: a+b, [f['TrackerStatus'] for f in gcp_frames])
            scan['TrackerPointing'] = reduce(lambda a,b: a+b, [f['TrackerPointing'] for f in gcp_frames])
            scan['OnlinePointingModel'] = gcp_frames[0]['OnlinePointingModel']
            scan['ACUStatus'] = gcp.ACUStatusVector([f['ACUStatus'] for f in gcp_frames])
            if 'ObservationID' in gcp_frames[0]:
                scan['ObservationID'] = gcp_frames[0]['ObservationID']
            # Concatenate the optical bench timestreams
            # Allow up to 10s gaps between frames
            benchpos = core.G3TimestreamMap.concatenate([f['BenchPosition'] for f in gcp_frames],
                                                        ts_interp_threshold=1000)
            benchcom = core.G3TimestreamMap.concatenate([f['BenchCommandedPosition'] for f in gcp_frames],
                                                        ts_interp_threshold=1000)
            benchzero = core.G3TimestreamMap.concatenate([f['BenchZeros'] for f in gcp_frames],
                                                         ts_interp_threshold=1000)
            scan['BenchPosition'] = benchpos
            scan['BenchCommandedPosition'] = benchcom
            scan['BenchZeros'] = benchzero
            out = self.buffer

        if frame.type == core.G3FrameType.Scan:
            self.buffer = [frame]
            # See out settings above
        elif frame.type == core.G3FrameType.EndProcessing:
            out.append(frame)
            self.buffer = []
        else:
            if len(self.buffer) == 0:    
                # Just pass frames before the first scan
                out = [frame]
            else:
                self.buffer.append(frame)
            
        return out

@core.indexmod
class ConsolidateHousekeepingToScanFrame(object):
    '''
    Move housekeeping data to Scan frames, optionally deleting the original Housekeeping frames. Analog to DfMuxCollator.
    '''

    def __init__(self, drop_hk_frames = True, tolerance = 20*core.G3Units.s):
        '''
        drop_hk_frames [True]: if True, delete housekeeping frames once merged
    
        tolerance [20 seconds]: The official housekeeping for a scan is the first housekeeping frame that arrives between
        the start time of the scan and tolerance later, or the most recent prior housekeeping data if none arrives in that
        window.
        '''

        self.drop_hk_frames = drop_hk_frames
        self.tolerance = tolerance
        self.hk = {}
        self.buffer = []
        self.waiting = False
        self.scan_start = None

    def FlushQueue(self, out):
        hkdat = dfmux.DfMuxHousekeepingMap()
        for k in self.hk.keys():
            hkdat[k] = self.hk[k]
        self.buffer[0]['DfMuxHousekeeping'] = hkdat
        out += self.buffer
        self.buffer = []
        self.waiting = False

    def __call__(self, frame):
        out = []
        if frame.type == core.G3FrameType.Housekeeping:
            timestamp = frame['DfMuxHousekeeping'].values()[0].timestamp

            if self.waiting and timestamp.time > self.scan_start + self.tolerance: 
                # We were waiting for HK data, but this is newer than we wanted.
                # Flush the old stuff and move on.
                self.FlushQueue(out)

            for k in frame['DfMuxHousekeeping'].keys():
                self.hk[k] = frame['DfMuxHousekeeping'][k]
            if self.waiting:
                # We were waiting for HK data and this is it. Flush the queue.
                self.FlushQueue(out)

            if not self.drop_hk_frames: # NB: we can't get here and still be queueing
                out.append(frame)
            return out

        if frame.type == core.G3FrameType.EndProcessing:
            # Last call. Use whatever we have.
            if self.waiting:
                self.FlushQueue(out)
            out.append(frame)
            return out

        if frame.type == core.G3FrameType.Scan:
            # Last call. Use whatever we have.
            if self.waiting:
                self.FlushQueue(out)
            assert(len(self.buffer) == 0)
            self.buffer.append(frame)
            self.scan_start = frame['RawTimestreams_I'].start.time
            self.waiting = True
            return out
        
        if self.waiting:
            self.buffer.append(frame)
        else:
            out.append(frame)
        return out
        
