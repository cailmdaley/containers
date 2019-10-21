from .GCPScanFlags import InsertGCPScanBoundaries, MoveGCPScanDataToScanFrame, ConsolidateHousekeepingToScanFrame, LimitScanLengths
from .opticalbench import TrimBenchPosition

from spt3g import core, dfmux



# from .pointing import NaiveBoresightPointing
from spt3g.pointing.pointing_for_mapmaking import NaiveBoresightPointing



@core.pipesegment
def BuildScanFramesFromRawData(pipe, flac_compress=True, drop_fast_frames=True, max_scan_length=100000):
    '''
    Takes a GcpSlow+Timepoint stream of data and converts it into scan
    frames. gcp.ARCExtract must have been run on the data first.
    '''
    pipe.Add(InsertGCPScanBoundaries)
    # Limit scan lengths to avoid using all the memory during noise stares
    if max_scan_length != 0:
        pipe.Add(LimitScanLengths, N=max_scan_length)
    pipe.Add(MoveGCPScanDataToScanFrame)
    pipe.Add(lambda f: f.type != core.G3FrameType.GcpSlow)
    pipe.Add(dfmux.DfMuxCollator, flac_compress=flac_compress, drop_timepoints=drop_fast_frames)
    pipe.Add(ConsolidateHousekeepingToScanFrame, drop_hk_frames=drop_fast_frames)
    pipe.Add(NaiveBoresightPointing)
    pipe.Add(TrimBenchPosition)
    
@core.pipesegment
def BuildScanFramesFromARCFile(pipe, flac_compress=True, drop_fast_frames=True, max_scan_length=100000):
    '''
    Takes a GcpSlow stream of data and converts it into scan
    frames. gcp.ARCExtract must have been run on the data first.
    '''
    def fake_timepoints(frame):
        if frame.type != core.G3FrameType.GcpSlow:
            return
        tp = core.G3Frame(core.G3FrameType.Timepoint)
        tp['EventHeader'] = core.G3Time(frame['antenna0']['tracker']['utc'][0][0])
        return [tp, frame]
    pipe.Add(fake_timepoints)
    pipe.Add(InsertGCPScanBoundaries)
    # Limit scan lengths to avoid using all the memory during noise stares
    if max_scan_length != 0:
        pipe.Add(LimitScanLengths, N=max_scan_length)
    pipe.Add(MoveGCPScanDataToScanFrame, drop_scans_with_no_gcp_data = False)
    pipe.Add(lambda f: f.type != core.G3FrameType.GcpSlow)
    pipe.Add(lambda f: f.type != core.G3FrameType.Timepoint)
    pipe.Add(NaiveBoresightPointing, bolots = None)
    pipe.Add(TrimBenchPosition, trim_key = 'RawBoresightAz')
    

