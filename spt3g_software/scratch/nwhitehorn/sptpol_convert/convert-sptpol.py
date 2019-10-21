import numpy, sys, os
from spt3g import core, gcp, std_processing, dfmux, coordinateutils, sptpol
import spt3g.pointing.offline_pointing as op
from spt3g.coordinateutils.coordsysmodules import FillCoordTransRotations

# Usage: convert-sptpol-data.py start-time stop-time output-dir input-files
# Converts SPTpol data read from ARC files into data that is supposed to look
# like SPT3G data fresh off the satellite.

#core.set_log_level(core.G3LogLevel.LOG_DEBUG)

# Make numpy problems fail loudly
#numpy.seterr(all='raise')

# In case you get warnings from astropy about timekeeping:
from astropy.utils import iers
#iers.IERS.iers_table = iers.IERS_A.open(iers.IERS_A_URL)

pipe = core.G3Pipeline()
print(sys.argv[4])
pipe.Add(gcp.ARCFileReader, filename=sys.argv[4:])

# For 100d data
#pipe.Add(core.G3Reader, filename='sptpol100d-wiring.g3')

start_t = core.G3Time(sys.argv[1])
stop_t = core.G3Time(sys.argv[2])
pipe.Add(lambda fr: fr.type != core.G3FrameType.GcpSlow or (fr['array']['frame']['utc'] >= start_t and fr['array']['frame']['utc'] < stop_t))

pipe.Add(gcp.ARCExtract)
pipe.Add(gcp.UnpackSPTpolHKData)
pipe.Add(gcp.GCPMuxDataDecoder)
pipe.Add(std_processing.BuildScanFramesFromRawData, flac_compress=True)

# Drop scans with no set feature bits
pipe.Add(lambda fr: fr.type != core.G3FrameType.Scan or
  'analyze' in fr['GCPFeatureBits'])

# Drop 1-sample scans that can occur from the interaction of the cut in
# DropOrphanSamples and the 10 ms offset between the time of the first
# tracker sample and the notional time of the frame. These will later
# have infinite sample rates and are guaranteed to be thrown away anyway
# since these scans occur only immediately before we start scanning.
pipe.Add(lambda fr: 'RawBoresightAz' not in fr or len(fr['RawBoresightAz']) > 1)

# Drop metadata that does not apply to any scans
pipe.Add(core.DropOrphanMetadata)

# Pretend f6 is set on RCW38 fast points, at least for the purpose of
# setting the source name. Guess what is and isn't a fast point by checking
# if the observation lasted more than 30 minutes for lack of better
# options.
def fabricate_f6_source_name(frame):
    if 'SourceName' in frame and frame['SourceName'] == 'RCW38' and (stop_t.time - start_t.time) > 30*core.G3Units.minute:
        del frame['SourceName']
        frame['SourceName'] = 'RCW38-pixelraster'
pipe.Add(fabricate_f6_source_name)

# Add Observation frames
class fabricate_obsframe(object):
    def __init__(self):
        self.obsframe_out = False
        self.buffer = []
    def __call__(self, frame):
        if frame.type == core.G3FrameType.EndProcessing:
            return
        if self.obsframe_out:
            frame['ObservationID'] = self.obsid
            return
        if 'SourceName' not in frame:
            self.buffer.append(frame)
            return []
        f = core.G3Frame(core.G3FrameType.Observation)
        f['SourceName'] = frame['SourceName']
        start = core.G3Time(sys.argv[1])
        f['ObservationStart'] = start
        stop = core.G3Time(sys.argv[2])
        f['ObservationStop'] = start
        self.obsid = int((start.time - core.G3Time('20170101_000000').time)/core.G3Units.s)
        f['ObservationID'] = self.obsid
        out = [f] + self.buffer + [frame]
        self.buffer = []
        self.obsframe_out = True
        return out
pipe.Add(fabricate_obsframe)

# Apply online pointing
pipe.Add(op.CorrectBoresightPointing, model='OnlinePointingModel',
         flags=['az_tilts', 'el_tilts', 'flexure','collimation',
                'refraction']) #, 'thermolin'])

# Make calibrator data look like SPT3G
pipe.Add(sptpol.DiscretizeCalibratorTimestream)

# First-pass pointing
pipe.Add(coordinateutils.azel.LocalToAstronomicalPointing,
         az_timestream='OnlineBoresightAz', el_timestream='OnlineBoresightEl',
         ra_timestream='OnlineBoresightRa', dec_timestream='OnlineBoresightDec')
pipe.Add(FillCoordTransRotations,
         transform_store_key = 'OnlineRaDecRotation',
         bs_az_key = 'RawBoresightAz', bs_el_key = 'RawBoresightEl',
         bs_ra_key = 'OnlineBoresightRa', bs_dec_key = 'OnlineBoresightDec',
         do_bad_transform = True)

pipe.Add(core.Dump)
def fcons(frame, seqno):
    outdir = os.path.join(sys.argv[3], frame['SourceName'], str(frame['ObservationID']))
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    return os.path.join(outdir, '%04d.g3' % seqno)
pipe.Add(core.G3MultiFileWriter, filename=fcons, size_limit=1024*1024*1024)

pipe.Run()
