import numpy, sys
from spt3g import core, gcp, std_processing, dfmux, coordinateutils, sptpol
from spt3g.coordinateutils.coordsysmodules import FillCoordTransRotations
from spt3g.pointing.offset_estimation import CalcAzElToRaDecOffsetFromBoresight
import spt3g.pointing.offline_pointing as op

# Usage: convert-sptpol-data.py start-time stop-time output/%04u.g3
# Converts SPTpol data read from ARC files into data that is supposed to look
# like SPT3G data fresh off the satellite.

#core.set_log_level(core.G3LogLevel.LOG_DEBUG)

# In case you get warnings from astropy about timekeeping:
from astropy.utils import iers
#iers.IERS.iers_table = iers.IERS_A.open(iers.IERS_A_URL)

pipe = core.G3Pipeline()
pipe.Add(std_processing.ARCTimerangeReader, start_time=sys.argv[1], stop_time=sys.argv[2])
#pipe.Add(core.Dump)
pipe.Add(gcp.ARCExtract)
pipe.Add(gcp.UnpackSPTpolHKData)
pipe.Add(gcp.GCPMuxDataDecoder)
pipe.Add(std_processing.BuildScanFramesFromRawData, flac_compress=True)

# Add Observation frames
obsframe_out = False
def fabricate_obsframe(frame):
    global obsframe_out
    if obsframe_out:
        return
    if 'SourceName' not in frame:
        return
    f = core.G3Frame(core.G3FrameType.Observation)
    f['SourceName'] = frame['SourceName']
    start = core.G3Time(sys.argv[1])
    f['ObservationStart'] = start
    f['ObservationID'] = int((start.time - core.G3Time('20170101_000000').time)/core.G3Units.s)
    return [f, frame]
pipe.Add(fabricate_obsframe)

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

# Apply online pointing
pipe.Add(op.CorrectBoresightPointing, online=True,
         flags=['az_tilts', 'el_tilts', 'flexure','collimation',
                'refraction', 'thermolin'])

# Make calibrator data look like SPT3G
pipe.Add(sptpol.DiscretizeCalibratorTimestream)

# First-pass pointing
pipe.Add(coordinateutils.azel.LocalToAstronomicalPointing,
         az_timestream='RawBoresightAz', el_timestream='RawBoresightEl',
         ra_timestream='BoresightRa', dec_timestream='BoresightDec')

pipe.Add(FillCoordTransRotations,
         transform_store_key = 'OnlineRaDecRotation',
         bs_az_key = 'RawBoresightAz', bs_el_key = 'RawBoresightEl',
         bs_ra_key = 'OnlineBoresightRa', bs_dec_key = 'OnlineBoresightDec',
         do_bad_transform = True)

pipe.Add(coordinateutils.azel.LocalToAstronomicalPointing,
         az_timestream='OnlineBoresightAz', el_timestream='OnlineBoresightEl',
         ra_timestream='OnlineBoresightRa', dec_timestream='OnlineBoresightDec')

pipe.Add(CalcAzElToRaDecOffsetFromBoresight,
         correct_bs_pointing_flags = ['az_tilts', 'el_tilts', 'flexure','collimation',
                                      'refraction', 'thermolin'],
         bs_az_key = 'RawBoresightAz', bs_el_key='RawBoresightEl',
         offset_az_key='OffsetBoresightAz', offset_el_key='OffsetBoresightEl', 
         offset_ra_key='OnlineOffsetRa', offset_dec_key='OnlineOffsetDec')

pipe.Add(FillCoordTransRotations,
         transform_store_key = 'OfflineRaDecRotation',
         bs_az_key = 'RawBoresightAz', bs_el_key='RawBoresightEl',
         bs_ra_key = 'OfflineBoresightRa', bs_dec_key='OfflineBoresightDec',
         offset_az_key='OffsetBoresightAz', offset_el_key='OffsetBoresightEl', 
         offset_ra_key='OfflineOffsetRa', offset_dec_key='OfflineOffsetDec')         

pipe.Add(core.Dump)
pipe.Add(core.G3MultiFileWriter, filename=sys.argv[3], size_limit=1024*1024*1024)

pipe.Run()
