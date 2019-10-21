#!/usr/bin/env python

from spt3g import core, dfmux, gcp, std_processing, coordinateutils
from spt3g.pointing import CalculateCoordTransRotations
import argparse, os, shutil, glob
import pandas
import subprocess as sp
import socket

core.set_log_level(core.G3LogLevel.LOG_INFO, 'G3Reader')

parser = argparse.ArgumentParser(description=
  'Merge GCP data with raw bolometer data for one observation, convert to scans, and compress')
parser.add_argument('output', metavar='/path/to/output',
                    help='Base for output file names (/path/to/output). Files will be placed at /path/to/output/field_name/observation_id/0000.g3, etc.')
parser.add_argument('gcparchives', metavar='gcparchives',
                    help='Path to directory with GCP archive files')
parser.add_argument('input', metavar='/path/to/input',
                    help='Path to directory with raw timepoint G3 files')

parser.add_argument('-v', '--verbose', action='store_true',
                    help='Verbose mode (print all frames)')
parser.add_argument('--max-file-size', default=1024,
                    help='Maximum file size in MB (default 1024)')
parser.add_argument('--start-time', help='Observation start time.', required=True)
parser.add_argument('--stop-time', help='Observation stop time.', required=True)
parser.add_argument('--user-obs', default=False, action='store_true',
                    help="Create the observation even if GCP claims one does not exist. "
                    "The source name will be set to 'userdata' and the observation ID "
                    "will be fixed to the given start time")
parser.add_argument('--scratch', default=False, action='store_true',
                    help="Copy input files to the condor scratch directory, and copy "
                    "completed jobs back to the output directory when complete.")
args = parser.parse_args()

# parse time range
start_t = core.G3Time(args.start_time)
start_stamp = os.path.join(args.input, start_t.GetFileFormatString() + '.g3')
stop_t = core.G3Time(args.stop_time)
stop_stamp = os.path.join(args.input, stop_t.GetFileFormatString() + '.g3')

# select only files that contain data within the start/stop times
files = sorted(glob.glob(os.path.join(args.input, '*.g3')))
if len(files) and files[-1] < stop_stamp:
    raise RuntimeError(
        "Rsync from daqcontrol is not up to date.  Please wait to rerun this "
        "observation until the rsync process has caught up to {}.".format(str(stop_t))
    )
while len(files) and files[-1] >= stop_stamp:
    files = files[:-1]
while len(files) >= 2 and files[1] <= start_stamp:
    files = files[1:]

if not len(files):
    raise RuntimeError("No G3 Timepoint files found for this observation!")

scratch = None
if args.scratch:
    scratch = os.getenv('_CONDOR_SCRATCH_DIR')

if scratch:
    remote_input = args.input
    remote_output = args.output

    # create temporary directory structure on local disks
    args.input = os.path.join(scratch, 'raw')
    try:
        os.makedirs(args.input)
    except OSError:
        pass
    args.output = os.path.join(scratch, 'output')
    try:
        os.makedirs(args.output)
    except OSError:
        pass

    # copy input files over
    sp.check_call(
        'rsync -aviP {} {}'.format(' '.join(files), os.path.join(args.input, '.')),
        shell=True,
    )
    files = [os.path.join(args.input, os.path.basename(f)) for f in files]

# Begin processing
pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename=files)

# Don't even bother with data outside our window of interest
print(str(start_t), str(stop_t))
pipe.Add(lambda fr: fr.type != core.G3FrameType.Timepoint or
         (fr['EventHeader'] >= start_t and fr['EventHeader'] < stop_t))

# Fill missing timepoint frames with a frame full of NaNs
pipe.Add(dfmux.DataQualityTools.FillMissingTimepointFrames)

# Make sure that the raw data are not garbage
class FrameCounter(object):
    def __init__(self, threshold=0.5):
        self.threshold = threshold
        self.count_missing = 0
        self.count_total = 0

    def __call__(self, frame):
        if frame.type == core.G3FrameType.Timepoint:
            if frame.get('GarbageData', False):
                self.count_missing += 1
            self.count_total += 1
        elif frame.type == core.G3FrameType.EndProcessing:
            if not self.count_total:
                core.log_fatal("No timepoint data found!")
            m = self.count_missing / float(self.count_total)
            if m > self.threshold:
                core.log_fatal("Missing {:.1f}% of timepoint data!".format(m * 100))
        return
fc = FrameCounter()
pipe.Add(fc)

# Decode GCP data
pipe.Add(std_processing.ARCInterposer, basedir=args.gcparchives)
pipe.Add(gcp.ARCExtract)

# Skip GCP frames outside of our window of interest
pipe.Add(lambda fr: fr.type != core.G3FrameType.GcpSlow or
         (fr['array']['frame']['utc'] < stop_t))

# Drop pointless duplicate metadata
pipe.Add(core.DeduplicateMetadata)

# Build to scans
pipe.Add(std_processing.BuildScanFramesFromRawData, flac_compress=True, 
         max_scan_length=50000)

if args.user_obs:
    # Inject a fake source name and observation ID for the requested data
    def make_user_obs(frame):
        if frame.type != core.G3FrameType.Scan:
            return

        del frame['SourceName']
        frame['SourceName'] = 'userdata'

        del frame['ObservationID']
        frame['ObservationID'] = std_processing.utils.time_to_obsid(start_t)

        return frame

    pipe.Add(make_user_obs)

# Drop scans with no set feature bits. NB: should never happen! Assert?
pipe.Add(lambda fr: fr.type != core.G3FrameType.Scan or
  'analyze' in fr['GCPFeatureBits'] or args.user_obs)

# Drop scans with no set source name.  NB: should never happen! Assert?
pipe.Add(lambda fr: fr.type != core.G3FrameType.Scan or
  (fr['SourceName'] != '' and fr['SourceName'] != '(none)' and
   pandas.notnull(fr['SourceName'])))

# Drop 1-sample scans that can occur from the interaction of the 1-second
# granularity of GCP's analyze bit and the 10 ms offset between the time of
# the first tracker sample and the notional time of the frame. These will later
# have infinite sample rates and are guaranteed to be thrown away anyway
# since these scans occur only immediately before we start scanning.
pipe.Add(lambda fr: 'RawBoresightAz' not in fr or len(fr['RawBoresightAz']) > 1)

# Drop metadata that does not apply to any scans
pipe.Add(core.DropOrphanMetadata)

# Recreate the observation header. This needs some data from the scan
# frames, so deal with buffering a bit.
class BuildObservationFrame(object):
    def __init__(self):
        self.buffer = []
        self.curfield = None
        self.obsid = None
    def __call__(self, frame):
        if frame.type == core.G3FrameType.EndProcessing:
            return
        if self.curfield is None and self.obsid is None:
            if not 'SourceName' in frame:
                self.buffer.append(frame)
                return []
            else:
                f = core.G3Frame(core.G3FrameType.Observation)
                f['SourceName'] = frame['SourceName']
                f['ObservationID'] = frame['ObservationID']
                self.curfield = frame['SourceName']
                self.obsid = frame['ObservationID']
                f['ObservationStart'] = start_t
                f['ObservationStop'] = stop_t

                # check for IERS updates
                import astropy.time
                t = astropy.time.Time(stop_t.mjd, format='mjd')
                try:
                    t.ut1
                    iers_ok = True
                except:
                    iers_ok = False

                if not iers_ok:
                    from astropy.utils import iers
                    iers.conf.auto_max_age = None
                    msg = (
                        "Unable to update IERS table for scanifier job %s/%d, "
                        "Using extrapolated polar motion parameters instead. "
                        "Run 'check_iers' when the satellite is up to fix this."
                    ) % (f['SourceName'], f['ObservationID'])
                    core.log_warn(msg)
                    from spt3g.std_processing.transfer_tools import send_cron_email
                    send_cron_email(os.path.abspath(__file__), msg)

                # Grab the first commanded bench position as the
                # nominal position for this observation.  This is
                # safe, because the commanded bench position never
                # changes during an ob.
                bench_commanded_pos = core.G3MapDouble()
                bench_zeros = core.G3MapDouble()
                for key in frame['BenchCommandedPosition'].keys():
                    bench_commanded_pos[key] = frame['BenchCommandedPosition'][key][0]
                    bench_zeros[key] = frame['BenchZeros'][key][0]
                f['BenchCommandedPosition'] = bench_commanded_pos
                f['BenchZeros'] = bench_zeros

                out = [f]

                out += self.buffer
                self.buffer = []
                out.append(frame)
                return out
        else:
            return

pipe.Add(BuildObservationFrame)

# The nominal bench position is now in the Observation frame, so drop
# it from Scan frames.
pipe.Add(core.Delete, keys=['BenchCommandedPosition', 'BenchZeros'],
         type=core.G3FrameType.Scan)
# More metadata might be bogus now, so clean again
pipe.Add(core.DropOrphanMetadata)

# Add Online pointing to scans.
pipe.Add(
    CalculateCoordTransRotations,
    raw_az_key='RawBoresightAz',
    raw_el_key='RawBoresightEl',
    output='OnlineBoresight',
    transform_store_key='OnlineRaDecRotation',
    model='OnlinePointingModel',
    flags=['az_tilts', 'el_tilts', 'flexure', 'collimation', 'refraction']
    #, 'thermolin'] # Thermoline broken as of 4/13/17
)

if args.verbose:
    pipe.Add(core.Dump)

# Cache nominal calibration, if present, so that it can get teed out below.
# Also drop it so it doesn't end up in the normal data stream. This uses
# an obs ID-indexed dictionary of calibration frames to avoid assumptions
# about relative processing times of frames between two unrelated modules
# in the chain.
class calcache(object):
    def __init__(self):
        self.cache = []
        self.cal = {}
        self.lastcal = None
        self.curobs = None
    def __call__(self, fr):
        if fr.type == core.G3FrameType.EndProcessing:
            return

        if fr.type == core.G3FrameType.Calibration:
            self.cal[self.curobs] = fr
            self.lastcal = fr
            return []

        if fr.type == core.G3FrameType.Scan:
            # Data here, flush cache -- if we had no cal frame,
            # we aren't getting one now.
            out = self.cache + [fr]
            self.cache = []
            return out

        if fr.type == core.G3FrameType.Observation:
            self.curobs = (fr['SourceName'], fr['ObservationID'])
            self.cal[self.curobs] = self.lastcal

        # Cache observation + all other types so that the file namer
        # below will have calibration data when it needs them
        self.cache.append(fr)
        return []
obs_cals = calcache()
pipe.Add(obs_cals)

# Count scan frames to make sure data exist for this time range.
class ScanCounter(object):
    def __init__(self):
        self.count = 0

    def __call__(self, frame):
        if frame.type == core.G3FrameType.Scan:
            self.count += 1
        elif frame.type == core.G3FrameType.EndProcessing:
            if self.count == 0:
                core.log_fatal('No data found!')
        return
sc = ScanCounter()
pipe.Add(sc)

# The "filename_constructor" both formats file names and creates output
# directories when observations change. When it does the latter,
# tees out the nominal calibration data inherited from the DAQ if present
# to a file named nominal_online_cal.g3 in the current observation directory.
class filename_constructor(object):
    def __init__(self):
        self.seqno = 0
        self.seen = []
    def __call__(self, frame, seqno):
        obsdir = os.path.join(args.output, frame['SourceName'])
        obsdir = os.path.join(obsdir, str(frame['ObservationID']))
        if (frame['SourceName'], frame['ObservationID']) not in self.seen:
            core.log_notice('Making new data directory %s' % obsdir,
              unit='DataOrganization')
            if os.path.exists(obsdir):
                core.log_notice('Found old data directory %s, deleting' % obsdir,
                  unit='DataOrganization')
                shutil.rmtree(obsdir)
            os.makedirs(obsdir)
            self.seen.append((frame['SourceName'], frame['ObservationID']))
            self.seqno = 0

            # New directory, tee out cal data if present
            cur_obs = (frame['SourceName'], frame['ObservationID'])
            if cur_obs in obs_cals.cal and obs_cals.cal[cur_obs] is not None:
                w = core.G3Writer(os.path.join(obsdir, 'nominal_online_cal.g3'))
                w(obs_cals.cal[cur_obs])
                w(core.G3Frame(core.G3FrameType.EndProcessing))
                del obs_cals.cal[cur_obs]

        path = os.path.join(obsdir, '%04d.g3' % self.seqno)
        while os.path.exists(path):
            self.seqno += 1
            path = os.path.join(obsdir, '%04d.g3' % self.seqno)
        core.log_info('Starting new output file %s' % path,
          unit='DataOrganization')
        self.seqno += 1
        return path

# Finally, write data to disk
pipe.Add(core.G3MultiFileWriter, filename=filename_constructor(),
  size_limit=args.max_file_size*1024*1024,
  divide_on=[core.G3FrameType.Observation])

pipe.Run()

if scratch:
    # copy output files back to network storage
    try:
        sp.check_call(
            'rsync -aviP {} {}'.format(os.path.join(args.output, '*'), remote_output),
            shell=True, stderr=sp.STDOUT,
        )
    except sp.CalledProcessError as e:
        sp.check_call(
            'rsync -aviP {} {}'.format(os.path.join(args.output, '*'), '/scratch'),
            shell=True,
        )
        hostname = socket.gethostname()
        raise OSError(
            "Error transferring completed job to {}, stored temporarily in {}:/scratch.\n{}"
            .format(remote_output, hostname, e.output)
        )
