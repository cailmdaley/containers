#!/usr/bin/env python

from spt3g import core, gcp, std_processing
from spt3g.std_processing.status_db import ScanifyDatabase
import os
import argparse
import pandas

# This script reads GCP ARC files, looking for the start and stop times
# of observations that begin with a given time range.

# The end of processing is always signalled by halt_processing(), issued
# on the first GCP frame following the end of the last observation starting
# in the provided time window. Absent such a call, the end of input data
# (for example, if the ARC file reader runs out of ARC files) will cause
# an exception to be thrown. As such, the output of the script is guaranteed
# correct and complete if (and only if) it runs to completion and returns 0.

core.set_log_level(core.G3LogLevel.LOG_INFO, 'ARCTimerangeReader')

parser = argparse.ArgumentParser(description='Add observations from ARC files to database')
parser.add_argument('gcparchives', metavar='gcparchives',
  help='Path to directory with GCP archive files')

parser.add_argument('-v', '--verbose', action='store_true',
  help='Verbose mode (print all frames)')
parser.add_argument('--first-obs', default=None,
  help="If set, will process only observations starting after the time specified. "
  "If a path to a file, read the timestamp from the file. "
  "If 'latest', and --database is given, set to the stop time of the most recent observation. ")
parser.add_argument('--last-obs', default=None,
  help='If set, will process only observations starting before the time specified. '
  'Will process data after this time if it is part of a continuing observation.')
parser.add_argument('--database', default = None,
  help='Processed observations are registered in this file, so that the relevant '
  'autoprocessing can be run later.')
parser.add_argument('--overwrite', default=False, action='store_true',
  help='Overwrite observations that already exist in the database')
parser.add_argument('--check', default=False, action='store_true',
                    help='Update timestamps and set to rerun any observations found')
args = parser.parse_args()

# Handle start/stop times
first_obs = None
last_obs = None
first_obs_file = None
if args.first_obs is not None:
    if os.path.exists(str(args.first_obs)):
        first_obs_file = args.first_obs
        first_obs = core.G3Time(open(args.first_obs, 'r').read().strip())
    elif args.first_obs == 'latest':
        if args.database is not None:
            db = ScanifyDatabase(filename=args.database)
            first_obs = core.G3Time(str(db.get_latest_stop())) + core.G3Units.s
            db.close()
        else:
            raise ValueError("Database required if first_obs=latest")
    else:
        first_obs = core.G3Time(args.first_obs)
if args.last_obs is not None:
    last_obs = core.G3Time(args.last_obs)

def write_first_obs(timestamp):
    if first_obs_file is None:
        return
    with open(first_obs_file, 'w') as f:
        ts = str(timestamp - 10 * core.G3Units.s)
        print('Recording {}: {}'.format(first_obs_file, ts))
        f.write(ts + '\n')

# Begin processing
pipe = core.G3Pipeline()

# Decode GCP data
print('Search range: {} to {}'.format(first_obs, last_obs))
pipe.Add(std_processing.ARCTimerangeReader, basedir=args.gcparchives,
         start_time=first_obs, stop_time=None, no_more_data_error=False)
pipe.Add(gcp.ARCExtract)

# Keep track of the most recent GcpSlow frame that has been read in
# Log its timestamp as a key in the EndProcessing frame,
# to be handled in ObsStopTimes
class TrackLastSlow(object):
    def __init__(self):
        self.last_slow = None
    def __call__(self, frame):
        if frame.type == core.G3FrameType.GcpSlow:
            self.last_slow = frame
        elif frame.type == core.G3FrameType.EndProcessing:
            frame['LastSample'] = self.last_slow['array']['frame']['utc']
        return frame
pipe.Add(TrackLastSlow)

# Drop data samples when GCP is just staring into space
class DropOrphanSamples(object):
    def __init__(self):
        self.in_obs = False
    def __call__(self, frame):
        if frame.type == core.G3FrameType.GcpSlow:
            self.in_obs = (frame['SourceName'] != '' and frame['SourceName'] != '(none)'
                           and pandas.notnull(frame['SourceName'])
                           and 'analyze' in frame['GCPFeatureBits'])
            if last_obs is not None and frame['array']['frame']['utc'] >= last_obs \
               and not self.in_obs:
                write_first_obs(frame['array']['frame']['utc'])
                core.G3Pipeline.halt_processing()
            return self.in_obs
        if frame.type == core.G3FrameType.EndProcessing:
            if not self.in_obs:
                write_first_obs(frame['LastSample'])
            return frame
pipe.Add(DropOrphanSamples)

# Drop scans while slewing to a source, which get feature bits of ['analyze']
# and a source name of 'current'. Since we know 'analyze' is in the list from
# above, just make sure either there are more flags or the source name isn't
# 'current'
pipe.Add(lambda fr: fr.type != core.G3FrameType.GcpSlow or
  len(fr['GCPFeatureBits']) > 1 or fr['SourceName'] != 'current')

# Drop data from real observations (all data now) that started before our start
# time but have continued after our end time (analog to DropOrphanSamples, but
# for observations we just don't care about rather than just non-observing
# periods).
def DropThroughGoingData(frame):
    if frame.type != core.G3FrameType.GcpSlow:
        return

    # The observation scope starts at a time recorded by GCP in the obs_id
    # register, which we know to be related to the number of seconds since
    # the 3G epoch (Jan. 1, 2017). This is used to know if the first slow
    # frame seen is actually from the beginning of the observation or just
    # the first one we got. This should be the *only* place that knows about
    # this definition.
    id_equivalent_start_time = core.G3Time('20170101_000000') + \
                               frame['ObservationID']*core.G3Units.s

    if first_obs is None or id_equivalent_start_time >= first_obs:
        return True

    # We're in an observation from before our window now. Before dropping
    # the data, see if we also want to halt processing (no more useful data
    # coming).
    if last_obs is not None and frame['array']['frame']['utc'] >= last_obs:
        print('Stopping on ' + str(frame))
        core.G3Pipeline.halt_processing()

    # And maybe we transitioned smoothly to past-relevant observations...
    if last_obs is not None and id_equivalent_start_time >= last_obs:
        print('Stopping on ' + str(frame))
        core.G3Pipeline.halt_processing()

    return False
pipe.Add(DropThroughGoingData)

# Extract observation parameters from GCP frames and prefix them with an
# Observation frame whenever the observation changes. 
class CollectObservationIds(object):
    def __init__(self):
        self.curfield = None
        self.obsid = None
    def __call__(self, frame):
        if 'SourceName' not in frame:
            return

        if frame['SourceName'] != self.curfield or \
           frame['ObservationID'] != self.obsid:
            fr = core.G3Frame(core.G3FrameType.Observation)
            fr['SourceName'] = frame['SourceName']
            fr['ObservationID'] = frame['ObservationID']
            fr['ObservationStart'] = frame['array']['frame']['utc'] - \
                                    core.G3Units.s
            print('Found observation %s/%d' % (fr['SourceName'], fr['ObservationID']))
            self.curfield = frame['SourceName']
            self.obsid = frame['ObservationID']
            return [fr, frame]

pipe.Add(CollectObservationIds)

# Grab observation stop times. Implicitly drops slow frames.
class ObsStopTimes(object):
    def __init__(self):
        self.obs_frame = None
        self.last_slow = None
    def __call__(self, fr):
        # If processing ends before an observation is complete
        # (e.g. if we've reached the most recent arcfile),
        # then return without recording the observation stop
        if fr.type == core.G3FrameType.EndProcessing:
            if self.last_slow is not None and \
               fr['LastSample'] == self.last_slow['array']['frame']['utc']:
                if self.obs_frame is not None:
                    write_first_obs(self.obs_frame['ObservationStart'])
                return fr
        out = []
        if fr.type == core.G3FrameType.Observation or \
           fr.type == core.G3FrameType.EndProcessing:
            if self.obs_frame is not None:
                self.obs_frame['ObservationStop'] = \
                    self.last_slow['array']['frame']['utc']
                out.append(self.obs_frame)
            if fr.type == core.G3FrameType.Observation:
                self.obs_frame = fr
            else:
                out.append(fr)
        if fr.type == core.G3FrameType.GcpSlow:
            self.last_slow = fr
        return out
pipe.Add(ObsStopTimes)

# Cut to only Observation headers
pipe.Add(lambda fr: fr.type == core.G3FrameType.Observation)

# Cut 1-second observations that arise from a race condition in
# the way that GCP sets the analyze bit when an observation begins
def CheckObs(frame):
    if frame.type != core.G3FrameType.Observation:
        return
    if frame['ObservationStop'].time - frame['ObservationStart'].time <= core.G3Units.s:
        print("Dropping orphan observation %s/%d" % (frame['SourceName'], frame['ObservationID']))
        return False
pipe.Add(CheckObs)

if args.verbose:
    pipe.Add(core.Dump)

if args.database is not None:
    # Filename constructor has cached a list of observations
    class RecordObservations(object):
        def __init__(self):
            self.db = None

        def __del__(self):
            if self.db is not None:
                self.db.close()
                self.db = None

        def __call__(self, fr):
            if fr.type != core.G3FrameType.Observation:
                return

            # Only open the database if we find an observation
            if self.db is None:
                self.db = ScanifyDatabase(filename=args.database)

            source = fr['SourceName']
            obsid = fr['ObservationID']
            obs_start = pandas.Timestamp(fr['ObservationStart'].isoformat(), tz='UTC')
            obs_stop = pandas.Timestamp(fr['ObservationStop'].isoformat(), tz='UTC')
            rerun = False

            idx = self.db.match(source, obsid, return_index=True)
            if len(idx):
                if not args.overwrite and not args.check:
                    print('Observation %s/%d already registered in scanify database' % (source, obsid))
                    return

                if args.check:
                    # Check that observation start and stop times haven't changed
                    entry = self.db.data.loc[idx[0]]
                    if obs_start == entry['obs_start'] and obs_stop == entry['obs_stop']:
                        print('Observation %s/%d is up to date' % (source, obsid))
                        return
                    elif entry['status_scanify'] == 'success':
                        print('Observation %s/%d is already scanified, not updating.' % (source, obsid))
                        return
                    else:
                        print('Updating observation %s/%d' % (source, obsid))
                        rerun = True
            else:
                print('Registering observation %s/%d in scanify database' % (source, obsid))

            self.db.update(source, obsid, status_fullrate=None,
                           status_downsampled=None, commit=True,
                           obs_start=obs_start, obs_stop=obs_stop,
                           status_scanify='rerun' if rerun else None)

    record_observations = RecordObservations()
    pipe.Add(record_observations)

pipe.Run()

