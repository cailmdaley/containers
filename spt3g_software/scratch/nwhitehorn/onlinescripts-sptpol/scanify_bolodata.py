#!/usr/bin/env python

from spt3g import core, dfmux, gcp, std_processing
import socket, argparse, os

parser = argparse.ArgumentParser(description='Merge GCP data with raw bolometer data, convert to scans, and compress')
parser.add_argument('output', metavar='/path/to/output', help='Base for output file names (/path/to/output). Files will be placed at /path/to/output/field_name/observation_id/001.g3, etc.')
parser.add_argument('gcparchives', metavar='gcparchives', help='Path to directory with GCP archive files')
parser.add_argument('input', metavar='input', nargs='+', help='Input files to read')

parser.add_argument('-v', dest='verbose', action='store_true', help='Verbose mode (print all frames)')
parser.add_argument('--max_file_size', default=1024, dest='max_file_size', help='Maximum file size in MB (default 1024)')
args = parser.parse_args()

pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename=args.input)
pipe.Add(std_processing.ARCInterposer, basedir=args.gcparchives)
pipe.Add(gcp.ARCExtract)
pipe.Add(std_processing.BuildScanFramesFromRawData, flac_compress=True)

# Need to put files somewhere based on field name. This is only in scan frames
# and the first frame isn't a scan frame. Place Observation frames in the
# appropriate places (at the beginning) whenever the observation changes.
class CollectObservationIds(object):
	def __init__(self):
		self.buffer = []
		self.curfield = None
		self.obsid = None
	def __call__(self, frame):
		if not 'SourceName' in frame:
			if self.curfield is not None:
				frame['SourceName'] = self.curfield
				frame['ObservationID'] = self.obsid
				return frame
			else:
				self.buffer.append(frame)
				return []
		else:
			out = []
			if 'ObservationID' not in frame:
				frame['ObservationID'] = 'NA'
			if frame['SourceName'] != self.curfield or frame['ObservationID'] != self.obsid:
				f = core.G3Frame(core.G3FrameType.Observation)
				f['SourceName'] = frame['SourceName']
				f['ObservationID'] = frame['ObservationID']
				out.append(f)
			self.curfield = frame['SourceName']
			self.obsid = frame['ObservationID']
			for f in self.buffer:
				f['SourceName'] = self.curfield
				f['ObservationID'] = self.obsid
			out += self.buffer
			self.buffer = []
			out.append(frame)
			return out
pipe.Add(CollectObservationIds)

if args.verbose:
	pipe.Add(core.Dump)

class filename_constructor(object):
	def __init__(self):
		self.seqno = 1
	def __call__(self, frame, seqno):
		dir = os.path.join(args.output, frame['SourceName'])
		dir = os.path.join(dir, frame['ObservationID'])
		if not os.path.exists(dir):
			core.log_notice('Making new data directory %s' % dir, unit='DataOrganization')
			os.makedirs(dir)
			self.seqno = 1
		path = os.path.join(dir, '%04d.g3' % self.seqno)
		while os.path.exists(path):
			self.seqno += 1
			path = os.path.join(dir, '%04d.g3' % self.seqno)
		core.log_info('Starting new output file %s' % path, unit='DataOrganization')
		self.seqno += 1
		return path

pipe.Add(core.G3MultiFileWriter, filename=filename_constructor(), size_limit=args.max_file_size*1024*1024, divide_on=[core.G3FrameType.Observation])

pipe.Run()

