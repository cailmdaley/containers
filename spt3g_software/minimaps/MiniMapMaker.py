from spt3g import core, coordinateutils
import numpy

@core.indexmod
class MiniMapMaker(object):
	'''
	Jason Gallichio-style "mini-map" maker used to look at sky arcs.
	Generates maps for each bolometer by binning every scan in azimuth
	and then placing subsequent scans in vertical series. This is useful
	for generating maps of large-scale patterns from a single observation
	without worrying about holes in the map from the scan pattern.
	'''
	def __init__(self, el_rows, az_bins, pointing='BoresightAz', bolots='RawTimestreams_I', detectors=[]):
		self.bolots = bolots
		self.detectors = detectors
		self.pointing = pointing
		self.az_bins = az_bins
		self.map_template = coordinateutils.FlatSkyMap(
			len(az_bins) - 1,
			el_rows,
			res=1,
			x_res=numpy.abs(numpy.diff(az_bins)[0]),
			coord_ref=coordinateutils.MapCoordReference.Local,
		)
		self.scan_row = 0
		self.maps = {}
		self.field = None
	def __call__(self, frame):
		if frame.type == core.G3FrameType.EndProcessing:
			map_frame = core.G3Frame(core.G3FrameType.Map)
			map_frame['SourceName'] = self.field
			for m in self.maps.keys():
				map_frame[m] = self.maps[m]
			return [map_frame, frame]

		if frame.type != core.G3FrameType.Scan or 'Turnaround' in frame:
			return [frame]

		out = []
		if self.field is not None and self.field != frame['SourceName']:
			map_frame = core.G3Frame(core.G3FrameType.Map)
			map_frame['SourceName'] = self.field
			for m in self.maps.keys():
				map_frame[m] = self.maps[m]
			self.maps = {}
			self.scan_row = 0
			out.append(map_frame)
		if self.field != frame['SourceName']:
			self.field = frame['SourceName']

		az_pointing = frame[self.pointing]
		bolos = self.detectors
		if len(bolos) == 0:
			bolos = list(frame[self.bolots].keys())
		for bolo in bolos:
			if bolo not in self.maps:
				self.maps[bolo] = coordinateutils.FlatSkyMap(self.map_template)
			m = self.maps[bolo]
			stripe = numpy.histogram(az_pointing, weights=frame[self.bolots][bolo], bins=self.az_bins)[0]/numpy.histogram(az_pointing, bins=self.az_bins)[0]
			numpy.asarray(m)[self.scan_row] = stripe
		self.scan_row += 1
		out.append(frame)

		return out

@core.indexmod
class ObsThumbMapper(object):
	'''
	Slow python mapmaker that works like the MiniMap maker but force aligns
	all scans to approximate RA rather than Az. Also bins in El rather than
	appending scan lines.
	'''
	def __init__(self, el_bins, az_bins, az='BoresightAz', el='BoresightEl', bolots='RawTimestreams_I'):
		self.bolots = bolots
		self.az = az
		self.el = el
		self.el_bins = el_bins
		self.az_bins = az_bins
		self.map_template = coordinateutils.FlatSkyMap(
			len(az_bins) - 1,
			len(el_bins) - 1,
			res=numpy.abs(numpy.diff(el_bins)[0]),
			x_res=numpy.abs(numpy.diff(az_bins)[0]),
			coord_ref=coordinateutils.MapCoordReference.Local,
			proj=coordinateutils.MapProjection.ProjCAR,
		)
		self.first_el = None
		self.maps = {}
		self.field = None
	def __call__(self, frame):
		if frame.type == core.G3FrameType.EndProcessing:
			map_frame = core.G3Frame(core.G3FrameType.Map)
			for m in self.maps.keys():
				mp = coordinateutils.FlatSkyMap(self.map_template)
				numpy.asarray(mp)[:] = self.maps[m][0]/self.maps[m][1]
				map_frame[m] = mp
			return [map_frame, frame]

		if frame.type != core.G3FrameType.Scan or 'Turnaround' in frame:
			return [frame]

		out = []
		if self.field is not None and self.field != frame['SourceName']:
			map_frame = core.G3Frame(core.G3FrameType.Map)
			map_frame['SourceName'] = self.field
			for m in self.maps.keys():
				mp = coordinateutils.FlatSkyMap(self.map_template)
				numpy.asarray(mp)[:] = self.maps[m][0]/self.maps[m][1]
				map_frame[m] = mp
			self.maps = {}
			out.append(map_frame)
			self.field = frame['SourceName']
			self.first_el = 0

		az = frame[self.az]
		el = frame[self.el]
		if self.first_el is None:
			self.first_el = numpy.min(el)
		for bolo in frame[self.bolots].iteritems():
			if bolo[0] not in self.maps:
				self.maps[bolo[0]] = (numpy.zeros((len(self.az_bins)-1, len(self.el_bins)-1)), numpy.zeros((len(self.az_bins)-1, len(self.el_bins)-1)))
			m = self.maps[bolo[0]]
			m[0][:] += numpy.histogram2d(az - numpy.min(az), el - self.first_el, (self.az_bins, self.el_bins), weights=bolo[1])[0]
			m[1][:] += numpy.histogram2d(az - numpy.min(az), el - self.first_el, (self.az_bins, self.el_bins))[0]
		out.append(frame)

		return out


