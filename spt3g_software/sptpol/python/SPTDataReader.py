from sptpol_software.data.readout import SPTDataReader as SPTPolDataReader
from sptpol_software.util.files import readSptData
from sptpol_software.util.time import SptDatetime
from spt3g import core, calibration
import sptpolflaggingutils
import numpy

class SPTDataReaderCore(object):
	def __init__(self, ignore_ts_flags = [], invert_ts_flags = [],
		     ignore_bolo_flags = ['has_time_const', 'good_angle_fit', 
					  'good_xpol_fit'],
		     invert_bolo_flags = ['has_pointing', 'has_polcal'],
                     enforce_partner_good = True ):
		self.ignore_ts_flags_ = ignore_ts_flags
		self.ignore_bolo_flags_ = ignore_bolo_flags
		self.invert_ts_flags_ = invert_ts_flags
		self.invert_bolo_flags_ = invert_bolo_flags
		self.hwm_emitted = False
		self.enforce_partner_good_ = enforce_partner_good
	def __call__(self, frame):
		try:
			scan = next(self.scan_iter)
		except StopIteration:
			return []

		frame = core.G3Frame(core.G3FrameType.Scan)
		ts = core.G3TimestreamMap()
		for i in range(len(self.data.observation.bolo_id_ordered)):
			if self.units == core.G3TimestreamUnits.Counts:
				# Use 3G-type bolometer timestream values where the 24-bit samples are right-aligned instead of left-aligned
				t = core.G3Timestream(numpy.asarray(numpy.asarray(self.data.bolodata_array[i][scan.scan_slice], dtype='int32') >> 8, dtype='float64'))
			else:
				t = core.G3Timestream(self.data.bolodata_array[i][scan.scan_slice])
			t.units = self.units
			ts[self.data.observation.bolo_id_ordered[i]] = t
		try:
			ts.start = core.G3Time(scan.start_time.formatDate('uni') * core.G3Units.s)
			ts.stop = core.G3Time(scan.stop_time.formatDate('uni') * core.G3Units.s)
		except AttributeError:
			# Some IDFs are missing this for reasons
			ts.start = core.G3Time(SptDatetime(self.data.antenna.track_utc[scan.start_index]).formatDate('uni') * core.G3Units.s)
			ts.stop = core.G3Time(SptDatetime(self.data.antenna.track_utc[scan.stop_index]).formatDate('uni') * core.G3Units.s)
		frame['CalTimestreams'] = ts
		frame['BoresightAz'] = core.G3Timestream(self.data.antenna.track_actual[0][scan.scan_slice] * core.G3Units.deg)
		frame['BoresightEl'] = core.G3Timestream(self.data.antenna.track_actual[1][scan.scan_slice] * core.G3Units.deg )

		frame['BoresightRa'] = core.G3Timestream(self.data.antenna.ra[scan.scan_slice] * core.G3Units.deg)
		frame['BoresightDec'] = core.G3Timestream(self.data.antenna.dec[scan.scan_slice] * core.G3Units.deg)

		frame['TimestreamWeights'] = core.G3MapDouble()

		bolo_ids = self.data.observation.bolo_id_ordered
		bolo_partner_ids = {}
		if hasattr(self.data, 'detector_parameters'):
			partner_index = self.data.detector_parameters.index_of_pol_partner
			for i in range(len( bolo_ids )):
				bolo_partner_ids[bolo_ids[i]] = bolo_ids[partner_index[i]]
		if hasattr(self.data, 'timestream_weights'):
			tmp_tswgt = self.data.timestream_weights
			for i in range(len(bolo_ids)):
				frame['TimestreamWeights'][bolo_ids[i]] = tmp_tswgt[i]
                
		outqueue = []
		if self.hwm_emitted is False:
			self.hwm_emitted = True
			calframe = core.G3Frame(core.G3FrameType.Calibration)
			cal = calibration.BolometerPropertiesMap()
			nomcal = calibration.BolometerPropertiesMap()
			
			calframe['TimeConst'] = core.G3MapDouble()
                
			for i, bid in enumerate(bolo_ids):
				bcal = calibration.BolometerProperties()
				bcal.physical_name = bid
				nombcal = calibration.BolometerProperties()
				nombcal.physical_name = bid

				bprops = self.data.rec[bid]
				bcal.x_offset = bprops.pointing_offset[1] * core.G3Units.arcmin
				nombcal.x_offset = bcal.x_offset
				bcal.y_offset = bprops.pointing_offset[0] * core.G3Units.arcmin
				nombcal.y_offset = bcal.y_offset
				bcal.pol_angle = bprops.pol_angle_deg * core.G3Units.deg
				nombcal.pol_angle = bprops.pol_angle_nominal_deg * core.G3Units.deg
				bcal.pol_efficiency = bprops.pol_eff
				nombcal.pol_efficiency = bprops.pol_eff_nominal
				bcal.band = bprops.band * core.G3Units.GHz
				nombcal.band = bcal.band

				cal[bid] = bcal
				nomcal[bid] = nombcal
				calframe['TimeConst'][bid] = bprops.time_const * core.G3Units.us

			calframe['BolometerProperties'] = cal
			calframe['NominalBolometerProperties'] = nomcal
			outqueue.append(calframe)
				
                #frame['ScanIsBad'] = scan['flags']
                
		frame['ScanNumber']  = -1
		frame['ObservationName'] = "SHAMESHAMESHAMESHAME"
                
		if hasattr(self.data, 'detector_parameters'):
			#import pdb, rlcompleter; pdb.Pdb.complete = rlcompleter.Completer(locals()).complete; pdb.set_trace()
			

			ts_flag_names = dict(self.data.timestream_flag_names)
			bolo_flag_names = dict(self.data.bolometer_flag_names)

			ts_flags = self.data.timestream_flags
			bolo_flags= self.data.bolometer_flags

			constant_flags = sptpolflaggingutils.get_constant_flagging(ts_flag_names, ts_flags, 
										   bolo_flag_names, bolo_flags,
										   bolo_ids, bolo_partner_ids,
										   ignore_ts_flags = self.ignore_ts_flags_, 
										   invert_ts_flags = self.invert_ts_flags_,
										   ignore_bolo_flags = self.ignore_bolo_flags_,
										   invert_bolo_flags = self.invert_bolo_flags_,
										   enforce_partner_good = self.enforce_partner_good_)
			scan_flags = sptpolflaggingutils.get_scan_flags(scan.is_bad_channel, bolo_ids, 
									enforce_partner_good = self.enforce_partner_good_, 
									bolo_partner_ids = bolo_partner_ids)

			scan_flags.update(constant_flags)



			frame['Flags'] = sptpolflaggingutils.convert_flag_dic_to_g3map(scan_flags)
		outqueue.append(frame)

		return outqueue
                
@core.indexmod
class SPTDataReader(SPTDataReaderCore):
	def __init__(self, start_date, stop_date, **kwargs):
		'''See SPTDataReader from sptpol_software for keyword arguments. Start and stop dates are retained for readData()'''
		SPTDataReaderCore.__init__(self)
		self.data = SPTPolDataReader(experiment='SPTpol', start_date=start_date, stop_date=stop_date, **kwargs)
		self.start_date = start_date
		self.stop_date = stop_date
	def readData(self, **kwargs):
		'''See SPTDataReader.readData for keyword arguments. Start and stop dates are implicit from construction'''
		self.data.readData(self.start_date, self.stop_date, **kwargs)
		self.scan_iter = iter(self.data.scan)
		if self.data.bolodata.values()[0].units == 'Amperes':
			self.units = core.G3TimestreamUnits.Current
		elif self.data.bolodata.values()[0].units == 'K_cmb' or self.data.bolodata.values()[0].units == 'K_CMB':
			self.units = core.G3TimestreamUnits.Tcmb
		elif self.data.bolodata.values()[0].units == 'Watts':
			self.units = core.G3TimestreamUnits.Power
		else:
			self.units = core.G3TimestreamUnits.Counts


@core.indexmod
class IDFReader(SPTDataReaderCore):
	def __init__(self, filename):
		SPTDataReaderCore.__init__(self)
		self.data = readSptData(filename)
		self.scan_iter = iter(self.data.scan)
		if self.data.bolodata.values()[0].units == 'Amperes':
			self.units = core.G3TimestreamUnits.Current
		elif self.data.bolodata.values()[0].units == 'K_cmb' or self.data.bolodata.values()[0].units == 'K_CMB':
			self.units = core.G3TimestreamUnits.Tcmb
		elif self.data.bolodata.values()[0].units == 'Watts':
			self.units = core.G3TimestreamUnits.Power
		else:
			self.units = core.G3TimestreamUnits.Counts
		
