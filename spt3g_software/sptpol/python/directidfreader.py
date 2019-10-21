import os
from copy import copy
from functools import reduce
from . import sptpolflaggingutils
import numpy as np

from spt3g import core, calibration, dfmux
from spt3g.core import indexmod
def shitty_con_fn_to_observation_name(fn):
    bn = os.path.basename(fn)
    return bn[:bn.rfind('.')]




class LoadBoloPropsFromIDF(object):
    def __init__(self, fn):
        self.fn = fn
        self.is_sent = False

    def __call__(self, frame):
        if self.is_sent or frame.type != core.G3FrameType.Wiring:
            return
        self.is_sent = True
        import h5py
        idf_file = h5py.File(self.fn, 'r')
        wm = frame['WiringMap']
        id_mapper = frame['PhysicalBoloIDs']
        id_mapper_inv = {}
        for k in id_mapper.keys():
            id_mapper_inv[id_mapper[k]] = k

        band = float(idf_file['/band'].value) * core.G3Units.GHz

        bolo_ids = []
        for x in idf_file["/data/observation/bolo_id_ordered"][:]:
            if isinstance(x, str):
                bolo_ids.append(x)
            else:
                bolo_ids.append(x.decode())

        #/data/observation/bolo_watts_to_kcmb
        #/data/observation/bolo_adc_to_amps
        #/data/observation/bolo_adc_to_watts



        bprops = calibration.BolometerPropertiesMap()
        for i in range(len(bolo_ids)):
            bprop = calibration.BolometerProperties()
            bprop.x_offset = float(idf_file['/data/detector_parameters/pointing_offset_x'][i]*core.G3Units.arcmin)
            bprop.y_offset = float(idf_file['/data/detector_parameters/pointing_offset_y'][i]*core.G3Units.arcmin)
            bprop.pol_angle = float(idf_file['/data/detector_parameters/pol_angle'][i]*core.G3Units.rad)
            bprop.pol_efficiency = float(idf_file['/data/detector_parameters/pol_eff'][i])

            bprop.band = band
            bprop_key = str(id_mapper_inv[bolo_ids[i]])
            #core.log_warn('adding ', bprop_key)
            bprops[bprop_key] = bprop

        cframe = core.G3Frame(core.G3FrameType.Calibration)
        cframe['BolometerProperties'] = bprops

        return [frame, cframe]


class LoadPolAngsFromIDF(object):
    def __init__(self, fn):
        self.fn = fn

        self.cal_frame = None
        self.wiring_frame = None

    def __call__(self, frame):
        if frame.type == core.G3FrameType.Calibration:
            self.cal_frame = frame

        if frame.type == core.G3FrameType.Wiring:
            self.wiring_frame = frame
        if frame.type != core.G3FrameType.Wiring and frame.type != core.G3FrameType.Calibration:
            return
        elif self.cal_frame == None or self.wiring_frame == None:
            return []
        wm = self.wiring_frame['WiringMap']
        print(frame)
        print(self.wiring_frame)
        id_mapper = self.wiring_frame['PhysicalBoloIDs']
        id_mapper_inv = {}
        for k in id_mapper.keys():
            id_mapper_inv[id_mapper[k]] = str(k)

        import h5py
        idf_file = h5py.File(self.fn, 'r')

        bolo_ids = map(str, idf_file["/data/observation/bolo_id_ordered"][:])          
        bprops = copy(self.cal_frame['BolometerProperties'])

        

        for k in bprops.keys():
            bprops[k].pol_angle = 0
            bprops[k].pol_eff = 0
            bprops[k].band = 0
        for i, bid_phys in enumerate(bolo_ids):

            if bid_phys in id_mapper_inv:
                bid = id_mapper_inv[bid_phys]
            else:
                print(bid_phys)
                continue

            bprop = bprops[bid]
            bprop.pol_angle = float(idf_file['/data/detector_parameters/pol_angle'][i]*core.G3Units.rad)
            bprop.pol_efficiency = float(idf_file['/data/detector_parameters/pol_eff'][i])

            bprop.band = 1
            bprops[bid] = bprop


        del self.cal_frame['BolometerProperties']
        self.cal_frame['BolometerProperties'] = bprops
        return [self.wiring_frame, self.cal_frame]


def GenerateWiringMapFromSptpolHwm(hwm):
    import socket, struct
    wiring_map = dfmux.DfMuxWiringMap()
    for k in hwm.channels:
        chan_map = dfmux.DfMuxChannelMapping()
        bid = hwm(k,'bolo_id')
        if bid == None or bid == '':
            continue
        ip = hwm(k, 'dfmux_ip')
        int_ip = struct.unpack("!i", socket.inet_aton(ip))[0] 
        chan_map.board_ip = int_ip
        chan_map.board_serial = int(hwm(k, 'dfmux_id'))
        chan_map.board_slot = int(hwm(k, 'dfmux_id'))
        chan_map.crate_serial = 0
        chan_map.module = int(hwm(k, 'module')) - 1
        chan_map.channel = int(hwm(k, 'chan_num')) - 1
        wiring_map[bid] = chan_map
    return wiring_map

@indexmod
class DirectIdfReader(object):
    def __init__(self, filename, preload_data = False, 
                 include_turn_arounds = False, load_bolo_data = True,
                 ignore_ts_flags = [], invert_ts_flags = [],
                 ignore_bolo_flags = ['has_time_const', 'good_angle_fit', 
                                      'good_xpol_fit'],
                 invert_bolo_flags = ['has_pointing', 'has_polcal'],
                 enforce_partner_good = True,
                 number_of_scans_to_read = -1,

                 hwm_path = None
             ):
        from spt3g import core, todfilter, calibration, util

        import h5py
        assert(os.path.exists(filename))
        self.enforce_partner_good_ = enforce_partner_good
        self.ignore_ts_flags_ = ignore_ts_flags
        self.invert_ts_flags_ = invert_ts_flags
        self.ignore_bolo_flags_ = ignore_bolo_flags
        self.invert_bolo_flags_ = invert_bolo_flags
        
        self.sent_first_frame_ = False
        self.shitty_obs_name_ = shitty_con_fn_to_observation_name(filename)

        self.load_bolo_data_ = load_bolo_data


        if hwm_path != None:
            from pywtl.common.DFML.HWM.HardwareMap import HardwareMap
            hwm = HardwareMap(hwm_path)
            self.wiring_map_cache = GenerateWiringMapFromSptpolHwm(hwm)
        else:
            self.wiring_map_cache = None
 
        dims = core.IntVector()       
        rfffff = util.open_hdf5_read(filename)
        bitfield = util.read_hdf5_bitfield(rfffff, "/data/scan/is_bad_channel", dims)
        util.close_hdf5_file(rfffff)
        self.scan_ts_flag_information_ = np.array(bitfield).reshape(dims)


        self.idf_file_ = h5py.File(filename, 'r')
        self.band = float(self.idf_file_['/band'].value) * core.G3Units.GHz

        self.scan_starts_ = self.idf_file_["/data/scan/start_index"][:]
        self.scan_stops_ = self.idf_file_["/data/scan/stop_index"][:]
        self.scan_flags_ = self.idf_file_["/data/scan_flags"][:]

        self.bolo_ids_ = map(str, self.idf_file_["/data/observation/bolo_id_ordered"][:])          
        partner_index = self.idf_file_['/data/detector_parameters/index_of_pol_partner'][:] 
        self.bolo_partner_ids_ = {}
        for i in range(len( self.bolo_ids_ )):
            self.bolo_partner_ids_[self.bolo_ids_[i]] = self.bolo_ids_[partner_index[i]]

        self.current_scan_ = 0
        self.preload_data_ = preload_data
        if preload_data and load_bolo_data:
            self.bolo_data_storage_ = self.idf_file_['/data/bolodata_array'][:,:].astype('float64')
        self.number_of_scans_to_read_ = number_of_scans_to_read

        self.include_turn_arounds = include_turn_arounds
        self.is_turn_around_scan = False

    def __del__(self):
        if hasattr(self, "idf_file_"):
            self.idf_file_.close()

    def __call__(self, frame):
        if self.current_scan_ >= len(self.scan_starts_):
            return []

        from spt3g import core, todfilter, calibration

        frames = []
        #first send configuration frame if we need to
        if not self.sent_first_frame_:
            self.sent_first_frame_ = True
            oframe = core.G3Frame(core.G3FrameType.Observation)
            oframe['ObservationNumber'] = core.G3Int(int(''.join(self.shitty_obs_name_.split('_')[2:4])))
            oframe['SourceName'] = self.idf_file_['/data/header/source'].value
            frames.append(oframe)

            cframe = core.G3Frame(core.G3FrameType.Calibration)
            idf_to_frame_keys = [ ('/data/detector_parameters/pointing_offset_x', 'PointingXOffset'),
                                   ('/data/detector_parameters/pointing_offset_y', 'PointingYOffset'),
                                   ('/data/detector_parameters/pol_angle', 'PolarizationAngle'),
                                   ('/data/detector_parameters/pol_eff', 'PolarizationEfficiency'),
                                   ('/data/detector_parameters/time_const', 'TimeConst'),
                               ]
            amps_to_kcmb = core.G3MapDouble()

            bprops = calibration.BolometerPropertiesMap()
            taus = core.G3MapDouble()
            for i in range(len(self.bolo_ids_)):
                bprop = calibration.BolometerProperties()
                bprop.x_offset = float(self.idf_file_['/data/detector_parameters/pointing_offset_x'][i]*core.G3Units.arcmin)
                bprop.y_offset = float(self.idf_file_['/data/detector_parameters/pointing_offset_y'][i]*core.G3Units.arcmin)
                bprop.pol_angle = float(self.idf_file_['/data/detector_parameters/pol_angle'][i]*core.G3Units.rad)
                bprop.pol_efficiency = float(self.idf_file_['/data/detector_parameters/pol_eff'][i])

                bprop.wafer_id = self.bolo_ids_[i].split('.')[0]
                bprop.pixel_id = self.bolo_ids_[i][:-2]                
                bprop.physical_name = self.bolo_ids_[i]
                bprop.band = self.band

                tau = self.idf_file_['/data/detector_parameters/time_const'][i]*core.G3Units.ms
                taus[self.bolo_ids_[i]] = tau
                bprops[self.bolo_ids_[i]] = bprop

                alpha = self.idf_file_['/data/observation/bolo_watts_to_kcmb'][i]
                beta = self.idf_file_['/data/observation/bolo_adc_to_amps'][i]
                lam = self.idf_file_['/data/observation/bolo_adc_to_watts'][i]
                amps_to_kcmb[self.bolo_ids_[i]] = ( (alpha * lam / beta) )

            cframe['BolometerProperties'] = bprops
            cframe['TimeConst'] = taus
            cframe['id'] = 'mainconfig'

            cframe['AmpsToKcmb'] = amps_to_kcmb
    

            frames.append(cframe)

            if not self.wiring_map_cache is None:
                wframe = core.G3Frame(core.G3FrameType.Wiring)
                wframe['WiringMap'] = self.wiring_map_cache
                frames.append(wframe)

            ts_flag_names = dict(self.idf_file_['timestream_flag_names'])
            for k in ts_flag_names:
                ts_flag_names[k] = ts_flag_names[k][()]

            bolo_flag_names = dict(self.idf_file_['bolometer_flag_names'])
            for k in bolo_flag_names:
                bolo_flag_names[k] = bolo_flag_names[k][()]

            self.constant_flags = sptpolflaggingutils.get_constant_flagging(ts_flag_names, self.idf_file_['/data/timestream_flags'][:],
                                                        bolo_flag_names, self.idf_file_['/data/bolometer_flags'][:],
                                                        self.bolo_ids_, self.bolo_partner_ids_,
                                                        ignore_ts_flags = self.ignore_ts_flags_,
                                                        invert_ts_flags = self.invert_ts_flags_,
                                                        ignore_bolo_flags = self.ignore_bolo_flags_,
                                                        invert_bolo_flags = self.invert_bolo_flags_,
                                                        enforce_partner_good = self.enforce_partner_good_)

        #now fill the scan frame

        #makes certain to increment the current scan variable if there are flagged scans

        if (self.current_scan_ >= len(self.scan_flags_) or
            (self.current_scan_ > self.number_of_scans_to_read_  and self.number_of_scans_to_read_ > 0)):
            return []

        if (self.current_scan_ == len(self.scan_flags_)-1 and self.is_turn_around_scan):
            return []



        if self.include_turn_arounds and self.is_turn_around_scan:
            start_index = self.scan_stops_[self.current_scan_ - 1]
            stop_index = self.scan_starts_[self.current_scan_]
        else:
            start_index = self.scan_starts_[self.current_scan_]
            stop_index = self.scan_stops_[self.current_scan_]

        sframe = core.G3Frame(core.G3FrameType.Scan)
        ob_start_time = core.G3Time()
        ob_stop_time = core.G3Time()
        ob_start_time.mjd = self.idf_file_['/data/header/start_date/date'].value
        ob_stop_time.mjd = self.idf_file_['/data/header/stop_date/date'].value
        ob_n_samples = self.idf_file_['/data/header/n_samples'].value
        sample_dt = (ob_stop_time.time-ob_start_time.time)/ob_n_samples
        start = core.G3Time(ob_start_time.time + sample_dt*start_index)
        stop = core.G3Time(ob_start_time.time + sample_dt*stop_index)
        
        if self.load_bolo_data_:
            ts = core.G3TimestreamMap()
            for i in range(len(self.bolo_ids_)):
                if self.preload_data_:
                    uninit = self.bolo_data_storage_[i,start_index:stop_index]
                else:
                    uninit = self.idf_file_['/data/bolodata_array'][i,start_index:stop_index].astype('float64')
                t = core.G3Timestream(uninit)
                t.units = core.G3TimestreamUnits.Tcmb
                ts[self.bolo_ids_[i]] = t
            ts.start = start
            ts.stop = stop
            sframe['CalTimestreams'] = ts

        sframe['ScanIsBad'] = int(self.scan_flags_[self.current_scan_])

        sframe['BoresightRa'] = core.G3Timestream(self.idf_file_['/data/antenna/ra'][start_index:stop_index].astype('float64') * core.G3Units.deg)
        sframe['BoresightRa'].start = start
        sframe['BoresightRa'].stop = stop
        sframe['BoresightDec'] = core.G3Timestream(self.idf_file_['/data/antenna/dec'][start_index:stop_index].astype('float64') * core.G3Units.deg)
        sframe['BoresightDec'].start = start
        sframe['BoresightDec'].stop = stop

        sframe['BoresightAz'] = core.G3Timestream(self.idf_file_['/data/antenna/track_actual/data'][0][start_index:stop_index].astype('float64') * core.G3Units.deg)
        sframe['BoresightAz'].start = start
        sframe['BoresightAz'].stop = stop
        sframe['BoresightEl'] = core.G3Timestream(self.idf_file_['/data/antenna/track_actual/data'][1][start_index:stop_index].astype('float64') * core.G3Units.deg)
        sframe['BoresightEl'].start = start
        sframe['BoresightEl'].stop = stop

        sframe['TimestreamWeights'] = core.G3MapDouble()
        tmp_tswgt = self.idf_file_['/data/timestream_weights/'][:] / 1e6
        for i in range(len(self.bolo_ids_)):
            sframe['TimestreamWeights'][self.bolo_ids_[i]] = tmp_tswgt[i]
        sframe["ScanNumber"] = int(self.current_scan_)

        scan_flags =  sptpolflaggingutils.get_scan_flags( self.scan_ts_flag_information_[ self.current_scan_, :],
                                                          self.bolo_ids_, 
                                                          enforce_partner_good = self.enforce_partner_good_, 
                                                          bolo_partner_ids = self.bolo_partner_ids_)
        scan_flags.update(self.constant_flags)
        sframe['Flags'] = sptpolflaggingutils.convert_flag_dic_to_g3map(scan_flags)
        
        
        if ((self.include_turn_arounds and not self.is_turn_around_scan) or
            (not self.include_turn_arounds)):
            self.current_scan_ += 1

        sframe['Turnaround'] = self.is_turn_around_scan
        if self.include_turn_arounds:
            self.is_turn_around_scan = not self.is_turn_around_scan

        frames.append(sframe)

        return frames


class MultipleIdfReader(object):
    def __init__(self, fn_lst, preload_data = False, 
                 include_turn_arounds = False, 
                 load_bolo_data = True,       
                 ignore_ts_flags = [], invert_ts_flags = [],
                 ignore_bolo_flags = ['has_time_const', 'good_angle_fit', 
                                      'good_xpol_fit'],
                 invert_bolo_flags = ['has_pointing', 'has_polcal'],
                 enforce_partner_good = True):

        self.enforce_partner_good_ = enforce_partner_good
        self.ignore_ts_flags_ = ignore_ts_flags
        self.invert_ts_flags_ = invert_ts_flags
        self.ignore_bolo_flags_ = ignore_bolo_flags
        self.invert_bolo_flags_ = invert_bolo_flags

        assert(len(fn_lst) > 0)
        self.fn_lst_ = fn_lst
        self.idf_index_ = 0

        self.preload_data_ = preload_data 
        self.include_turn_arounds_ = include_turn_arounds
        self.load_bolo_data_ = load_bolo_data
        self.current_idf_reader_ = DirectIdfReader(
            self.fn_lst_[self.idf_index_], preload_data=self.preload_data_,
            include_turn_arounds=self.include_turn_arounds_,
            load_bolo_data=self.load_bolo_data_, 
            ignore_ts_flags = self.ignore_ts_flags_, 
            invert_ts_flags = self.invert_ts_flags_,
            ignore_bolo_flags = self.ignore_bolo_flags_,
            invert_bolo_flags = self.invert_bolo_flags_,
            enforce_partner_good = self.enforce_partner_good_)
    def __call__(self, frame):
        od = self.current_idf_reader_(frame)
        if not od:
            print('advancing', self.idf_index_)
            self.idf_index_ += 1
            if self.idf_index_ < len(self.fn_lst_):
                self.current_idf_reader_ = DirectIdfReader(
                    self.fn_lst_[self.idf_index_], preload_data=self.preload_data_,
                    include_turn_arounds=self.include_turn_arounds_,
                    load_bolo_data=self.load_bolo_data_)
                print('done advancing')
                return self.current_idf_reader_(frame)
            else:
                return od 
        else:
            return od



