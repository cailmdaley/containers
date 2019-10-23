import datetime as dt
from spt3g import core, std_processing
import numpy as np
from spt3g.util import extractdata


def dt_to_spt(dt_time):
    gcp_fmt = '%y%m%d %H:%M:%S'
    return core.G3Time(dt_time.strftime(gcp_fmt))


def spt_to_dt(spt_time):
    g3_fmt = '%d-%b-%Y:%H:%M:%S'
    return dt.datetime.strptime(spt_time.Summary().split('.')[0], g3_fmt)

def movingAvg(data, window):
    data = np.convolve(data, np.ones((window,))/window, mode='same')
    return data



class Cycle:
    spt_fmt = '%y%m%d %H:%M:%S'
    def __init__(self, num=0, cadence='Q', dither=0, 
                 cycle_start='', 
                 cycle_end = '', 
                 schedule_end = ''):
        
        self._cycle_num = num
        self._cadence = cadence.upper()
        self._dither = dither
        self._cycle_start = cycle_start
        self._cycle_end = cycle_end
        self._schedule_end = schedule_end

        
    @property
    def cycle_num(self):
        return self._cycle_num
    @cycle_num.setter
    def cycle_num(self, num):
        self._cycle_num = num

   
    @property
    def cycle_start(self):
        return self._cycle_start
    @cycle_start.setter
    def cycle_start(self, time_string):
        self._cycle_start = time_string

    @property
    def dt_cycle_start(self):
        return dt.datetime.strptime(self.cycle_start, self.spt_fmt)

    @property
    def spt_cycle_start(self):
        return core.G3Time(self.cycle_start)

    @property
    def cycle_end(self):
        return self._cycle_end
    @cycle_end.setter
    def cycle_end(self, time_string):
        self._cycle_end = time_string

    @property
    def dt_cycle_end(self):
        return dt.datetime.strptime(self.cycle_end, self.spt_fmt)

    @property
    def spt_cycle_end(self):
        return core.G3Time(self.cycle_end)

    @property
    def schedule_end(self):
        return self._schedule_end
    @schedule_end.setter
    def schedule_end(self, time_string):
        self._schedule_end = time_string

    @property
    def dt_schedule_end(self):
        return dt.datetime.strptime(self.schedule_end, self.spt_fmt)
    
    @property
    def spt_schedule_end(self):
        return core.G3Time(self.schedule_end)

    @property
    def hold_time(self):
        return self.dt_uc_blown - self.dt_cycle_end

    @property
    def total_time(self):
        return self.dt_uc_blown - self.dt_cycle_start

    @property
    def efficiency(self):
        return self.hold_time/self.total_time
    
    @property
    def cycle_time(self):
        return self.dt_cycle_end - self.dt_cycle_start

    @property
    def cadence(self):
        return self._cadence
    @cadence.setter
    def cadence(self, cadence):
        self._cadence = cadence.upper()

    @property
    def dither(self):
        return self._dither
    @dither.setter
    def dither(self, dither):
        self._dither = dither


    @property
    def uc_blown(self):
        return self._uc_blown
    @uc_blown.setter
    def uc_blown(self,uc_blown):
        self._uc_blown = uc_blown

    @property
    def ic_blown(self):
        return self._ic_blown
    @ic_blown.setter
    def ic_blown(self,ic_blown):
        self._ic_blown = ic_blown

    @property
    def dt_uc_blown(self):
        return dt.datetime.strptime(self._uc_blown, self.spt_fmt)

    @property
    def dt_ic_blown(self):
        return dt.datetime.strptime(self._ic_blown, self.spt_fmt)

    @property
    def spt_ic_blown(self):
        return core.G3Time(self._ic_blown)

    @property
    def spt_uc_blown(self):
        return core.G3Time(self._uc_blown)


def archive_extract(keys, start_time, stop_time):
    spt_fmt = '%y%m%d %H:%M:%S'
    arcdir = '/spt_data/arc/'

    if isinstance(start_time, str):
        try:
            start_time = dt.datetime.strptime(start_time, spt_fmt)
        except:
            print('Start time must be formatted like the logtable or a datetime object')
        
        start_time = dt_to_spt(start_time)
    
    if isinstance(stop_time, str):
        try:
            stop_time = dt.datetime.strptime(stop_time, spt_fmt)
        except:
            print('Stop time must be formatted like the logtable or a datetime object')
        
        stop_time = dt_to_spt(stop_time)
    
    if isinstance(start_time, dt.datetime):
        start_time = dt_to_spt(start_time)
    if isinstance(stop_time, dt.datetime):
        stop_time = dt_to_spt(stop_time)
    
    mult_accumulator = extractdata.MultiAccumulator(keys)
    pipe1 = core.G3Pipeline()
    pipe1.Add(std_processing.ARCTimerangeReader, start_time = start_time,
              stop_time = stop_time, basedir=arcdir)
    pipe1.Add(mult_accumulator)
    pipe1.Run()
    return mult_accumulator.extract_values()
