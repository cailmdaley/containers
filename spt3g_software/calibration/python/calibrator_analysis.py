import numpy
from spt3g import core, dfmux, todfilter
from .elnod_analysis import RotateIQ

def taper_window(length= 500, taper_length=100):
    window= numpy.ones(length)
    dist_to_edge=numpy.array([min(a) for a in zip(numpy.arange(length)-0, length-numpy.arange(length)-1)])
    ind_taper=numpy.where(dist_to_edge<=taper_length)
    window[ind_taper]= 0.5*window[ind_taper]- 0.5*numpy.cos(dist_to_edge[ind_taper]*numpy.pi/taper_length)
    norm_factor= numpy.mean(window)
    return numpy.array(window), norm_factor

@core.indexmod
class CalibratorRotateIQ(RotateIQ):
    '''
    Rotate I and Q data by the optimal responsitivity phase angle.
    If q_rotate_key is not given, immediately returns the input I
    data if the Q data have been stripped. If destructive is set to
    True, deletes the input timestreams. If normalize is set to
    True, normalize the eigenvector and ensure that it preserves
    the sign of I on rotation.

    Sam Guns, Jan 29 2019. Copied from elnod_analysis.py and changed default key names.
    ASR, Feb 7 2019. Refactored to avoid code duplication.
    '''

    def __init__(self,
                 i_vector_key='CalibratorResponseI',
                 q_vector_key='CalibratorResponseQ',
                 i_data_key='RawTimestreams_I',
                 q_data_key='RawTimestreams_Q',
                 i_rotate_key='RawTimestreams',
                 q_rotate_key=None,
                 destructive=False,
                 normalize=True):

        super(CalibratorRotateIQ, self).__init__(
            i_vector_key=i_vector_key,
            q_vector_key=q_vector_key,
            i_data_key=i_data_key,
            q_data_key=q_data_key,
            i_rotate_key=i_rotate_key,
            q_rotate_key=q_rotate_key,
            destructive=destructive,
            normalize=normalize)

@core.indexmod
def CreateCalibratorOn(fr, InputTimestreams='CalTimestreams',
                       Output='CalibratorOn', 
                       FlagKey='CalibratorOnIsFake'):
    '''
    If calibrator sync signal is missing, synthesize a fake one from
    bolometer data. This is inefficient (executes many steps that will
    be repeated later in the pipeline) and is only to be used in an
    emergency when the sync signal is broken or didn't get written to
    file.
    '''

    if fr.type != core.G3FrameType.Scan:
        return

    stds = {}
    for name in fr[InputTimestreams].keys():
        if numpy.isfinite(numpy.sum(fr[InputTimestreams][name])): 
            if numpy.max(numpy.abs(fr[InputTimestreams][name])) > 0.:
                stds[name] = numpy.std(fr[InputTimestreams][name])
            else:
                stds[name] = numpy.nan
    medstds = numpy.nanmedian(list(stds.values()))
    avg = numpy.zeros(len(fr[InputTimestreams][name]))
    for name in fr[InputTimestreams].keys():
        if numpy.isfinite(numpy.sum(fr[InputTimestreams][name])): 
            if numpy.max(numpy.abs(fr[InputTimestreams][name])) > 0. and stds[name] < medstds*10.:
                avg += fr[InputTimestreams][name]
    avg2 = avg - numpy.mean(avg)
    import scipy.fftpack
    sd_est = scipy.fftpack.fft(avg2 * numpy.hanning(len(avg2)))
    asd_est = numpy.abs(sd_est)
    asd_freqs = scipy.fftpack.fftfreq(n = len(avg2),
                                    d =1.0/(fr[InputTimestreams][fr[InputTimestreams].keys()[0]].sample_rate/core.G3Units.Hz))
    inds = numpy.where(asd_freqs < 2)
    asd_est[inds] = 0
    cal_freq_ind = numpy.where(asd_est == numpy.max(asd_est))
    cal_freq = float(asd_freqs[cal_freq_ind]*core.G3Units.Hz)
    sd_fund = sd_est.copy()
    sd_fund[:] = 0.
    sd_fund[cal_freq_ind] = sd_est[cal_freq_ind]
    fake_cal_sine = numpy.real(scipy.fftpack.ifft(sd_fund))
    fake_cal_sine -= numpy.min(fake_cal_sine)
    fake_cal_sine /= numpy.max(fake_cal_sine)
    fake_cal = numpy.zeros(len(fake_cal_sine))
    fake_cal[numpy.where(fake_cal_sine > 0.5)] = 1.
    fake_cal[numpy.where(fake_cal_sine <= 0.5)] = 0.

    cal = core.G3Timestream(fake_cal)
    ts1 = fr[InputTimestreams][fr[InputTimestreams].keys()[0]]
    cal.start = ts1.start
    cal.stop = ts1.stop

    fr[Output] = cal
    fr[FlagKey] = True

@core.indexmod
def ExtractCalibratorPhase(fr, Calibrator='CalibratorOn',
                           InputTimestreams='CalTimestreams',
                           Output='CalibratorPhase', PhaseHeadroom=2.1):
    '''
    Writes an integer into the frame containing the number of samples
    by which the calibrator phase in the average detector timestream is
    offset from the calibrator timestream. Searches up to PhaseHeadroom times
    the apparent half-width of the calibrator pulses.
    '''

    if fr.type != core.G3FrameType.Scan:
        return

    calthere = True
    if Calibrator not in fr:# or (numpy.asarray(fr[Calibrator]) == 0).all():
        # make fake sync signal if none present
        calthere = False
        if Calibrator in fr:
            del fr[Calibrator]
        createcal = CreateCalibratorOn(fr,InputTimestreams=InputTimestreams,Output=Calibrator)

    cal = fr[Calibrator]

    phasing = []
    ntransitions = numpy.sum(numpy.abs(numpy.diff(cal)))
    if ntransitions < 2: # If no transitions, skip this scan
        return

    phaselength = int(PhaseHeadroom*len(cal)*1./ntransitions)
    stds = {}
    for name in fr[InputTimestreams].keys():
        if numpy.isfinite(numpy.sum(fr[InputTimestreams][name])): 
            if numpy.max(numpy.abs(fr[InputTimestreams][name])) > 0.:
                stds[name] = numpy.std(fr[InputTimestreams][name])
            else:
                stds[name] = numpy.nan
    medstds = numpy.nanmedian(list(stds.values()))
    avg = numpy.zeros(len(cal))
    for name in fr[InputTimestreams].keys():
        if numpy.isfinite(numpy.sum(fr[InputTimestreams][name])): 
            if numpy.max(numpy.abs(fr[InputTimestreams][name])) > 0. and stds[name] < medstds*10.:
                avg += fr[InputTimestreams][name]

    internaldata = avg[phaselength:-phaselength]

    # Search for the maximum S/N offset. Is there a reason we want fractional
    # offsets here?
    for i in range(0, phaselength):
        calsub = cal[phaselength + i:-phaselength + i]
        phasing.append(numpy.mean(internaldata[calsub == 1]) - numpy.mean(internaldata[calsub == 0]))
    phase = int(numpy.argmax(phasing))

    fr[Output] = phase

@core.indexmod
def ExtractCalibratorFrequency(fr, calibrator='CalibratorOn',
                               output='CalibratorFrequency',
                               use_fft_method = True):
    '''
    Extracts the calibrator frequency.

    if use_fft_method:
        estimates the ASD of the calibrator signal and finds the maximum
    else:
        counts the step transitions the calibrator signal does and estimates
        the frequency from that.
    '''
    if fr.type != core.G3FrameType.Scan or calibrator not in fr:
        return
    cal = fr[calibrator]
    if (numpy.asarray(cal) == 0).all():
        core.log_fatal('Cannot determine calibrator frequency from empty timestream.')
    if use_fft_method:
        import scipy.fftpack
        asd_est = numpy.abs(scipy.fftpack.fft(cal * numpy.hanning(len(cal))))
        asd_freqs = scipy.fftpack.fftfreq(n = len(cal),
                                          d =1.0/(cal.sample_rate/core.G3Units.Hz))
        inds = numpy.where(asd_freqs < 2)
        asd_est[inds] = 0
        cal_freq = float(asd_freqs[numpy.argmax(asd_est)] * core.G3Units.Hz)
    else:
        ntransitions = numpy.sum(numpy.abs(numpy.diff(cal)))
        cal_freq = float(ntransitions/2.0 * cal.sample_rate)
    fr[output] = cal_freq

@core.indexmod
def FitCalibratorResponse(frame, Calibrator='CalibratorOn',
                          InputTimestreams='CalTimestreams',
                          Phase='CalibratorPhase', Output='CalibratorResponse', 
                          taper_length=100):
    '''
    Fit for the response to the calibrator at some previously-determined phase
    in all detectors. Output is a mapping from detector ID to a floating point
    number in the same units as the input timestreams with the best-fit
    peak-to-peak amplitude of each detector's response to the calibrator.
    '''
   
    if frame.type != core.G3FrameType.Scan or Phase not in frame:
        return

    cal = frame[Calibrator]
    phase = frame[Phase]

    calresponse = core.G3MapDouble()
    for bolo,ts in frame[InputTimestreams].iteritems():
        if phase != 0:
            subts = ts[:-phase]
        else:
            subts = ts
        len_data= len(subts)
        window, norm_factor=taper_window(len_data, taper_length)
        subts=subts-numpy.mean(subts)
        subts=subts*window/norm_factor
        calresponse[bolo] = numpy.mean(subts[cal[phase:] == 1.]) - numpy.mean(subts[cal[phase:] == 0.])

    frame[Output] = calresponse

@core.indexmod
def CalibratorResponseNoise(frame, Calibrator='CalibratorOn',
                            InputTimestreams='CalTimestreams',
                            Phase='CalibratorPhase',
                            Output='CalibratorResponseNoise',
                            Chunks=10, taper_length=100):
    '''
    Fit for the response to the calibrator at some previously-determined phase
    in all detectors. Output is a mapping from detector ID to a floating point
    number in the same units as the input timestreams with the best-fit
    peak-to-peak noise amplitude of each detector's response to the calibrator.

    Operates by doing the normal fit (FitCalibratorResponse) on pieces of the
    data (Chunks), looking at the spread, and dividing by sqrt(Chunks)
    '''
   
    if frame.type != core.G3FrameType.Scan or Phase not in frame:
        return

    cal = frame[Calibrator]
    phase = frame[Phase]

    chunklen = len(cal) // Chunks

    calresponse = {}
    for bolo in frame[InputTimestreams].keys():
        calresponse[bolo] = []
    for chunk in range(Chunks):
        start = chunklen*chunk
        stop = chunklen*(chunk + 1)
        for bolo,ts in frame[InputTimestreams].iteritems():
            subts = ts[start:stop-phase]
            len_data=len(subts)
            window, norm_factor=taper_window(len_data, taper_length)
            subts= subts-numpy.mean(subts)
            subts=subts*window/norm_factor
            subcal = cal[start+phase:stop]
            calresponse[bolo].append(numpy.mean(subts[subcal == 1.]) - numpy.mean(subts[subcal == 0.]))

    outcalresponse = core.G3MapDouble()
    for bolo in frame[InputTimestreams].keys():
        outcalresponse[bolo] = numpy.std(calresponse[bolo])/numpy.sqrt(Chunks)

    frame[Output] = outcalresponse

@core.indexmod
def CalibratorResponseSN(frame, Input='CalibratorResponse',
                         Noise='CalibratorResponseNoise',
                         Output='CalibratorResponseSN'):
    '''
    Calculate the signal-to-noise of the calibrator response given the
    response and a noise estimate calculated previously.
    '''

    if frame.type != core.G3FrameType.Scan or Input not in frame:
        return

    calresponse = frame[Input]
    calnoise = frame[Noise]

    sn = core.G3MapDouble()
    for bolo in calresponse.keys():
        if numpy.isfinite(calresponse[bolo]) and numpy.isfinite(calnoise[bolo]) and calnoise[bolo] != 0.:
            sn[bolo] = numpy.abs(calresponse[bolo])/calnoise[bolo]
        else:
            sn[bolo] = numpy.nan

    frame[Output] = sn

@core.pipesegment
def CalibratorResponseRfrac(pipe, InputTimestreams='CalTimestreams',
                            Output='CalibratorResponseRfrac'):

    res_ts_key = '_ResTimestreams'
    pipe.Add(dfmux.ConvertTimestreamUnits, Input=InputTimestreams,
             Output=res_ts_key, Units=core.G3TimestreamUnits.Resistance)

    @core.scan_func_cache_data(bolo_props='BolometerProperties',
                               wiring_map='WiringMap', hk='DfMuxHousekeeping')
    def CalRespRfracInner(frame, bolo_props=None, wiring_map=None, hk=None):
        if frame.type != core.G3FrameType.Scan:
            return

        rtsm = frame[res_ts_key]
        frame[Output] = core.G3MapDouble()

        for k in rtsm.keys():
            chan_hk = dfmux.HousekeepingForBolo(hk, wiring_map, k)
            Rnorm = chan_hk.rnormal
            res = todfilter.get_mean(rtsm[k])
            frame[Output][k] = res/Rnorm

    pipe.Add(CalRespRfracInner)
    pipe.Add(core.Delete, keys=[res_ts_key], type=core.G3FrameType.Scan)

@core.pipesegment
def AnalyzeCalibratorData(pipe, Output='CalibratorResponse',
                          Input='CalTimestreams', Input_Q = None,
                          DropNoise=True, NoiseChunks=10, Do_IQ_Rotation = False):
    '''
    Analyze the calibrator data in the frame. The main response dictionary will
    be in Output, with the signal-to-noise estimate in OutputSN. If
    DropNoise is False, the noise estimate will stored as OutputNoise. Since
    this can be calculated trivially from the response and the signal-to-noise,
    DropNoise is True by default.
    '''
    from spt3g import todfilter

    pipe.Add(CalibratorResponseRfrac, InputTimestreams=Input,
             Output=Output+'Rfrac')

    if Do_IQ_Rotation and Input_Q is None:
        core.log_fatal("Need to specify Q input timestream to do IQ rotation.")
    if Do_IQ_Rotation and Input_Q is not None:
        Input_I = Input +'I'
        Input_Rotated = Input + 'Rotated'
        Output_I = Output + 'I'
        Output_Q = Output + 'Q' 
        pipe.Add(todfilter.MaskedPolyHpf, poly_order=40, in_ts_map_key = Input_I, out_ts_map_key='PolyFiltered'+Input_I) 
        pipe.Add(todfilter.MaskedPolyHpf, poly_order=40, in_ts_map_key = Input_Q, out_ts_map_key='PolyFiltered'+Input_Q) 
        pipe.Add(ExtractCalibratorPhase, InputTimestreams='PolyFiltered'+Input_I, Output=Output_I + 'Phase')
        pipe.Add(FitCalibratorResponse, InputTimestreams='PolyFiltered'+Input_I, Output=Output_I, Phase=Output_I + 'Phase')
        pipe.Add(FitCalibratorResponse, InputTimestreams='PolyFiltered'+Input_Q, Output=Output_Q, Phase=Output_I + 'Phase') #use I phase for Q fitting
        pipe.Add(CalibratorRotateIQ, i_vector_key = Output_I, q_vector_key = Output_Q, i_data_key = Input_I, q_data_key = Input_Q, i_rotate_key = Input_Rotated) 

    if not Do_IQ_Rotation:
        Input_Rotated = Input

    pipe.Add(todfilter.MaskedPolyHpf, poly_order=40, in_ts_map_key=Input_Rotated, out_ts_map_key='PolyFiltered'+Input_Rotated)
    pipe.Add(ExtractCalibratorPhase, InputTimestreams='PolyFiltered'+Input_Rotated, Output=Output + 'Phase')
    pipe.Add(ExtractCalibratorFrequency, output = Output+'Frequency')
    pipe.Add(FitCalibratorResponse, InputTimestreams='PolyFiltered'+Input_Rotated, Output=Output, Phase=Output + 'Phase')
    pipe.Add(CalibratorResponseNoise, InputTimestreams='PolyFiltered'+Input_Rotated, Output=Output + 'Noise', Phase=Output + 'Phase', Chunks=NoiseChunks)
    pipe.Add(CalibratorResponseSN, Input=Output, Noise=Output + 'Noise', Output=Output + 'SN')
    todel = []
    if DropNoise:
        todel.append(Output + 'Noise')
    pipe.Add(core.Delete, keys=todel, type=core.G3FrameType.Scan)

@core.indexmod
class MakeCalibratorFrame(object):
    def __init__(self, cal_prefix='CalibratorResponse',
                 boresight_el_key='RawBoresightEl'):

        self._clear()
        self.cal_prefix = cal_prefix
        self.bs_el_key = boresight_el_key

    def _clear(self):
        self.bs_el = []
        self.cal_frames = []

    def __call__(self, frame):

        if frame.type == core.G3FrameType.Scan and self.cal_prefix in frame:

            cal_prefix_i = self.cal_prefix +'I'
            cal_prefix_q = self.cal_prefix +'Q'

            # accumulate calibrator frames
            cframe = core.G3Frame(core.G3FrameType.Calibration)
            cframe[self.cal_prefix] = frame[self.cal_prefix]
            cframe[self.cal_prefix + 'SN'] = frame[self.cal_prefix + 'SN']
            cframe[self.cal_prefix + 'Frequency'] = frame[self.cal_prefix + 'Frequency']
            cframe[self.cal_prefix + 'Phase'] = frame[self.cal_prefix + 'Phase']
            cframe[self.cal_prefix + 'Rfrac'] = frame[self.cal_prefix + 'Rfrac']
            self.cal_frames.append(cframe)
            self.bs_el.append(frame[self.bs_el_key].mean())

            if cal_prefix_i in frame:
                cframe[cal_prefix_i] = frame[cal_prefix_i]
            if cal_prefix_q in frame:
                cframe[cal_prefix_q] = frame[cal_prefix_q]

            return

        elif frame.type == core.G3FrameType.EndProcessing:

            # no calibrator data found
            if len(self.bs_el) == 0:
                return

            # there's only one calibrator frame in this dataset, so just return it
            if len(self.bs_el) == 1:
                cframe = self.cal_frames[0]
                self._clear()
                return [cframe, frame]

            # found more than one calibrator frame, so probably a cal-el dataset
            # so we create a frame with vectors for interpolation
            cframe = core.G3Frame(core.G3FrameType.Calibration)

            # reference elevations
            cframe[self.cal_prefix + 'El'] = core.G3VectorDouble(self.bs_el)

            # response data map
            # is there a more efficient way to do this without looping over keys?
            cal_data = core.G3MapVectorDouble()
            cal_data_sn = core.G3MapVectorDouble()
            for k in self.cal_frames[0][self.cal_prefix].keys():
                print(repr(k))
                cal_data[k] = core.G3VectorDouble(
                    [fr[self.cal_prefix][k] for fr in self.cal_frames])
                cal_data_sn[k] = core.G3VectorDouble(
                    [fr[self.cal_prefix + 'SN'][k] for fr in self.cal_frames])
            cframe[self.cal_prefix + 'Interp'] = cal_data
            cframe[self.cal_prefix + 'InterpSN'] = cal_data_sn

            # return the cal-el frame
            self._clear()
            return [cframe, frame]
