import numpy
from spt3g import core

@core.indexmod
def ExtractCalibratorPhase(fr, Calibrator='CalibratorOn', InputTimestreams='CalTimestreams', Output='CalibratorPhase', PhaseHeadroom=2.1):
    '''
    Writes an integer into the frame containing the number of samples
    by which the calibrator phase in the average detector timestream is
    offset from the calibrator timestream. Searches up to PhaseHeadroom times
    the apparent half-width of the calibrator pulses.
    '''

    if fr.type != core.G3FrameType.Scan or Calibrator not in fr:
        return

    cal = fr[Calibrator]

    phasing = []
    ntransitions = numpy.sum(numpy.abs(numpy.diff(cal)))
    if ntransitions < 2: # If no transitions, skip this scan
        return

# kludge for telescope-moving spikes at beginning and end
#    phaselength = int(PhaseHeadroom*len(cal)*1./ntransitions)
    phaselength = 100
    avg = numpy.nanmedian(fr[InputTimestreams].values(), axis=0)
    internaldata = avg[phaselength:-phaselength]

    # Search for the maximum S/N offset. Is there a reason we want fractional
    # offsets here?
    for i in range(0, phaselength):
        calsub = cal[phaselength + i:-phaselength + i]
        phasing.append(numpy.mean(internaldata[calsub == 1]) - numpy.mean(internaldata[calsub == 0]))
    phase = numpy.argmax(phasing)

    fr[Output] = phase
    print('Phase is: '+str(phase)+' samples.\n')

@core.indexmod
def MakeSinCosFromSync(frame, Calibrator='CalibratorOn'):
    '''
    Make sine and cosine from fundamental of calibrator sync signal.
    '''

    if frame.type != core.G3FrameType.Scan:
        return

    cal = frame[Calibrator]
    calarr = numpy.asarray(cal)
    npts = len(calarr)
    isodd = False
    npts2 = npts
    if numpy.mod(npts,2) == 1:
        isodd = True
        calarr = numpy.append(calarr,calarr[npts-1])
        npts2 += 1
    calf = numpy.fft.fft(calarr - numpy.mean(calarr))
    whmax = numpy.argmax(numpy.abs(calf[0:len(calf)/2]))
    cosf = calf.copy()
    cosf[:] = 0.
    cosf[whmax] = calf[whmax]
    cosf[npts2 - whmax] = calf[npts2 - whmax]
    coscal = numpy.real(numpy.fft.ifft(cosf))
    coscal /= numpy.max(numpy.abs(coscal))
    sinf = cosf.copy()
    sinf[whmax] = numpy.complex(numpy.imag(cosf[whmax]),-numpy.real(cosf[whmax]))
    sinf[npts2 - whmax] = numpy.conj(sinf[whmax])
    sincal = numpy.real(numpy.fft.ifft(sinf))
    sincal /= numpy.max(numpy.abs(sincal))

    return coscal, sincal

@core.indexmod
def FitCalibratorResponse(frame, Calibrator='CalibratorOn', InputTimestreams='CalTimestreams', Phase='CalibratorPhase', Output='CalibratorResponse', CalcPerBoloPhase=False, PerBoloPhase='PerBoloPhase'):
    '''
    Fit for the response to the calibrator at some previously-determined phase in all
    detectors. Output is a mapping from detector ID to a floating point number in the
    same units as the input timestreams with the best-fit peak-to-peak amplitude of
    each detector's response to the calibrator.
    '''
   
    if frame.type != core.G3FrameType.Scan or Phase not in frame:
        return

    cal = frame[Calibrator]
    phase = frame[Phase]
    calresponse = core.G3MapDouble()

    if CalcPerBoloPhase:           # only use this for cal sweep tau analysis
        coscal, sincal = MakeSinCosFromSync(frame, Calibrator=Calibrator)
        calphase = core.G3MapDouble()
        for bolo,ts in frame[InputTimestreams].iteritems():
            response_i = numpy.mean(ts*coscal)
            response_q = numpy.mean(ts*sincal)
            calresponse[bolo] = numpy.sqrt(response_i**2 + response_q**2)
            calphase[bolo] = numpy.arctan2(response_q, response_i)
        frame[Output] = calresponse
        frame[PerBoloPhase] = calphase
    else:
        for bolo,ts in frame[InputTimestreams].iteritems():
            if phase != 0:
                subts = ts[:-phase]
            else:
                subts = ts
            calresponse[bolo] = numpy.mean(subts[cal[phase:] == 1.]) - numpy.mean(subts[cal[phase:] == 0.])
            
        frame[Output] = calresponse

@core.indexmod
def CalibratorResponseNoise(frame, Calibrator='CalibratorOn', InputTimestreams='CalTimestreams', Phase='CalibratorPhase', Output='CalibratorResponseNoise', Chunks=10):
    '''
    Fit for the response to the calibrator at some previously-determined phase in all
    detectors. Output is a mapping from detector ID to a floating point number in the
    same units as the input timestreams with the best-fit peak-to-peak noise amplitude of
    each detector's response to the calibrator.

    Operates by doing the normal fit (FitCalibratorResponse) on pieces of the data
    (Chunks), looking at the spread, and dividing by sqrt(Chunks)
    '''
   
    if frame.type != core.G3FrameType.Scan or Phase not in frame:
        return

    cal = frame[Calibrator]
    phase = frame[Phase]

# kludge for telescope-moving spikes at beginning and end
#    chunklen = len(cal) / Chunks
    chunklen = (len(cal)-200) / Chunks

    calresponse = {}
    for bolo in frame[InputTimestreams].keys():
        calresponse[bolo] = []
    for chunk in range(Chunks):
#        start = chunklen*chunk
        start = chunklen*chunk + 100
        stop = chunklen*(chunk + 1) + 100
        for bolo,ts in frame[InputTimestreams].iteritems():
            subts = ts[start:stop-phase]
            subcal = cal[start+phase:stop]
            calresponse[bolo].append(numpy.mean(subts[subcal == 1.]) - numpy.mean(subts[subcal == 0.]))

    outcalresponse = core.G3MapDouble()
    for bolo in frame[InputTimestreams].keys():
        outcalresponse[bolo] = numpy.std(calresponse[bolo])/numpy.sqrt(Chunks)

    frame[Output] = outcalresponse

@core.indexmod
def CalibratorResponseSN(frame, Input='CalibratorResponse', Noise='CalibratorResponseNoise', Output='CalibratorResponseSN'):
    '''
    Calculate the signal-to-noise of the calibrator response given the response and a
    noise estimate calculated previously.
    '''

    if frame.type != core.G3FrameType.Scan or Input not in frame:
        return

    calresponse = frame[Input]
    calnoise = frame[Noise]

    sn = core.G3MapDouble()
    for bolo in calresponse.keys():
        sn[bolo] = numpy.abs(calresponse[bolo])/calnoise[bolo]

    frame[Output] = sn

@core.pipesegment
def AnalyzeCalibratorData(pipe, Output='CalibratorResponse', Input='CalTimestreams', DropNoise=True, NoiseChunks=10, CalcPerBoloPhase=False, PerBoloPhase='PerBoloPhase'):
    '''
    Analyze the calibrator data in the frame. The main response dictionary will be in
    Output, with the signal-to-noise estimate in OutputSN. If DropNoise is False,
    the noise estimate will stored as OutputNoise. Since this can be calculated
    trivially from the response and the signal-to-noise, DropNoise is True by default.
    '''

    pipe.Add(ExtractCalibratorPhase, InputTimestreams=Input, Output=Output + 'Phase')
    pipe.Add(FitCalibratorResponse, InputTimestreams=Input, Output=Output, Phase=Output + 'Phase', CalcPerBoloPhase=CalcPerBoloPhase, PerBoloPhase=PerBoloPhase)
    pipe.Add(CalibratorResponseNoise, InputTimestreams=Input, Output=Output + 'Noise', Phase=Output + 'Phase', Chunks=NoiseChunks)
    pipe.Add(CalibratorResponseSN, Input=Output, Noise=Output + 'Noise', Output=Output + 'SN')
    todel = [Output + 'Phase']
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

            # accumulate calibrator frames
            cframe = core.G3Frame(core.G3FrameType.Calibration)
            cframe[self.cal_prefix] = frame[self.cal_prefix]
            cframe[self.cal_prefix + 'SN'] = frame[self.cal_prefix + 'SN']
            self.cal_frames.append(cframe)
            self.bs_el.append(frame[self.bs_el_key].mean())
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
            for k in self.cal_frames[0][self.cal_prefix]:
                cal_data[k] = core.G3VectorDouble(
                    [fr[self.cal_prefix][k] for fr in self.cal_frames])
                cal_data_sn[k] = core.G3VectorDouble(
                    [fr[self.cal_prefix + 'SN'][k] for fr in self.cal_frames])
            cframe[self.cal_prefix + 'Interp'] = cal_data
            cframe[self.cal_prefix + 'InterpSN'] = cal_data_sn

            # return the cal-el frame
            self._clear()
            return [cframe, frame]
