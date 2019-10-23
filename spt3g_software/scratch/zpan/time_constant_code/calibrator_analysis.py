import numpy
from spt3g import core
from scipy.optimize import curve_fit

# This is the calibrator analysis code for cal sweeps
# The autoprocessing version loses out-of-phase parts,
# which is bad at high freqs.

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

    if fr.type != core.G3FrameType.Scan or Calibrator not in fr:
        return

    cal = fr[Calibrator]

    phasing = []
    ntransitions = numpy.sum(numpy.abs(numpy.diff(cal)))
    if ntransitions < 2: # If no transitions, skip this scan
        return

    phaselength = int(PhaseHeadroom*len(cal)*1./ntransitions)
    avg = numpy.nanmedian(fr[InputTimestreams].values(), axis=0)
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
    if use_fft_method:
        import scipy.fftpack
        asd_est = numpy.abs(scipy.fftpack.fft(cal *
                                           numpy.hanning(len(cal))))
        asd_freqs = scipy.fftpack.fftfreq(n = len(cal),
                                          d =1.0/(cal.sample_rate/core.G3Units.Hz))
        inds = numpy.where(asd_freqs < 2)
        asd_est[inds] = 0
        cal_freq = float(asd_freqs[ numpy.where(asd_est == numpy.max(asd_est))]*core.G3Units.Hz)
    else:
        ntransitions = numpy.sum(numpy.abs(numpy.diff(cal)))
        cal_freq = float(ntransitions/2.0 * cal.sample_rate)
    fr[output] = cal_freq

@core.indexmod
def ExtractAzEl(fr, OnlineBoresightAz='OnlineBoresightAz',
                OnlineBoresightEl='OnlineBoresightEl',
                output1='CalibratorResponseOnlineBoresightAzMedian',
                output2='CalibratorResponseOnlineBoresightElMedian'):
    '''
    Extracts the calibrator azimuth and elevation for elevation sweeps.
    '''
    if fr.type != core.G3FrameType.Scan or OnlineBoresightAz not in fr:
        return
    fr[output1]=numpy.median(numpy.asarray(fr[OnlineBoresightAz]))
    fr[output2]=numpy.median(numpy.asarray(fr[OnlineBoresightEl]))


@core.indexmod
def FitCalibratorResponse(frame, Calibrator='CalibratorOn',
                          InputTimestreams='CalTimestreams',
                          Phase='CalibratorPhase', Output='CalibratorResponse'):
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
        calresponse[bolo] = numpy.mean(subts[cal[phase:] == 1.]) - numpy.mean(subts[cal[phase:] == 0.])

    frame[Output] = calresponse

@core.indexmod
def CalibratorResponseNoise(frame, Calibrator='CalibratorOn',
                            InputTimestreams='CalTimestreams',
                            Phase='CalibratorPhase',
                            Output='CalibratorResponseNoise', Chunks=10):
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
        sn[bolo] = numpy.abs(calresponse[bolo])/calnoise[bolo]

    frame[Output] = sn


# below are some functions copied from sptpol_sofrware for the calibrator response analysis with individual phase information


def pad2power2len(n, factor=1.1):
    'return the power of two that is at least factor (10% by default) bigger.'
    return 2**int(numpy.ceil(numpy.log2(n*factor)))

    
def pad2power2(a, factor=1.1):
    '''
    return a version of 1D array a, 
    which is zero-padded at the front and back to a new length, 
    which is a power of two that is at least factor (10% by default) bigger.
    '''
    n = len(a)
    N = pad2power2len(n, factor)
    b = numpy.zeros(N)
    i = (N-n)//2  # starting index
    b[i:i+n] = a
    return b
    
def unpad2power2(b, n):
    '''
    undo pad2power2 by extracting the center of the 1D array b.
    n is the original length
    '''
    N = len(b)
    i = (N-n)//2  # starting index
    return b[i:i+n]


def do_padded_windowed_fft(a, window, npadded):
    '''
    Helper function to take windowed fft of a mean-subtracted version of a (usually caldata).
        a:           1D array like data.caldata.data
        window:      1D array the same length as a
        npadded:     number of samples in the padded array to be transformed.  (Should be a power of 2 for speed.)
    returns:
        fft:             1D complex array, which is real fft of a
    '''
    nomean = a - numpy.mean(a)
    padded = pad2power2(nomean*window)  # Multiply the window before the padding so that there's no sharp transition
    fft = numpy.fft.rfft(padded)
    assert(len(padded)==npadded)
    return fft
    


def frequency_and_ifft(frequencies, unmasked_fft, window, min_freq, linewidth):
    '''
    Helper function used on both cal and bolo fft.  Find weighted-average peak frequency and inverses the fft.
        frequencies:  list of fft frequencies
        unmasked_fft: real fft data, which we will copy and zero bins below min_freq and outside of linewidth beyond the peak frequency
        window:       a window function, used to determine the length of the ifft data, which then gets divided by this window.
        min_freq:     like 3.0, don't look for a peak below this
        linewidth:    look linewidth on either side of the peak frequency.
    return:
        fft:             1D complex array, which is a copy of unmasked_fft, but with everything away from peak set to 0
        ifft:            1D real array, which is inverse real fft of the masked fft
        frequency:       float, which is the weighted-average of the frequency hump.  Best extimate of actual frequency of calibrator.
    '''
    fft = unmasked_fft.copy()  # Don't destroy original.  So we can plot it.  Maybe we don't really need this.
    fft[frequencies<min_freq] = 0.  # Zero out low frequencies
    frequency_peak = frequencies[numpy.argmax(abs(fft))]
    
    fft[abs(frequencies-frequency_peak)>linewidth] = 0.  # Zero out all frequencies away from max.
    #fft[frequencies!=frequency_peak] = 0.  # Zero out all frequencies away from max.  Too narrow. Screws up phase
    weights = abs(fft)**2
    if sum(weights) > 0.:
        frequency = numpy.average(frequencies, weights=weights)  # a better estimate, weighted by the narrow power spectrum peak.
    else:
        frequency = 6.  # numpy.average throws ZeroDivisionError if Weights sum to zero.  Just pick something non-zero to return.
    ifft = unpad2power2(numpy.fft.irfft(fft), len(window))/window
    #ifft = numpy.fft.irfft(fft, n=len(window))/window
    return fft, ifft, frequency
    

def fitted_sine(t, a, frequency):
    '''
    Helper function.  Fits a cosine of given frequency*t to array a.
    returns fitted (cos,sin) arrays (same length as t and a).
    '''
    def cosfunc(t, phase):
        return numpy.cos(2.*numpy.pi*frequency*t + phase)
    def sinfunc(t, phase):
        return numpy.sin(2.*numpy.pi*frequency*t + phase)
    try:
        popt, pcov = curve_fit(cosfunc, t, a, p0=(0.0,))  # don't fit freq or amp.  just phase
    except RuntimeError as e:
        print(e)
        popt = [0]
    cal_phase = popt[0]
    return cosfunc(t, *popt), sinfunc(t, *popt)
    

def centnorm(a):
    '''
    Do a mean subtraction to center around zero, then normalize by the rms.
    A sin wave will come out going from -1 to +1
    '''
    return (a-numpy.mean(a))/numpy.std(a)/numpy.sqrt(2)


@core.indexmod
def MakeSinCosFromSync(frame, Calibrator='CalibratorOn'):
    '''
    Make sine and cosine from fundamental of calibrator sync signal.
    Output is a 
    '''

    if frame.type != core.G3FrameType.Scan:
        return

    cal = frame[Calibrator]
    calarr = numpy.asarray(cal)
    #window function for the calibration data
    window=numpy.hamming(len(calarr))
    #normalize the window function
    window/=numpy.std(window)
    #pad to powers of 2 to speed up the fft
    npadded=pad2power2len(len(window))
    #number of points in the real fft
    nrfft=npadded//2+1
    #frequency list
    frequencies=numpy.linspace(0.0, cal.sample_rate/core.G3Units.hz/2., num=nrfft, endpoint=True)
    #time stamps for the calibrator data array
    t=numpy.arange(len(calarr))/(cal.sample_rate/core.G3Units.hz)
    #take a fft for the calibration signal
    unmasked_cal_fft=do_padded_windowed_fft(calarr, window, npadded)
    cal_linewidth=0.05
    min_freq=3.
    #get the ifft for the calibrator data with components near the peak of the fft. 
    (cal_fft, cal_ifft, cal_freq)= frequency_and_ifft(frequencies, unmasked_cal_fft, window, min_freq, cal_linewidth)
    # fit to a cosine function to get the cos/sin components. 
    (coscal, sincal) = fitted_sine(t, centnorm(cal_ifft), cal_freq)
    coscal /= numpy.std(coscal)*2/numpy.sqrt(2)
    sincal /= numpy.std(sincal)*2/numpy.sqrt(2)
    return coscal, sincal

@core.indexmod
def FitCalibratorResponsePerBolo(frame, Calibrator='CalibratorOn',
                          InputTimestreams='CalTimestreams',
                          Output='CalibratorResponsePerBolo',
                          OutputI='CalibratorResponsePerBoloI',
                          OutputQ='CalibratorResponsePerBoloQ',
                          PerBoloPhase='CalibratorResponsePerBoloPhase', 
                          CalCos='CalCos', CalSin='CalSin'):
    '''
    Fit the detector response to a cosine wave and sine wave which are in 
    phase and out of phase to the calibrator LED signal.The LED being on 
    and off is not necessarily coincident with the power being turned on 
    and off. 
    '''
   
    if frame.type != core.G3FrameType.Scan:
        return

    cal = frame[Calibrator]
    calresponse = core.G3MapDouble()
    #This segment will fit for a cosine and a sine wave.
    coscal, sincal = MakeSinCosFromSync(frame, Calibrator=Calibrator)
    calphase = core.G3MapDouble()
    calresponseI = core.G3MapDouble()
    calresponseQ = core.G3MapDouble()
    for bolo,ts in frame[InputTimestreams].iteritems():
        zeromeants= numpy.array(ts)-numpy.mean(ts)
        response_i = numpy.mean(zeromeants*coscal)
        response_q = numpy.mean(zeromeants*sincal)
        calresponse[bolo] = numpy.sqrt(response_i**2 + response_q**2)
        #This number will be the phase of the time delay, in radians. 
        calphase[bolo] = numpy.abs(numpy.arctan2(response_i, response_q))
        calresponseI[bolo]= response_i
        calresponseQ[bolo]= response_q
        
    frame[Output] = calresponse
    frame[OutputI] = calresponseI
    frame[OutputQ] = calresponseQ
    #Here we take the absolute value because when fitting for a cosine function
    # the input could be both positive/negative. So sign is meaningless.
    frame[PerBoloPhase] = calphase
    #save the cos/sin waves extracted from the calibrator square wave
    frame[CalCos]=core.G3Timestream(coscal)
    frame[CalSin]=core.G3Timestream(sincal)

@core.indexmod
def CalibratorResponsePerBoloNoise(frame, Calibrator='CalibratorOn',
                            InputTimestreams='CalTimestreams', 
                            Output='CalibratorResponsePerBoloNoise', Chunks=10, 
                            CalCos='CalCos',CalSin='CalSin'):
    '''
    It's similar to CalibratorResponseNoise function, but for the per bolo method.
    '''
    if frame.type != core.G3FrameType.Scan:
        return
    cal = frame[Calibrator]
    coscal=numpy.asarray(frame[CalCos]) 
    sincal=numpy.asarray(frame[CalSin]) 
    chunklen = len(cal) // Chunks
    calresponse = {}
    calresponseI= {}
    calresponseQ= {}
    for bolo in frame[InputTimestreams].keys():
        calresponse[bolo] = []
        calresponseI[bolo]= []
        calresponseQ[bolo]= []
    for chunk in range(Chunks):
        start = chunklen*chunk
        stop = chunklen*(chunk + 1)
        for bolo,ts in frame[InputTimestreams].iteritems():
            subts = ts[start:stop]
            subcos=coscal[start:stop]
            subsin=sincal[start:stop]
            response_i=numpy.mean((numpy.array(subts)-numpy.mean(subts))*subcos)
            response_q=numpy.mean((numpy.array(subts)-numpy.mean(subts))*subsin)
            calresponse[bolo].append(numpy.sqrt(response_i**2 + response_q**2))
            calresponseI[bolo].append(response_i)
            calresponseQ[bolo].append(response_q)

    outcalresponse = core.G3MapDouble()
    outcalresponseI = core.G3MapDouble()
    outcalresponseQ = core.G3MapDouble()

    for bolo in frame[InputTimestreams].keys():
        outcalresponse[bolo] = numpy.std(calresponse[bolo])/numpy.sqrt(Chunks)
        outcalresponseI[bolo] = numpy.std(calresponseI[bolo])/numpy.sqrt(Chunks)
        outcalresponseQ[bolo] = numpy.std(calresponseQ[bolo])/numpy.sqrt(Chunks)
    frame[Output] = outcalresponse
    frame[Output+'I']=outcalresponseI
    frame[Output+'Q']=outcalresponseQ

@core.indexmod
def CalibratorResponsePerBoloSN(frame, Input='CalibratorResponsePerBolo',
                         Noise='CalibratorResponsePerBoloNoise',
                         Output='CalibratorResponsePerBoloSN'):
    '''
    Calculate the signal-to-noise of the calibrator response given the
    response and a noise estimate calculated previously.
    '''

    if frame.type != core.G3FrameType.Scan or Input not in frame:
        return

    calresponse = frame[Input]
    calresponseI = frame[Input+'I']
    calresponseQ = frame[Input+'Q']
    calnoise = frame[Noise]
    calnoiseI = frame[Noise+'I']
    calnoiseQ = frame[Noise+'Q']
    sn = core.G3MapDouble()
    snI = core.G3MapDouble()
    snQ = core.G3MapDouble()
    for bolo in calresponse.keys():
        sn[bolo] = numpy.abs(calresponse[bolo])/calnoise[bolo]
        snI[bolo] = numpy.abs(calresponseI[bolo])/calnoiseI[bolo]
        snQ[bolo] = numpy.abs(calresponseQ[bolo])/calnoiseQ[bolo]
    frame[Output] = sn
    frame[Output+'I'] = snI
    frame[Output+'Q'] = snQ

@core.pipesegment
def AnalyzeCalibratorData(pipe, Output='CalibratorResponse',
                          Input='CalTimestreams',
                          DropNoise=True, NoiseChunks=10):
    '''
    Analyze the calibrator data in the frame. The main response dictionary will
    be in Output, with the signal-to-noise estimate in OutputSN. If
    DropNoise is False, the noise estimate will stored as OutputNoise. Since
    this can be calculated trivially from the response and the signal-to-noise,
    DropNoise is True by default.
    '''

    pipe.Add(ExtractCalibratorFrequency, output = Output+'Frequency')
    pipe.Add(ExtractAzEl, output1=Output+'MedianOnlineBoresightAz', output2=Output+'MedianOnlineBoresightEl')
    #pipe.Add(ExtractCalibratorPhase, InputTimestreams=Input, Output=Output + 'Phase')
    #pipe.Add(FitCalibratorResponse, InputTimestreams=Input, Output=Output, Phase=Output + 'Phase')
    #pipe.Add(CalibratorResponseNoise, InputTimestreams=Input, Output=Output + 'Noise', Phase=Output + 'Phase', Chunks=NoiseChunks)
    #pipe.Add(CalibratorResponseSN, Input=Output, Noise=Output + 'Noise', Output=Output + 'SN')
    # For calculating the response, noise, SN, phase for each bolometer individually.
    pipe.Add(FitCalibratorResponsePerBolo, InputTimestreams=Input, Output=Output+'PerBolo', PerBoloPhase=Output+'PerBoloPhase')
    pipe.Add(CalibratorResponsePerBoloNoise, InputTimestreams=Input, Output=Output + 'PerBoloNoise', Chunks=NoiseChunks)
    pipe.Add(CalibratorResponsePerBoloSN, Input=Output+'PerBolo', Noise=Output + 'PerBoloNoise', Output=Output + 'PerBoloSN')
    #todel = [Output + 'Phase']
    #if DropNoise:
    #    todel.append(Output + 'Noise')
    #pipe.Add(core.Delete, keys=todel, type=core.G3FrameType.Scan)

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

        if frame.type == core.G3FrameType.Scan and self.cal_prefix+'PerBolo' in frame:

            # accumulate calibrator frames
            cframe = core.G3Frame(core.G3FrameType.Calibration)
            #cframe[self.cal_prefix] = frame[self.cal_prefix]
            #cframe[self.cal_prefix + 'SN'] = frame[self.cal_prefix + 'SN']
            cframe[self.cal_prefix + 'Frequency'] = frame[self.cal_prefix + 'Frequency']
            # save the per bolo response, SN, and phase.
            cframe[self.cal_prefix+'PerBolo'] = frame[self.cal_prefix+'PerBolo']
            cframe[self.cal_prefix+'PerBoloI'] = frame[self.cal_prefix+'PerBoloI']
            cframe[self.cal_prefix+'PerBoloQ'] = frame[self.cal_prefix+'PerBoloQ']
            cframe[self.cal_prefix + 'PerBoloSN'] = frame[self.cal_prefix + 'PerBoloSN']
            cframe[self.cal_prefix + 'PerBoloSNI'] = frame[self.cal_prefix + 'PerBoloSNI']
            cframe[self.cal_prefix + 'PerBoloSNQ'] = frame[self.cal_prefix + 'PerBoloSNQ']
            cframe[self.cal_prefix + 'PerBoloPhase'] = frame[self.cal_prefix +'PerBoloPhase']
            cframe[self.cal_prefix +'MedianOnlineBoresightEl']= frame[self.cal_prefix +'MedianOnlineBoresightEl']
            cframe[self.cal_prefix +'MedianOnlineBoresightAz']= frame[self.cal_prefix +'MedianOnlineBoresightAz']
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
