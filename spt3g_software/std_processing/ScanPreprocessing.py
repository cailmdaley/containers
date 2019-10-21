from spt3g import core, dfmux, todfilter, timestreamflagging, xtalk, calibration
from spt3g.calibration.apply_t_cal import ApplyTCalibration
from spt3g.mapmaker import CheckForCalFrames
import copy
import numpy

hz = core.G3Units.hz

@core.pipesegment
def CalibrateRawTimestreams(pipe, units=core.G3TimestreamUnits.Tcmb,
                            i_data_key='RawTimestreams_I',
                            q_data_key='RawTimestreams_Q',
                            invert_xtalk=True, allow_missing_xtalk_matrix=True,
                            output='CalTimestreams', keep_intermediate=False,
                            q_output = None, keep_original=True,
                            calibrator_iq_rotation=False
                            ):
    '''
    Convert raw timestreams into calibrated timestreams in the requested
    units, applying whatever calibration is needed to get there. Also applies
    the optimal I/Q rotation from the elnod analysis and crosstalk inversion
    if a crosstalk matrix is present in the data stream. Deletes all internally
    generated intermediate data products unless keep_intermediate is set to
    True.
    '''

    if q_output is None:
        q_rotate_key = None
    else:
        q_rotate_key = 'Q_RawTimestreamsRotated'
    if calibrator_iq_rotation:
        pipe.Add(calibration.calibrator_analysis.CalibratorRotateIQ,
          i_data_key=i_data_key, q_data_key=q_data_key,
          i_rotate_key='RawTimestreamsRotated', q_rotate_key=q_rotate_key,
          destructive=not keep_original)
    else:
        pipe.Add(calibration.elnod_analysis.RotateIQ,
          i_data_key=i_data_key, q_data_key=q_data_key,
          i_rotate_key='RawTimestreamsRotated', q_rotate_key=q_rotate_key,
          destructive=not keep_original)
    if units == core.G3TimestreamUnits.Tcmb:
        muxcal = core.G3TimestreamUnits.Power
        muxout = 'TimestreamsWatts'
    else:
        muxcal = units
        muxout = output
    if not keep_original:
        pipe.Add(core.Delete, keys=[i_data_key, q_data_key])

    if invert_xtalk:
        pipe.Add(xtalk.CrosstalkCleanedTimestreams,
          input='RawTimestreamsRotated', output=muxout, units=muxcal,
          ignore_missing=allow_missing_xtalk_matrix)
    else:
        pipe.Add(dfmux.ConvertTimestreamUnits, Input='RawTimestreamsRotated', 
          Output=muxout, Units=muxcal)
    
    if not keep_intermediate:
        pipe.Add(core.Delete, keys=['RawTimestreamsRotated'])
    if units == core.G3TimestreamUnits.Tcmb:
        pipe.Add(ApplyTCalibration, Input=muxout, Output=output)

    if not q_output is None:
        if units == core.G3TimestreamUnits.Tcmb:
            q_muxcal = core.G3TimestreamUnits.Power
            q_muxout = 'Q_TimestreamsWatts'
        else:
            q_muxcal = units
            q_muxout = q_output
        pipe.Add(dfmux.ConvertTimestreamUnits, Input='Q_RawTimestreamsRotated', 
                 Output=q_muxout, Units=q_muxcal)
        if units == core.G3TimestreamUnits.Tcmb:
            pipe.Add(ApplyTCalibration, Input=q_muxout, Output=q_output)

    if not keep_intermediate:
        delkeys = ['Q_RawTimestreamsRotated']
        if units == core.G3TimestreamUnits.Tcmb:
            delkeys.append(muxout)
            if not q_output is None:
                delkeys.append(q_muxout)
        pipe.Add(core.Delete, type=core.G3FrameType.Scan, keys=delkeys)

@core.pipesegment
def DropWasteFrames(pipe):
    '''
    Check that calibration frames arrive before scan frames
    Drop repeated wiring, observation, and calibration frames
    Drop turn-arounds
    '''
    pipe.Add(CheckForCalFrames)
    pipe.Add(core.DeduplicateMetadata)
    pipe.Add(lambda fr: not fr.get('Turnaround', False))

@core.usefulfunc
def InterpOverNans(frame, input_timestreams='RawTimestreams_I', 
                   input_q_timestreams='RawTimestreams_Q', 
                   badtimes='NanSampleTimes', max_nnan=50):
    '''
    Check for individual bolometer samples set to NaN (presumably by
    the DAQ when individual packets have been lost or are in the wrong
    order). If there are a few isolated cases, interpolate over them
    and record the timestamp of the samples interpolated over. If
    there are too many, bail with a warning. If there are none,
    proceed happily.  If the key `badtimes` already exists in a frame,
    assume this has already been done and do nothing.
    '''
    if frame.type != core.G3FrameType.Scan:
        return
    if badtimes in frame:
        return

    if input_timestreams in frame:
        tsmap = frame.pop(input_timestreams)
    else:
        return

    qthere = False
    if input_q_timestreams in frame:
        tsmap_q = frame.pop(input_q_timestreams)
        qthere = True

    nantimes = core.G3MapVectorTime()
    badbolos = []

    times = tsmap.times()
    for bolo in tsmap.keys():
        thists = numpy.asarray(tsmap[bolo])
        isf = numpy.isfinite(thists)
        if not isf.all():
            if len(thists) - isf.sum() > max_nnan:
                badbolos.append(bolo)
                continue
            nantimes[bolo] = [thistime for (thistime,thisisf) in zip(times,isf) if not thisisf]
            x1 = numpy.arange(len(thists))
            thists_interp = numpy.interp(x1,x1[isf],thists[isf])
            tsmap[bolo][:] = thists_interp
            if qthere:
                thists_q = numpy.asarray(tsmap_q[bolo])
                thists_q_interp = numpy.interp(x1,x1[isf],thists_q[isf])
                tsmap_q[bolo][:] = thists_q_interp

    if len(badbolos) > 0:
        core.log_warn(str(len(badbolos)) + " bolometers had more than " + str(max_nnan) + " NaN samples and were left alone.")

    frame[input_timestreams] = tsmap
    if qthere:
        frame[input_q_timestreams] = tsmap_q
    frame[badtimes] = nantimes
