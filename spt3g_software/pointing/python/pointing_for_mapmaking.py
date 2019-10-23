from spt3g import core, coordinateutils
from . import offline_pointing as op
import numpy
from copy import copy

@core.indexmod
def NaiveBoresightPointing(frame, output='RawBoresight', bolots='RawTimestreams_I',
                           interpolate_over_dropouts = True):
    '''
    Makes two timesteams (<output>Az, <output>El) by linear interpolation
    of the GCP tracker registers onto a bolometer timestream specified by
    the bolots argument.  Interpolates over pointing dropouts by default
    '''
    def interpolate_glitches(pointing_ts, sample_rate = 100., max_spd = 20.0):
        speed = abs(numpy.diff(pointing_ts / core.G3Units.deg) * sample_rate)
        spd_glitches = numpy.where(speed > 5.0)[0]
        if len(spd_glitches) == 0:
            return pointing_ts
        out_ts = core.DoubleVector(pointing_ts)
        not_glitches = list(range(len(out_ts)))
        pos_glitches = []
        for i_glitch in range(1, len(spd_glitches)):
            if spd_glitches[i_glitch] - spd_glitches[i_glitch - 1] != 1:
                continue
            pos_glitches.append(spd_glitches[i_glitch])
            not_glitches.remove(pos_glitches[-1])
        
        interpolatable_pointing = numpy.array([pointing_ts[i] for i in not_glitches])
        x = numpy.linspace(0, 100, len(pointing_ts))
        for glitch_index in pos_glitches:
            out_ts[int(glitch_index)] = numpy.interp(x[glitch_index],
                      x[not_glitches], interpolatable_pointing)
        return out_ts
    
    if frame.type != core.G3FrameType.Scan:
        return
    pointing_t = [t.time for t in frame['TrackerStatus'].time]
    
    az_input = frame['TrackerStatus'].az_pos
    el_input = frame['TrackerStatus'].el_pos

    # Unwrap az so that straight linear interpolation does the right thing
    # and the checks in interpolate_glitches() don't decide wrapping in az
    # is a pointing glitch and "interpolate" the data into madness.
    az_diff = numpy.diff(az_input)
    az_diff[az_diff > numpy.pi*core.G3Units.rad] -= 2*numpy.pi*core.G3Units.rad
    az_diff[az_diff < -numpy.pi*core.G3Units.rad] += 2*numpy.pi*core.G3Units.rad
    az_input = numpy.concatenate(([az_input[0]], az_input[0] + numpy.cumsum(az_diff)))

    if interpolate_over_dropouts:
        # Find pointing glitches (single, non-physical samples)
        # Since we know max velocity pretty well, use that
        sample_rate = float(len(pointing_t)) / (pointing_t[-1] - pointing_t[0]) * core.G3Units.s
        az_input = interpolate_glitches(az_input, sample_rate)
        el_input = interpolate_glitches(el_input, sample_rate)

    if bolots is not None:
        bolo_t = numpy.arange(frame[bolots].n_samples)/frame[bolots].sample_rate + frame[bolots].start.time
        az_output = core.G3Timestream(numpy.interp(bolo_t, pointing_t, az_input) % (2*numpy.pi*core.G3Units.rad))
        el_output = core.G3Timestream(numpy.interp(bolo_t, pointing_t, el_input))
        az_output.start = el_output.start = frame[bolots].start
        az_output.stop = el_output.stop = frame[bolots].stop
    else:
        start = core.G3Time()
        start.time = pointing_t[0]
        stop = core.G3Time()
        stop.time = pointing_t[-1]        
        az_output = core.G3Timestream(az_input)
        el_output = core.G3Timestream(el_input)
        az_output.start = el_output.start = start
        az_output.stop = el_output.stop = stop
    
    frame[output + 'Az'] = az_output
    frame[output + 'El'] = el_output

    
@core.indexmod
def CalculateLocalOffsetPointing(frame, source, trans_key,
                                 x_offset_key='AzOffset',
                                 y_offset_key='ElOffset',
                                 ts_ref_key='OnlineBoresightAz',
                                 max_throw=5 * core.G3Units.deg):
    '''
    Calculates timestreams of offsets in local coordinates, relative the position of
    a source in celestial coordinates.

    Processing Arguments:
        source [string] :
            Name of reference source against which to calculate local offsets
            Ra and Dec positions will be extracted from keys named
            <source>Ra and <source>Dec, respectively.
        max_throw [float] : maximum valid X-offset in the output data.
        ts_ref_key -> G3Timestream : a reference timestream whose start time is used
            to determine the position of the reference source in Celestial coordinates.
    Frame Input data:
        trans_key -> G3VectorQuat: The pointing transformation that maps (1,0,0)
            to Celestial coordinates.

    Frame Output data:
        x_offset_key -> G3Timestream:
        y_offset_key -> G3Timestream:
            Timestreams of X and Y offsets of the boresight relative to the position
            of the source.
    '''

    if trans_key not in frame:
        return

    # source position at the start of the scan, in local coordinates
    from spt3g.std_processing.sourcemaps import get_source_ra_dec
    source_pos = get_source_ra_dec(source, at_time=frame[ts_ref_key].start)
    source_q = coordinateutils.ang_to_quat(source_pos[0], source_pos[1])
    trans = frame[trans_key]
    source_xy = trans**-1 * source_q * trans

    # offset from source in local coordinates
    x, y = coordinateutils.quat_to_ang(source_xy)
    x = numpy.mod(-x / core.G3Units.rad + 2 * numpy.pi, 2 * numpy.pi)
    x[x > numpy.pi] -= 2 * numpy.pi
    x[x < -numpy.pi] += 2 * numpy.pi

    if max_throw is not None:
        xthrow = numpy.nanmax(numpy.abs(x))
        xmsg = "x throw %f < %f degrees" % (xthrow / core.G3Units.deg, max_throw / core.G3Units.deg)
        assert xthrow < max_throw, xmsg

        ythrow = numpy.nanmax(numpy.abs(y)) 
        ymsg = "y throw %f < %f degrees" % (ythrow / core.G3Units.deg, max_throw / core.G3Units.deg)
        assert ythrow < max_throw, ymsg

    ts_ref = frame[ts_ref_key]

    ts_x = core.G3Timestream(x)
    ts_x.start = ts_ref.start
    ts_x.stop = ts_ref.stop
    frame[x_offset_key] = ts_x

    ts_y = core.G3Timestream(y)
    ts_y.start = ts_ref.start
    ts_y.stop = ts_ref.stop
    frame[y_offset_key] = ts_y

    return frame


@core.pipesegment
def CalculateCoordTransRotations(
        pipe,
        raw_az_key='RawBoresightAz',
        raw_el_key='RawBoresightEl',
        output='OnlineBoresight',
        transform_store_key='OnlineRaDecRotation',
        model='OnlinePointingModel',
        flags=['az_tilts', 'el_tilts', 'flexure', 'collimation', 'refraction']
        #, 'thermolin'] # Thermolin broken as of 4/13/17
):
    '''
    Calculate the coordinate transformations from raw encoder Az/El coordinates
    to sky coordinates, including the effects of the pointing model.

    Arguments
    ---------
    raw_az_key, raw_el_key : string
        The frame keys that contain the raw encoder positions of the telescope.
    output : string
        The prefix to apply to all output fields.  The following fields will be
        created by this module:  <output>Az, <output>El, <output>Ra, <output>Dec.
    transform_store_key : string
        The key where the complete pointing rotation transiformation quaternions
        will be stored.
    model : string
        The pointing model to use to compute corrected coordinates.
    flags : list of strings
        A list of the pointing corrections to apply.
    '''

    # Add corrected pointing to scans.
    pipe.Add(op.CorrectBoresightPointing,
             raw_az_key=raw_az_key, raw_el_key=raw_el_key, output=output,
             model=model, flags=flags)

    # Add astronomical pointing
    pipe.Add(coordinateutils.azel.LocalToAstronomicalPointing,
             az_timestream=output+'Az', el_timestream=output+'El',
             ra_timestream=output+'Ra', dec_timestream=output+'Dec')

    # Compute a set of coordinates offset from boresight by a small constant
    # elevation.  These coordinates will be used to construct a rotation
    # quaternion that accounts for rotations of the boresight relative to the
    # local meridian.
    offset_input_az = '_Offset' + raw_az_key
    offset_input_el = '_Offset' + raw_el_key
    offset_output = '_Offset' + output
    offset_val = 5 * core.G3Units.arcmin

    def AddOffsetPointing(frame):
        if frame.type != core.G3FrameType.Scan:
            return
        az = frame[raw_az_key]
        el = copy(frame[raw_el_key])

        # with some slop, if we are near the pole the ell offset can do some weird things
        numpy.asarray(el)[el > (numpy.pi - 2 * offset_val)] -= 2 * offset_val
        el += offset_val
        frame[offset_input_az] = az
        frame[offset_input_el] = el

    pipe.Add(AddOffsetPointing)

    # Add corrected pointing to scans.
    pipe.Add(op.CorrectBoresightPointing,
             raw_az_key=offset_input_az, raw_el_key=offset_input_el,
             output=offset_output,
             model=model, flags=flags)

    # Add astronomical pointing
    pipe.Add(coordinateutils.azel.LocalToAstronomicalPointing,
             az_timestream=offset_output+'Az', el_timestream=offset_output+'El',
             ra_timestream=offset_output+'Ra', dec_timestream=offset_output+'Dec')

    # Compute coordinate rotation
    pipe.Add(coordinateutils.coordsysmodules.FillCoordTransRotations,
             transform_store_key=transform_store_key,
             bs_az_key=raw_az_key, bs_el_key=raw_el_key,
             bs_ra_key=output+'Ra', bs_dec_key=output+'Dec',
             offset_az_key=offset_input_az, offset_el_key=offset_input_el,
             offset_ra_key=offset_output+'Ra', offset_dec_key=offset_output+'Dec',
             do_bad_transform=False)

    # Remove intermediate keys
    pipe.Add(core.Delete, type=core.G3FrameType.Scan,
             keys=[offset_input_az, offset_input_el,
                   offset_output+'Az', offset_output+'El',
                   offset_output+'Ra', offset_output+'Dec'])

