from spt3g import core
from spt3g.todfilter import dftutils

def do_noise_aux(input_file, bins, output_calibration_file = None):
    frames = [f for f in core.G3File(input_file)]
    if len(frames) > 2:
        core.log_fatal("Too many frames in the noise stare")
    labels = ['noise_stare_av_psd_%.2fHz_%.2fHz'%( b[0]/core.G3Units.hz, b[1]/core.G3Units.hz) for b in bins]
    cal_frame = analyze_noise_data(frame['CalTimestreams_I'], bins, labels)
    if output_calibration_file != None:
        writer = core.G3Writer(output_calibration_file)
        writer(cal_frame)
    return cal_frame

def analyze_noise_data(ts_map,
                       bins,
                       labels):
    frame = core.G3Frame(core.G3FrameType.Calibration)
    psd_ests, freqs = dftutils.get_dft_of_ts_map(ts_map)
    for i, b in enumerate(bins):
        output_map = core.G3MapDouble()
        for k in psd_ests.keys():
            output_map[k] = dftutils.average_psd_over_bin( psd_ests[k], freqs,
                                                           b[0], b[1])
        frame[labels[i]] = output_map
    return frame
