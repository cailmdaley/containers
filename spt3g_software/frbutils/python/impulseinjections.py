from __future__ import print_function
import pickle, random, copy
from spt3g import core
import numpy as np
from spt3g.util import genericutils as GU
from spt3g.calibration import template_groups

class FastTransientSignalBeam(object):
    #because I was jamming in sptpol data res is in units of arcminutes,
    #the output has the appropriate use of G3Units.
    def __init__(self, 
                 actual_map, res, plot = True):
        from scipy.interpolate import interp2d

        smooth= False
        downsample = False

        if smooth:
            import scipy.signal#.convolve2d
            smoothing_size = 5
            t = 1 - np.abs(np.linspace(-1, 1, smoothing_size))
            kernel = t.reshape(smoothing_size, 1) * t.reshape(1, smoothing_size)
            kernel /= kernel.sum()   # kernel should sum to 1!  :) 
            actual_map = scipy.signal.convolve2d(actual_map, kernel, mode='same')
        #actual_map = ven.map[390:410, 390:410]
        if downsample:
            ds_fac = 2
            actual_map = GU.downsample_2d_array(actual_map, ds_fac = ds_fac)
            res *= ds_fac
        actual_map /= np.sum(actual_map) * res**2.0 * (1.0/60.0)**2 * ( np.pi / 180.0 )**2.0 #gives a 1 jy source with pixels in units of Jy/sr
        actual_map /= 396.3e6 # divide by (Jy/K/sr) conversion to make map in units of K EARLIER UNITS
        #actual_map /= 396.3e3 # divide by (Jy/K/sr) conversion to make map in units of K
        print(np.max(actual_map))
        if plot:
            import pylab as pl
            pl.clf()
            pl.imshow(actual_map, interpolation = 'nearest', cmap = pl.cm.gray)
            pl.title("Venus")
            pl.show()
        self.actual_map = actual_map
        self.n_pixels = np.size(self.actual_map)
        self.sky_coverage = self.n_pixels * res**2.0 * core.G3Units.arcmin**2.0
        self.n_x = actual_map.shape[1]
        self.n_y = actual_map.shape[0]
        self.f = interp2d(np.arange(self.n_x), np.arange(self.n_y), self.actual_map)
        
    def get_signal_k(self, signal_jy):
        x_place = np.array(random.random() * (self.n_x - 1))
        y_place = np.array(random.random() * (self.n_y - 1))
        return signal_jy * self.f(x_place, y_place)[0]

    def apply_signal_two_deg_s(self, signal_jy, signal  ):
        input_stream_sample_rate = 195312.5 #hz
        scan_speed = 120.0 * core.G3Units.arcmin #arcmin/sec
        sample_step = scan_speed / input_stream_sample_rate 

        res = (self.sky_coverage/self.n_pixels)**0.5

        y_place = np.array(random.random() * (self.n_y - 1))
        x_place = np.array(random.random() * (self.n_x - 1))

        s_len = len(signal)
        center_sample = np.where(signal ==np.max(signal))[0][0]

        scan_samples = x_place + sample_step/res * (np.arange(0, s_len) - center_sample)
        #import pdb; pdb.set_trace()
        #print('slen', s_len)
        return signal_jy * signal * (self.f(scan_samples, y_place)).flatten()


    def get_sky_coverage(self):
        return self.sky_coverage




def get_fake_frb_template(width, decay_fraction = 0.2):
    assert(0)
    a = np.zeros(2 * width+1)
    a[width] = 1
    for i in range(1,width):
        sgn = (-1)**(i)
        a[width-i] = abs(decay_fraction) * sgn / float(i)
        a[width+i] = abs(decay_fraction) * sgn / float(i) 
    return a


def apply_test_beam_to_in_stream(in_stream, use_scan_correction = False, x=None, y=None):
    source_model = lambda x, y: np.exp( -(x)**2.0/(2 * 0.44**2)) * np.exp( -(y)**2.0/(2 * 0.44**2)) #arcmin
    input_stream_sample_rate = 195312.5 #hz
    scan_speed = 120.0 #arcmin/sec
    sample_step = scan_speed / input_stream_sample_rate

    sample_scale = 1.0
    if x is None or y is None:
        x = (np.random.rand()-0.5)*2.0 * sample_scale
        y = (np.random.rand()-0.5)*2.0 * sample_scale
    if not use_scan_correction:
        return source_model(x,y) * in_stream
    else:
        xs = np.arange(len(in_stream)) * sample_step
        xs = xs + (x-xs[len(xs)//2])
        return source_model(xs,y) * in_stream
        

def create_input_stream( time_scale, sample_delay, curve_type, fluence):
    '''
    time_scale is in units of seconds
    
    sample_delay is in seconds

    curve_type 0 flat, 1 sin^2, 2: triangle, 3: gaussian

    fluence is in Jy ms
    '''
    #the timescale of 1 sample
    sample_rate = 195312.5 # Hz
    sample_timescale = 1.0/sample_rate * 1e3 # ms
    sf = (fluence/sample_timescale) 
    if curve_type == 0:
        ugh = np.zeros(int((sample_delay + time_scale) * sample_rate))
        offset = int(sample_delay * sample_rate)
        ugh[offset:] = 1
        ugh *= sf  / np.sum(ugh)
        return ugh
    elif curve_type == 1:
        ugh = np.zeros(int((sample_delay + time_scale * 2.0) * sample_rate))
        offset = int(sample_delay * sample_rate)
        npix = len(ugh[ offset : ] )
        ugh[ offset : ]  = np.sin(np.arange(npix) * (np.pi / (npix-1.0)) )**2.0
        ugh *= sf  / np.sum(ugh)
    elif curve_type == 2:
        #triangle
        ugh = np.zeros(int((sample_delay + time_scale * 2.0) * sample_rate))
        offset = int(sample_delay * sample_rate)
        npix = len(ugh[ offset : ] )
        ugh[ offset : ]  = (npix-1)/2.0 - np.abs(np.arange(npix) - (npix-1)/2.0)
        ugh *= sf  / np.sum(ugh)
    elif curve_type == 3:
        #gaussian
        ugh = np.zeros(int((sample_delay + time_scale * 3.0) * sample_rate))
        offset = int(sample_delay * sample_rate)
        npix = len(ugh[ offset : ] )
        ugh[ offset : ]  = np.exp(-1 * (np.arange(npix) - (npix-1)/2.0)**2.0 / (2 * ((npix-1) / 2.355 / 3.0)**2))
        ugh *= sf  / np.sum(ugh)
    else:
        raise RuntimeError("I DON'T KNOW WHAT YOU WANT")
    return ugh

def dfmux_impulse_response(input_stream):
    '''
    normalized such that if input_stream = np.zeros( infinity) + 1
       the output will be a stream of 1s
    '''
    cic2 = reduce(np.convolve, [np.ones(512) for x in range(8)])
    cic2resp = np.convolve( np.array([0] + list(cic2)),  input_stream)[::512] 
    fir = [-1, 0, 1, 0, -2, -2, 3, 5, -2, -9, -2, 12, 8, -14, -18, 11, 31, -2, -44, -16, 53, 44, -53, -81, 36, 122, 4, -157, -69, 174, 158, -157, -262, 94, 366, 27, -443, -208, 465, 436, -400, -686, 221, 916, 85, -1071, -511, 1089, 1027, -909, -1571, 486, 2054, 203, -2366, -1140, 2380, 2264, -1974, -3462, 1037, 4565, 516, -5351, -2740, 5532, 5664, -4708, -9321, 2210, 13798, 3477, -19197, -17433, 23236, 67107, 67107, 23236, -17433, -19197, 3477, 13798, 2210, -9321, -4708, 5664, 5532, -2740, -5351, 516, 4565, 1037, -3462, -1974, 2264, 2380, -1140, -2366, 203, 2054, 486, -1571, -909, 1027, 1089, -511, -1071, 85, 916, 221, -686, -400, 436, 465, -208, -443, 27, 366, 94, -262, -157, 158, 174, -69, -157, 4, 122, 36, -81, -53, 44, 53, -16, -44, -2, 31, 11, -18, -14, 8, 12, -2, -9, -2, 5, 3, -2, -2, 0, 1, 0, -1]
    ir = np.convolve(cic2resp, fir)[::2] / 6.189605749097244e+26
    #normalized response is: 6.189605749097244e+26
    return ir


def get_signal_response(  time_scale, sample_delay, curve_type, fluence, fts = None, 
                          include_scan = False ):
    if not include_scan:
        if fts is None:
            sig = dfmux_impulse_response( create_input_stream( time_scale, sample_delay, curve_type, fluence ))
        else:
            sig = fts.get_signal_k(1.0) * dfmux_impulse_response( create_input_stream( time_scale, sample_delay, curve_type, fluence ))
    else:
        assert(not fts is None)
        sig = dfmux_impulse_response(fts.apply_signal_two_deg_s(1.0, create_input_stream( time_scale, sample_delay, curve_type, fluence )))
    return sig, np.where( sig == np.max(sig))[0][0]


def add_fake_signal_to_timestream(timestream, timestream_2, 
                                  time_scale, curve_type, fluence, fts, inject_index,
                                  two_deg_s_speed = False):
    sample_delay = random.random() *  2.0 * 0.00524288

    sig, max_ind = get_signal_response(  time_scale, sample_delay, curve_type, fluence, fts, include_scan = two_deg_s_speed)

    if inject_index is None:
        inject_index = random.randint(0, len(timestream) - 1)
    start_index = inject_index - max_ind
    stop_index = len(sig) + start_index
    if start_index < 0:
        clip_low = abs(start_index)
        start_index = 0
    else:
        clip_low = 0

    if stop_index >= len(timestream):
        clip_high = len(sig) - (stop_index - len(timestream))
        stop_index = len(timestream) 
    else:
        clip_high = len(sig)
    timestream[start_index:stop_index] += sig[clip_low:clip_high]
    if not timestream_2 is None:
        timestream_2[start_index:stop_index] += sig[clip_low:clip_high]

    return inject_index, sig[max_ind]


@core.scan_func_cache_data(bolo_props = 'BolometerProperties')
def InjectFrbSignal(frame, ts_key, out_ts_key,
                    time_scale, curve_type, fluence, fts, 
                    inject_info_dets_key = 'InjectedDets', 
                    inject_info_ind_key = 'InjectedIndex',
                    inject_info_amp_key = 'InjectedAmp',
                    inject_index = None,
                    bolo_props = None,
                    two_deg_s_speed = False
                    ):
    if frame.type != core.G3FrameType.Scan:
        return

    pixel_to_bolo_map = template_groups.get_template_groups(
        bolo_props, per_band = False, per_pixel = True, include_keys = True)

    partner_map = {}
    for k in pixel_to_bolo_map.keys():
        if len(pixel_to_bolo_map[k]) > 2:
            ks = filter( lambda c: bolo_props[c].physical_name[-1] == 'X' or bolo_props[c].physical_name[-1] == 'Y', pixel_to_bolo_map[k])
        else:
            ks = pixel_to_bolo_map[k][0]
        if len(ks) < 2:
            continue
        partner_map[ks[0]] = ks[1]
        partner_map[ks[1]] = ks[0]

    ts_map = copy.copy(frame[ts_key])
    ts_lst = ts_map.keys()
    
    lnumber = 0
    keep_going = True #so we only get events where both detectors aren't flagged in case we haven't already done that
    while (keep_going):
        lnumber += 1
        id_x = ts_lst[random.randint(1,len(ts_lst))-1]

        if id_x in partner_map:
            id_y = partner_map[id_x]
        else:
            continue

        assert(lnumber < 10000)

        if (id_y in ts_lst):
            keep_going = False

    #print('\n'*3, 'Injecting', id_y, id_x, '\n'*3)
    scan_index, max_amp = add_fake_signal_to_timestream(ts_map[id_x], ts_map[id_y],
                                                        time_scale, curve_type, fluence, 
                                                        fts, inject_index, two_deg_s_speed = two_deg_s_speed)
    #print('injected at', scan_index, id_x, id_y)
    
    frame[out_ts_key] = ts_map
    frame[inject_info_ind_key] = int(scan_index)
    frame[inject_info_dets_key] = core.G3VectorString([id_x, id_y])
    frame[inject_info_amp_key] = max_amp

def LookForInjectedEvents(frame,
                          frb_event_key, slop = 2,
                          inject_info_ind_key = 'InjectedIndex',
                          inject_info_dets_key = 'InjectedDets'):
    if frame.type != core.G3FrameType.Scan:
        return
    if inject_info_dets_key not in frame or inject_info_ind_key not in frame:
        return
    if frb_event_key not in frame:
        return
    frb_events = copy.copy(frame[frb_event_key])
    inj_dets = frame[inject_info_dets_key]
    inj_ind = frame[inject_info_ind_key]
    
    for i, ev in enumerate(frb_events):
        if (ev.scan_index-slop <= inj_ind and inj_ind <= ev.scan_index+slop):
            ev_dets = [di.bid for di in ev.det_info]
            for inj_det in inj_dets:
                try:
                    ind = ev_dets.index(inj_det)
                except:
                    continue
                dinf = ev.det_info[ind]
                #print('setting was inj')
                dinf.was_injected = True
                ev.det_info[ind] = dinf
                #print('f', dinf.was_injected)
            frb_events[i] = ev
            #print('l',frb_events[i].det_info[0].was_injected)
    del frame[frb_event_key]
    frame[frb_event_key] = frb_events



if __name__ == '__main__':
    #AA filter plots
    if 0:
        import pylab as pl
        in_sig = create_input_stream( time_scale=0.001, sample_delay=0, curve_type=3, fluence=1)
        pl.clf()
        pl.plot(in_sig)
        pl.xlabel("Samples (Input Sample Rate)")
        pl.ylabel("Signal")
        pl.title("Digital Anti-aliasing Filter Input Waveform")
        pl.savefig('/home/nlharr/Thesis/FrbPlots/AntiAliasingInputWaveform.png')
        pl.clf()
        out_sigs = []
        for sd in np.arange(0.00524288*0, 0.00524288*1, 0.00524288/10.0):
            out_sigs.append(dfmux_impulse_response(create_input_stream( time_scale=0.001, sample_delay=sd, curve_type=3, fluence=1)))
        for out_s in out_sigs:
            pl.plot(out_s)

        pl.xlabel("Samples (Output Sample Rate)")
        pl.ylabel("Signal With Input Offsets")
        pl.title("Digital Anti-aliasing Filter Output Waveforms")
        pl.xlim([20,60])
        pl.savefig('/home/nlharr/Thesis/FrbPlots/AntiAliasingOutputWaveforms.png')

        pl.clf()
        out_maxes = []
        for sd in np.arange(0.00524288*0, 0.00524288*2, 0.00524288/4000.0):
            out_maxes.append(np.max(dfmux_impulse_response(create_input_stream( time_scale=0.001, sample_delay=sd, curve_type=3, fluence=1))))
        pl.hist(out_maxes, bins = 40)
        pl.title("Histogram of Maximum of AA Filter Output")
        pl.xlabel("Maximum of AA Filter Output Amplitude For Various Input Time Offsets")
        pl.savefig('/home/nlharr/Thesis/FrbPlots/AntiAliasingOutputAmplitudeHistogram.png')
    #curve shape plots
    if 0:
        import pylab as pl
        t_widths = np.arange(1e-3, 40e-3, 0.00524288/10.0)
        pl.clf()

        for ct in [0,1,2,3]:
            out_maxes = []
            for tw in t_widths: 
                local_maxes = []
                for sd in np.arange(0.00524288*0, 0.00524288*1, 0.00524288/30.0):
                    local_maxes.append(np.max(dfmux_impulse_response(create_input_stream( time_scale=tw, sample_delay=sd, curve_type=ct, fluence=1))))
                out_maxes.append(max(local_maxes))
            pl.loglog(t_widths, out_maxes)
        pl.title("Anti-aliasing Filter Response For Various Waveforms")
        pl.xlabel("Input Signal FWHM (s)")
        pl.ylabel("Anti-aliasing Filter Response Maximum")
        pl.ion()
        pl.legend(['Square', 'Sin', 'Triangle', 'Gaussian'])
        pl.savefig('/home/nlharr/Thesis/FrbPlots/AntiAliasingOutputAmplitudeCurveVariations.png')
    if 0:
        if 1:
            map_fn = '/home/nlharr/frb_side_products/venus_beam_150_2013_noXtalk.pkl'
            ven = pickle.load(open(map_fn))
            res = ven.reso_arcmin
            actual_map = ven.map[385:415, 385:415]
            print(res)
        else:
            #map_fn = '/home/nlharr/spt_code/spt3g_software/frbutils/python/beams/0537-441_map_20120923_065033_150ghz.hdf5'
            map_fn = '/home/nlharr/spt_code/spt3g_software/frbutils/python/beams/0537-441_map_20130621_181809_150ghz.hdf5'
            import sptpol_software.observation.sky
            m = sptpol_software.observation.sky.read(map_fn)
            m = m.removeWeight()['T']
            res = m.reso_arcmin
            # 317 317 center
            actual_map = m.map[303:333, 303:333]
            1/0
        fts = FastTransientSignalBeam(actual_map, res, plot = False)
        if 1:
            pickle.dump(fts, open('/spt/public/nlharr/frb_side_products/fast_transient_signal_conversion.pkl', 'w'))
        if 0:
            while 1:
                print(fts.get_signal_k(1))
    #scan plots
    if 0:
        import pylab as pl
        n_sims = 20000
        ts_ms = 1
        for ts_ms in [1,2,4,8]:
            print("doing ts_ms", ts_ms)
            pl.clf()
            time_scale = ts_ms * 0.001        
            maxes = []
            for i in range(n_sims):
                if (i%1000 == 0):
                    print(i)
                in_stream = create_input_stream( time_scale=time_scale, sample_delay=np.random.rand()*0.00524288, curve_type=3, fluence=1)
                o_stream = dfmux_impulse_response(apply_test_beam_to_in_stream(in_stream, False))
                maxes.append(np.max(o_stream))
            bins = pl.hist(maxes, alpha=0.7, bins = 25)
            o_maxes = []
            for i in range(n_sims):
                in_stream = create_input_stream( time_scale=time_scale, sample_delay=np.random.rand()*0.00524288, curve_type=3, fluence=1)
                o_stream = dfmux_impulse_response(apply_test_beam_to_in_stream(in_stream, True))
                o_maxes.append(np.max(o_stream))
            pl.hist(o_maxes, bins=bins[1], alpha=0.7)
            pl.xlabel('Maxima of AA Filtered Output')
            pl.title('Histogram of Monte Carlo Simulations of Output Amplitude, %dms FWHM'%ts_ms)
            pl.legend(["One Point Sample Beam", "Scanned Beam"])
            pl.savefig('/home/nlharr/Thesis/FrbPlots/AAOutputBeamScanCorrectionFineBins%dms.png'%ts_ms)

    if 0:
        fns = ['fast_transient_signal_conversion_ven.pkl', 
               'fast_transient_signal_conversion_quas.pkl']
        sigies = []
        for fn in fns:
            fts = pickle.load(open(fn))
            sigy = []
            for i in xrange(int(6e5)):
                sigy.append(fts.get_signal_k(1))
            sigies.append(sigy)

    #plot the decimation filter impulse response
    if 0:
        import matplotlib.pyplot as plt
        #overplots the impulses
        sr = 25e6/2**17.
        sds = np.arange( 0, 1.05/192., 0.1/sr)

        mv = 0
        for s in sds:
            sr, offset = get_signal_response(time_scale= 0.00001, sample_delay = s, curve_type = 3, fluence = 1)
            maxval = np.max(sr)
            mv = maxval if maxval > mv else mv
        print(mv)
        for s in sds:
            sr, offset = get_signal_response(time_scale= 0.00001, sample_delay = s, curve_type = 3, fluence = 1)
            plt.plot(sr[offset-40:]/mv)
            
        plt.xlim([30,50])
        plt.ylabel("Decimated Impulse Response (Normalized)")
        plt.xlabel("Sample (190.7 Hz)")
        plt.title("Decimation Filter Impulse Responses")
        plt.show()
            
    #plot the point spread function
    if 0:
        import matplotlib.pyplot as plt
        fts = pickle.load(open('/home/nlharr/frb_side_products/fast_transient_signal_conversion.pkl'))
        plt.imshow(fts.actual_map, interpolation = 'nearest', cmap='gray')
        plt.title("Detector Point Spread Function Estimated from Venus")
        plt.xlabel("Pixel (0.1 arcmin)")
        plt.ylabel("Pixel (0.1 arcmin)")
        plt.show()

    #plot the distribution of peak heights
    if 0:
        import matplotlib.pyplot as plt

        mvals = []
        for i in range(3000):
            sample_delay = random.random() *  0.00524288 * 2
            sig, __ = get_signal_response(  time_scale = 0.002, sample_delay = sample_delay, 
                                        curve_type = 3, fluence = 1, fts = None )

            mval = np.max(sig)
            mvals.append(mval)
        plt.hist(mvals)
        plt.show()
            
