import numpy as np
from spt3g import core

class AnalyzeNoise(object):
    '''
    Class for computing NEP, NEI, and NET. This is designed for running on
    noise stares, but it could be run on other types of observations in 
    principle.
    '''
    def __init__(self,
                 freq_ranges=[(1*core.G3Units.Hz, 2*core.G3Units.Hz),
                              (3*core.G3Units.Hz, 5*core.G3Units.Hz),
                              (10*core.G3Units.Hz, 15*core.G3Units.Hz),
                              (30*core.G3Units.Hz, 40*core.G3Units.Hz)],
                 noise_types=('NEP', 'NEI', 'NET'),
                 nep_input='',
                 nei_input='',
                 net_input='',
                 min_cal_sn=None):
        self.freq_ranges = freq_ranges
        self.noise_types = noise_types
        self.inputs = {'NEP': nep_input,
                       'NEI': nei_input,
                       'NET': net_input}
        self.min_cal_sn = min_cal_sn

        # There is a difference of a factor of sqrt(2) in the units of
        # K / sqrt(Hz) and K * sqrt(sec) because of the fact that 1/2 sec of
        # integration time corresponds to 1 Hz of bandwidth due to the 
        # Nyquist theorem. The units package in spt3g_software does not observe
        # this convention and treats K / sqrt(Hz) == K * sqrt(sec). This factor
        # corrects for this convention so that you can blindly apply units
        # from spt3g_software in the usual way without thinking about it.
        self.nyquist_factor = {'NEP': 1,
                               'NEI': 1,
                               'NET': 1./np.sqrt(2)}

        # Unit conversions in spt3g_software always transform timestreams into
        # RMS currents and powers. For power (and by extension, Kcmb) this makes
        # sense. When we refer to NEP, we are referring to the noise-equivalent
        # power *from the sky*. Since the thermal reseponse of the TES
        # is only sensitive to RMS fluctuations in electrical power, it makes
        # sense to use RMS. When we report NEI, we are typically referring the
        # noise to the SQUID input coil. The bandwidth of the SQUID input coil
        # is very high, so it is sensitive to all frequencies of interest, so
        # it makes sense to not convert to RMS. 
        self.rms_factor = {'NEP': 1,
                           'NEI': np.sqrt(2),
                           'NET': 1}

        self.normalization_str = {'NEP': 'rms',
                                  'NEI': 'amplitude',
                                  'NET': 'rms'}

    def __call__(self, frame):
        if frame.type == core.G3FrameType.Calibration and \
           "CalibratorResponseSN" in frame:
            self.cal_sn = frame["CalibratorResponseSN"]

        if frame.type == core.G3FrameType.Scan:
            from spt3g.todfilter.dftutils import get_psd_of_ts_map
            for noise_type in self.noise_types:
                frame['{}_normalization_convention'.format(noise_type)] = \
                    self.normalization_str[noise_type]
                if self.inputs[noise_type] in frame.keys():
                    # compute the PSD for this input type
                    psds, freqs = get_psd_of_ts_map(frame[self.inputs[noise_type]], pad=False)

                    for frange in self.freq_ranges:
                        low_f = frange[0]
                        high_f = frange[1]

                        # apply calibrator S/N cut for NET calculation only,
                        # if set in constructor
                        if noise_type == 'NET' and \
                           self.min_cal_sn is not None:
                            bolos_to_save = [bolo for bolo in frame[self.inputs[noise_type]].keys() \
                                             if bolo in self.cal_sn and \
                                             self.cal_sn[bolo] > self.min_cal_sn]
                            if 'MinCalibratorSNforNET' not in frame.keys():
                                frame['MinCalibratorSNforNET'] = self.min_cal_sn
                        else:
                            bolos_to_save = frame[self.inputs[noise_type]].keys()

                        # integrates noise power between low_f and high_f
                        net = core.G3MapDouble()
                        low_idx = np.where(freqs>low_f)[0][0]
                        hi_idx = np.where(freqs<high_f)[0][-1]
                        bandwidth = (high_f - low_f)

                        for bolo in bolos_to_save:
                            noise_power = np.trapz(np.asarray(psds[bolo])[low_idx:hi_idx],
                                                   x = freqs[low_idx:hi_idx])
                            net[bolo] = self.rms_factor[noise_type] * \
                                self.nyquist_factor[noise_type] / \
                                np.sqrt(bandwidth / noise_power)
                        frame['{}_{:.1f}Hz_to_{:.1f}Hz'.format(noise_type,
                                                               low_f / core.G3Units.Hz,
                                                               high_f / core.G3Units.Hz)] = net
                        

class MakeNoiseFrame(object):
    def __init__(self, noise_prefixes=('NEP', 'NEI', 'NET')):
        self.noise_prefixes = noise_prefixes
        
    def __call__(self, frame):
        self.nframe = core.G3Frame()
        for field in frame.keys():
            for prefix in self.noise_prefixes:
                if prefix in field:
                    self.nframe[field] = frame[field]
        
        if frame.type == core.G3FrameType.EndProcessing:
            return [self.nframe, frame]
        if frame.type == core.G3FrameType.Scan:
            return self.nframe
        else:
            return False
