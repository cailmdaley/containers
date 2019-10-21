import numpy
from spt3g import core

def extract_source_integral_from_sptsz_map(path, center=None):
    '''
    Returns the calibration numbers to get K_CMB hard-coded below
    by analyzing an SPT-SZ source map. This uses a (hard-coded) 4x4 arcmin
    box. Do not change this -- it *must* match (well, weakly must) the box
    size in fit_fluxandpointing.py and opacity.py.

    The optional center argument is the pixel position (y, x) around which to
    construct the box.  If none is given, the center is chosen to be at the
    hottest pixel in the map.
    '''

    import astropy.io.fits as fits

    f = fits.open(path)[0]
    if center is None:
        center = numpy.unravel_index(numpy.argmax(f.data), f.data.shape)
        # print('Found map peak at pixel {}.'.format(center))
    assert(f.header['CUNIT1'].strip() == 'deg')
    assert(f.header['TUNIT'].strip() == 'K_CMB')

    mapres = numpy.abs(f.header['CD1_1']*core.G3Units.deg)

    box = 2*core.G3Units.arcmin # +/- 2 arcmin
    chunk = f.data[center[0] - int(box/mapres):center[0] + int(box/mapres),
      center[1] - int(box/mapres):center[1] + int(box/mapres)]*(mapres**2)

    return numpy.sum(chunk)

# The following numbers are in K_cmb*(angular units)**2 and are
# integrated fluxes of the listed source[s] as extracted by the
# above function from SPT-SZ data.
kcmb_conversion_factors = {
  'RCW38': {
    90.0*core.G3Units.GHz: 4.0549662e-07*core.G3Units.K,
    150.0*core.G3Units.GHz: 2.5601153e-07*core.G3Units.K,
    220.0*core.G3Units.GHz: 2.8025804e-07*core.G3Units.K,
  },
  'MAT5A': {
    90.0*core.G3Units.GHz: 2.5738063e-07*core.G3Units.K, # center (608, 555)
    150.0*core.G3Units.GHz: 1.7319235e-07*core.G3Units.K,
    220.0*core.G3Units.GHz: 2.145164e-07*core.G3Units.K,
  },
}

@core.indexmod
class ApplyTCalibration(object):
    '''
    Apply most recent sky calibration data (or elevation-interpolated
    calibration) and calibrator response to convert timestreams in watts
    into timestreams in sky-calibrated units (K_cmb).
    '''
    def __init__(self, Calibrator='CalibratorResponse',
      Source=['RCW38', 'MAT5A'], Input='TimestreamsWatts',
      Output='CalTimestreams', BoresightEl='RawBoresightEl',
      OpacityCorrection=True, InKCMB=True, interp=False):
        '''
        Converts the timestreams in the map Input (in watts) into a new
        timestream map named Output in sky-calibrated units.

        Calibrator and SkyCal are the names of the keys of the appropriate
        calibration quantities in the calibration frame and can usually be
        left at their defaults.
        '''
        self.input = Input
        self.output = Output
        self.interp = interp
        self.calkey_interp = Calibrator + 'Interp'
        self.calkey_el = Calibrator + 'El'
        self.calkey = Calibrator
        self.caldata = None
        self.caldata_el = None
        self.source = Source
        self.skydata = None
        self.bp = None
        self.skytransmission = None
        self.elkey = BoresightEl

        self.kcmb = InKCMB
        self.kcmbconversions = None

        self.opacitycorrection = OpacityCorrection

    def __call__(self, frame):
        if frame.type == core.G3FrameType.Calibration:
            if isinstance(self.source, list):
                source = None
                for src in self.source:
                    if src + 'FluxCalibration' in frame:
                        source = src
                        break
            else:
                source = self.source
            if self.interp:
                if self.calkey_interp in frame:
                    self.caldata = frame[self.calkey_interp]
                    self.caldata_el = frame[self.calkey_el]
            else:
                if self.calkey in frame:
                    self.caldata = frame[self.calkey]
            if source is not None and source + 'FluxCalibration' in frame:
                self.skydata = frame[source + 'FluxCalibration']
            if 'BolometerProperties' in frame:
                self.bp = frame['BolometerProperties']
            if self.opacitycorrection and source is not None and source + 'SkyTransmission' in frame:
                self.skytransmission = {}
                for bolo,props in self.bp:
                    if props.band <= 0 or not numpy.isfinite(props.band):
                        continue
                    strband = str(int(props.band/core.G3Units.GHz))
                    if strband in frame[source + 'SkyTransmission']:
                        self.skytransmission[bolo] = \
                          frame[source + 'SkyTransmission'][strband]
                    else:
                        # XXX Warn about this?
                        self.skytransmission[bolo] = 1.
            if self.kcmb and source is not None and source + 'IntegralFlux' in frame:
                self.kcmbconversions = {}
                for bolo,fl in frame[source + 'IntegralFlux'].iteritems():
                    if bolo not in self.bp or self.bp[bolo].band <= 0 or \
                      not numpy.isfinite(self.bp[bolo].band):
                        continue
                    band = self.bp[bolo].band
                    if band not in kcmb_conversion_factors[source]:
                        core.log_error('Detector %s has an observing band %f '
                          'with no K_cmb calibration' %
                        (bolo, band/core.G3Units.GHz), unit='ApplyTCalibration')
                    self.kcmbconversions[bolo] = \
                      fl/kcmb_conversion_factors[source][band]

        elif frame.type == core.G3FrameType.Scan:
            if self.skydata is None:
                core.log_fatal('Calibration not found for sources {}'.format(self.source))
            assert(self.caldata is not None)
            assert(self.skydata is not None)

            if self.interp:
                el = frame[self.elkey].mean()

            out = core.G3TimestreamMap()
            inmap = frame[self.input]
            for bolo,in_ts in inmap.iteritems():
                assert(in_ts.units == core.G3TimestreamUnits.Power)
                out_ts = core.G3Timestream(in_ts)
                if bolo not in self.caldata:
                    continue
                if self.interp:
                    cal = numpy.interp(el, self.caldata_el, self.caldata[bolo])
                else:
                    cal = self.caldata[bolo]
                if bolo not in self.skydata:
                    continue
                cal *= self.skydata[bolo]
                if self.opacitycorrection:
                    if bolo not in self.skytransmission:
                        continue
                    cal *= self.skytransmission[bolo]
                if self.kcmb:
                    if bolo not in self.kcmbconversions:
                        continue
                    cal *= self.kcmbconversions[bolo]
                if cal != 0:
                    out_ts[:] /= cal
                else:
                    out_ts[:] = numpy.nan
                out_ts.units = core.G3TimestreamUnits.Tcmb
                out[bolo] = out_ts

            frame[self.output] = out
