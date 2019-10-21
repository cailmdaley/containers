#!/usr/bin/env python

import unittest
import time, numpy
from spt3g.core import G3Timestream, G3Units, G3Time

from spt3g import util 

class TestPSD(unittest.TestCase):
    def test_psd(self):
        sample_rate = 250*G3Units.Hz
        frequency = 25*G3Units.Hz

        samples = numpy.arange(100000)
        tsdata = numpy.sin(2*numpy.pi*frequency/sample_rate*samples)
        ts = G3Timestream(tsdata)

        # Test with implicit sample rate
        ts.start = G3Time(0 * G3Units.s)
        ts.stop = G3Time(len(samples)/sample_rate)
        psd, freq = util.timestream_psd(ts)
        f_max = freq[numpy.where(psd==max(psd))][0]
        self.assertAlmostEqual(f_max/frequency, 1.0, places=4)

if __name__ == '__main__':
    unittest.main()
