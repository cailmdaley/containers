#!/usr/bin/env python 

from spt3g import core
from spt3g import todfilter

import unittest
import numpy as np

class TestPsdNorm(unittest.TestCase):
    def test_cpp_psd_parsevals(self):
        ts_len = 101
        m = core.G3TimestreamMap()
        m['a'] = core.G3Timestream(np.random.normal(size = ts_len))
        m['a'].start = core.G3Time(0)
        m['a'].stop = core.G3Time(1)

        window_func = core.G3VectorDouble( np.zeros(ts_len) + 1.0)
        psd_data = todfilter.get_psd_pybinding_grr( m, ts_len, window_func, 1.0)
        self.assertTrue( np.sum(  m['a'] * m['a'] ) - np.sum(  psd_data['a'] ) < 1e-5)

        psd_data = todfilter.get_psd_pybinding_grr( m, ts_len+40, window_func, 1.0)
        self.assertTrue( np.sum(  m['a'] * m['a'] ) - np.sum(  psd_data['a'] ) < 1e-5)

        ts_len = 102
        m['a'] = core.G3Timestream(np.random.normal(size = ts_len))
        m['a'].start = core.G3Time(0)
        m['a'].stop = core.G3Time(1)

        window_func = core.G3VectorDouble( np.zeros(ts_len) + 1.0)
        psd_data = todfilter.get_psd_pybinding_grr( m, ts_len, window_func, 1.0)
        self.assertTrue( abs(np.sum(  m['a'] * m['a'] ) - np.sum(  psd_data['a'] )) < 1e-5)

        psd_data = todfilter.get_psd_pybinding_grr( m, ts_len+40, window_func, 1.0)
        self.assertTrue( abs(np.sum(  m['a'] * m['a'] ) - np.sum(  psd_data['a'] )) < 1e-5)

        todfilter.dftutils.get_psd_of_ts_map(m)


if __name__ == '__main__':
    unittest.main()
