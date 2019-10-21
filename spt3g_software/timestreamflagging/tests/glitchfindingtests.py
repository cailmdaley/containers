#!/usr/bin/env python

import unittest

from spt3g import core
from spt3g.timestreamflagging import glitchfinding as GF
import numpy as np

class TestGlitchFinding(unittest.TestCase):
    def test_find_simple(self):
        ts = core.G3Timestream(np.zeros(100))
        ts[70] = 100
        for kernel_size in range(4,11):
            self.assertEqual( GF.get_num_glitches(ts, thresholds = [5], 
                                                  kernel_size = kernel_size)[0],  
                              1)

        ts[30] = 100
        for kernel_size in range(4,11):
            self.assertEqual( GF.get_num_glitches(ts, thresholds = [5], 
                                                  kernel_size = kernel_size)[0],  
                              2)

    def test_add_glitches(self):
        frame = core.G3Frame(core.G3FrameType.Scan)
        ts = core.G3Timestream(np.zeros(100))
        ts[70] = 100
        ts_map = core.G3TimestreamMap()
        ts_map['A'] = ts
        
        frame['ts'] = ts_map
        thresholds = [5,500]
        GF.AddNumGlitches(frame, thresholds,  'ts')
        for i in range(2):
            self.assertTrue(frame['GlitchesThresholds'][i] == thresholds[i])
        self.assertTrue(frame['GlitchesNumberOf']['A'][0] == 1)
        self.assertTrue(frame['GlitchesNumberOf']['A'][1] == 0)

    def test_glitch_flagging(self):
        frame = core.G3Frame(core.G3FrameType.Scan)
        ts_map = core.G3TimestreamMap()

        ts = core.G3Timestream(np.zeros(100))
        ts[70] = 100
        ts_map['A'] = ts

        ts = core.G3Timestream(np.zeros(100))
        ts[70] = 100
        ts[30] = 100
        ts_map['B'] = ts
        
        frame['ts'] = ts_map
        thresholds = [5,500]

        GF.AddNumGlitches(frame, thresholds,  'ts')
        
        GF.FlagGlitches(frame, {5: 2, 500:1})
 
        self.assertTrue( frame['Flags'].keys()[0] == 'B')



if __name__ == '__main__':
    unittest.main()



