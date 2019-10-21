#!/usr/bin/env python
import unittest
from spt3g import core, timestreamflagging
import numpy as np

class TestFlagging(unittest.TestCase):
    def test_add_flag(self):
        frame = core.G3Frame()
        timestreamflagging.add_flag(frame, 'FlagsKey', 'flag_1', ['A'])
        timestreamflagging.add_flag(frame, 'FlagsKey', 'flag_2', ['A', 'B'])
        timestreamflagging.add_flag(frame, 'FlagsKey', 'flag_3', ['B'])
        timestreamflagging.add_flag(frame, 'FlagsKey', 'flag_3', ['A', 'B'])
        
        flags = frame['FlagsKey']
        self.assertTrue('A' in flags)
        self.assertTrue('B' in flags)
        for fs in ['flag_1', 'flag_2', 'flag_3']:
            self.assertTrue(fs in flags['A'])
        for fs in [ 'flag_2', 'flag_3']:
            self.assertTrue(fs in flags['B'])

    def test_ts_flagging(self):
        frame = core.G3Frame(core.G3FrameType.Scan)
        timestreamflagging.add_flag(frame, 'FlagsKey', 'flag_1', ['A'])
        ts_map = core.G3TimestreamMap()
        ts_map['A'] = core.G3Timestream()
        ts_map['B'] = core.G3Timestream()

        frame['in_ts'] = ts_map
        timestreamflagging.RemoveFlaggedTimestreams(frame, 'in_ts', 'FlagsKey','out_ts')
        self.assertTrue('A' in frame['in_ts'])
        self.assertFalse('A' in frame['out_ts'])
        self.assertTrue('B' in frame['in_ts'])
        self.assertTrue('B' in frame['out_ts'])
    


if __name__ == '__main__':
    unittest.main()



