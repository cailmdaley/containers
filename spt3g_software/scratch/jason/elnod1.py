'''
Play with an average sptpol elnod
'''

import numpy
import ipdb
from spt3g import core, gcp, std_processing, dfmux
import scipy.signal

def examine(frame):
    ipdb.set_trace()

def printtime(frame):
    if frame.type == core.G3FrameType.GcpSlow:
        print(frame['receiver']['frame']['utc'].GetFileFormatString())

pipe = core.G3Pipeline()
pipe.Add(std_processing.ARCTimerangeReader, start_time='06-Jun-2013:18:09:53', stop_time='06-Jun-2013:18:10:36')
pipe.Add(gcp.ARCExtract)
pipe.Add(gcp.UnpackSPTpolHKData)
pipe.Add(gcp.GCPMuxDataDecoder)
#pipe.Add(std_processing.BuildScanFramesFromRawData, flac_compress=True)
pipe.Add(core.Dump)
#pipe.Add(printtime)
#pipe.Add(core.G3Writer, filename='rcw38-21Apr2015.g3')

pipe.Run(profile=True)

'''
KeyError                                  Traceback (most recent call last)
/home/jason/spt3g_software/scratch/jason/elnod1.py in <module>()
     25 #pipe.Add(core.G3Writer, filename='rcw38-21Apr2015.g3')
     26 
---> 27 pipe.Run(profile=True)
     28 

/home/jason/spt3g_software/build/spt3g/core/modconstruct.pyc in Process(self, fr)
    116     class PyCallObjModule(G3Module):
    117         def Process(self, fr):
--> 118             return pycallable(fr)
    119 
    120     return PyCallObjModule()

/home/jason/spt3g_software/build/spt3g/gcp/ARCHKExtractor.py in __call__(self, f)
     57                     modhk.nuller_gain = hkregs['nullerGain'][i][j*2 + k]
     58                     modhk.demod_gain = hkregs['demodGain'][i][j*2 + k]
---> 59                     modhk.carrier_railed = hkregs['mezz' + mezz + 'Car' + mod + 'Overload'][i]
     60                     modhk.demod_railed = hkregs['mezz' + mezz + 'Adc' + mod + 'Overload'][i]
     61                     modhk.nuller_railed = hkregs['mezz' + mezz + 'Nul' + mod + 'Overload'][i]

KeyError: 'Invalid key'
'''

