import numpy
from spt3g import core, gcp, std_processing, dfmux
import scipy.signal

#  31-Oct-2015:01:03:08   31-Oct-2015:01:35:30
cena_obs = [['28-Oct-2015:17:16:35',   '28-Oct-2015:17:48:55'],
            ['29-Oct-2015:11:33:33',   '29-Oct-2015:12:05:56'],
            ['31-Oct-2015:01:03:08',   '31-Oct-2015:01:35:30'],
            ['01-Nov-2015:14:26:01',   '01-Nov-2015:14:58:22'],
            ['03-Nov-2015:02:22:00',   '03-Nov-2015:02:54:21']]

#28-Oct-2015
for cenaob in cena_obs:
    pipe = core.G3Pipeline()
    pipe.Add(std_processing.ARCTimerangeReader, start_time=cenaob[0], stop_time=cenaob[1])
    pipe.Add(gcp.ARCExtract)
    pipe.Add(gcp.UnpackSPTpolHKData)
    pipe.Add(gcp.GCPMuxDataDecoder)
    pipe.Add(std_processing.BuildScanFramesFromRawData, flac_compress=True)
    pipe.Add(core.Dump)
    pipe.Add(core.G3Writer, filename='cena-%s.g3'%str(cenaob[0][:11]) )
    pipe.Run(profile=True)
