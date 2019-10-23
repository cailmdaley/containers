from spt3g import core, dfmux, calibration, std_processing
import sys

pipe = core.G3Pipeline()
pipe.Add(core.G3Reader(sys.argv[1]))
pipe.Add(dfmux.FixedLengthScans, N=100000)
pipe.Add(dfmux.DfMuxCollator)
pipe.Add(std_processing.ConsolidateHousekeepingToScanFrame)
pipe.Add(dfmux.ConvertTimestreamUnits, Input='RawTimestreams_I')
pipe.Add(core.Dump)
pipe.Add(calibration.calibrator_analysis.AnalyzeCalibratorData)
pipe.Add(core.G3Writer(sys.argv[2]))

pipe.Run()
