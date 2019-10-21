import numpy, sys
from spt3g import core, gcp, std_processing

# Usage: grab_gcp_data.py start-time stop-time arc-directory output.g3
# Reads pointing data from ARC files for tilt measurements.

pipe = core.G3Pipeline()
pipe.Add(std_processing.ARCTimerangeReader, start_time=sys.argv[1], stop_time=sys.argv[2], basedir=sys.argv[3])
pipe.Add(gcp.ARCExtract)
pipe.Add(core.Dump)
pipe.Add(core.G3Writer, filename=sys.argv[4])

pipe.Run(profile=True)
