import numpy as np
import pdb
import pickle
from spt3g import core, std_processing, gcp
from sptpol_software.util.tools import struct

pipe = core.G3Pipeline()
pipe.Add(std_processing.ARCTimerangeReader, start_time=start_time, stop_time=stop_time, basedir=arcdir)
pipe.Add(gcp.ARCExtract)

d = struct()
d.linsens_avg_l1 = []
d.linsens_avg_l2 = []
d.linsens_avg_r1 = []
d.linsens_avg_r2 = []
d.scu_temp = []

def fillPointingStructure(fr):
    if fr.type != core.G3FrameType.GcpSlow:
        return
    d.linsens_avg_l1.append(fr['TrackerPointing'].linsens_avg_l1)
    d.linsens_avg_l2.append(fr['TrackerPointing'].linsens_avg_l2)
    d.linsens_avg_r1.append(fr['TrackerPointing'].linsens_avg_r1)
    d.linsens_avg_r2.append(fr['TrackerPointing'].linsens_avg_r2)
    d.scu_temp.append(fr['TrackerPointing'].scu_temp)
    
pipe.Add(fillPointingStructure)
pipe.Run()

d.linsens_avg = np.median(np.array([np.array(d.linsens_avg_l1).flatten(), 
                                    np.array(d.linsens_avg_l2).flatten(), 
                                    np.array(d.linsens_avg_r1).flatten(), 
                                    np.array(d.linsens_avg_r2).flatten()]),axis=1)
d.scu_temp = np.median(np.array(d.scu_temp),axis=0)

thermolin = np.hstack([d.scu_temp,d.linsens_avg])

