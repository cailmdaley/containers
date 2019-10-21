import os,pickle
import numpy as np
from spt3g import core, std_processing

def find_nearest_scan(pydfmux_date='',scan_type='calibrator'):
    datadir='/spt/data/bolodata/downsampled'
    scan_types=os.listdir(datadir)
    if scan_type not in scan_types:
        print scan_types
        raise TypeError('Requested Scan Type does not exist')
    scan_dir=os.path.join(datadir,scan_type)
    scans=np.sort(np.int32(np.array(os.listdir(scan_dir))))
    pydfmux_obsdate=std_processing.time_to_obsid(pydfmux_date)
    ddate=scans-pydfmux_obsdate
    
    nearest_scan=np.where(ddate >0)[0]
    if len(nearest_scan) ==0:
        raise ValueError('No Data Found After Given Date')
    nearest_scan=scans[np.min(nearest_scan)]
    return nearest_scan
    
