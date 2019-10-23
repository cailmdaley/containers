from spt3g import core, gcp
import numpy as np
from functools import reduce
import operator
import glob
import os.path

class FieldExtractor:
    def __init__(self, fields, size=0):
        self.fields = fields
        self.data = [np.zeros(size) for field in self.fields]
        self.jEntry = 0

    def __call__(self, frame):
        if frame:
            for jfield, field in enumerate(self.fields):
                datum = reduce(operator.getitem, field, frame)
                if self.jEntry >= len(self.data[jfield]):
                    self.data[jfield] = np.append(self.data[jfield], datum)
                else:
                    self.data[jfield][self.jEntry] = datum
            self.jEntry += 1                


class Downsampler:
    def __init__(self, downsample_factor):
        self.downsample_factor = downsample_factor
        self.jFrame = 0

    def __call__(self, frame):
        if frame:
            self.jFrame += 1
            if self.jFrame == self.downsample_factor:
                self.jFrame = 0
                return frame
            else:
                return []


def extract_arc_data(fields, time_start, time_stop, downsample_factor=1,
                     arcfile_path='/spt/data/arc/'):
    '''
    Function for extracting data from ARC files in a manner that sucks minimally.
    
    Parameters
    ----------
    fields : Python list of tuples
        Tuple of strings or integers that index the quantity to extract from
        the ARC file. For example, the UC stage temperature lives at:
            frame['array']['cryo']['temperature'][0][0]
        To extract this data from the ARC file, one would set
        fields = [('array', 'cryo', 'temperature', 0, 0)]
    time_start : G3Time
        Extract data starting at this time (UTC).
    time_stop : G3Time
        Stop extracting data at this time (UTC).
    downsample_factor : int, default=1
        Factor by which to downsample data
    arcfile_path : str, default='/spt/data/arc/'
        Path to ARC file data. Default is the location on amundsen/scott.

    Returns
    -------
    data : Python list of numpy arrays
        List of arrays containing the requested data.
    '''
    
    # figure out files that match requested times
    arcfilelist = np.sort(glob.glob('{}/*.dat'.format(arcfile_path)))
    filetime_list = [core.G3Time(os.path.basename(fname).rstrip('.dat'))
                     for fname in arcfilelist]
    ind_inrange = [i for i in range(len(filetime_list))
                   if filetime_list[i]>time_start and filetime_list[i]<time_stop] 
    if ind_inrange[0] != 0:
        ind_inrange.insert(0, ind_inrange[0]-1)
    fnames = arcfilelist[ind_inrange]

    field_extractor = FieldExtractor(fields, np.floor(1000*len(fnames) / downsample_factor - 1))

    # extract the data
    for fname in fnames:
        pipe = core.G3Pipeline()
        pipe.Add(gcp.ARCFileReader(fname))
        pipe.Add(Downsampler, downsample_factor=downsample_factor)
        pipe.Add(field_extractor)
        pipe.Run()

    return field_extractor.data
