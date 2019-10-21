from spt3g import core, gcp, std_processing
import numpy as np
import pylab as py
import pickle as pk
import os
import glob

def logTimeToObsID(logtime):
    '''
    Convert logtime (from GCP log) into ObsID.
    '''
    t = core.G3Time(logtime)
    x = int((t.time - core.G3Time('20170101_000000').time)/core.G3Units.s)

    return x

def ObsIDtoLogTime(obsID):
    '''
    Convert an obsID into a start time.
    '''

    t = core.G3Time('20170101_000000') + obsID*core.G3Units.s
    return t.Summary()


def grabTimestreamFromScanified(obsID_location, timestream_key=None):
    '''
    Extract timestream from bolo timestream_key for a given obsID 
    location.
    '''
    files = np.sort(glob.glob(os.path.join(obsID_location,'*.g3')))

    x = []
    ts = []

    for d in files:         
        print d
        this_file = core.G3File(d)
        if timestream_key != None:
            x += [f['RawTimestreams_I'][timestream_key] 
                  for f in this_file if f.type == core.G3FrameType.Scan]
        else:
            print 'Whoa there!  You need to specify a bolo key!'
            return

    for i in range(len(x)):
        ts += list(np.array(x[i]))

    return np.array(ts)


def grabElnodFromG3File(file_location, timestream_key=None):
    '''                                                                                                              
    Extract timestream from bolo timestream_key for a given file
    location.                                                                                                        
    '''
    x = []
    ts = []
    this_file = core.G3File(file_location)
    f = this_file.next()
    if timestream_key != None:
        x += f['ElnodEigenvalueDominantVectorI'][timestream_key]
    else:
        print 'Whoa there!  You need to specify a bolo key!'
        return

    for i in range(len(x)):
        ts += list(np.array(x[i]))

    return np.array(ts)


def grabRawAzElFromScanified(obsID_location):
    '''                                                                                                              
    Extract raw Az/El for a given obsID location.                                                                                                        
    '''
    files = np.sort(glob.glob(os.path.join(obsID_location,'*.g3')))

    az = []
    el = []

    for d in files:
        print d
        this_file = core.G3File(d)
        az += [f['Boresight'][timestream_key]
                  for f in this_file if f.type == core.G3FrameType.Scan]

        if timestream_key != None:
            x += [f['RawTimestreams_I'][timestream_key]
                  for f in this_file if f.type == core.G3FrameType.Scan]
        else:
            print 'Whoa there!  You need to specify a bolo key!'
            return

    for i in range(len(x)):
        ts += list(np.array(x[i]))

    return np.array(ts)



def grabARC(start_time, stop_time, filename, basedir='/spt_data/arc'):
    pipe = core.G3Pipeline()
    pipe.Add(std_processing.ARCTimerangeReader, start_time=start_time, 
             stop_time=stop_time, basedir=basedir)
    pipe.Add(gcp.ARCExtract)
    pipe.Add(core.G3Writer, filename=filename)
    pipe.Run(profile=False)

    return filename

