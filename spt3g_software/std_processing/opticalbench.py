import numpy as np
from spt3g import core

@core.indexmod
def UpsampleBenchPosition(frame, bolots='RawTimestreams_I'):
    bolo_t = numpy.arange(frame[bolots].n_samples)/frame[bolots].sample_rate + frame[bolots].start.time
    slow_benchpos = frame.pop(['BenchPosition'], None)
    benchpos = core.G3TimestreamMap()
    bench_t = numpy.array([t.time for t in frame['BenchPosition'][0].times()])
    for key in slow_benchpos.keys():
        benchpos[key] = numpy.interp(bolo_t, bench_t, slow_benchpos[key])
    frame['BenchPosition'] = benchpos

def TrimBenchPosition(frame, trim_key = 'RawTimestreams_I', bench_key = 'BenchPosition'):
    '''
    Trim the stored bench positions so that they begin and end at the same
    time as `bolots`.  By default, trims on bolometer timestreams, but 
    could trim on pointing timestreams.
    '''
    if frame.type != core.G3FrameType.Scan:
        return
    trim_start = frame[trim_key].start.time
    trim_stop = frame[trim_key].stop.time
    benchtimes = np.array([t.time for t in frame[bench_key].times()])
    i_start = 0
    i_stop = len(benchtimes) - 1
    for i, t in enumerate(benchtimes):
        # Note that we may include one sample immediately before and
        # after the trimming timestream.  This allows for upsampling later.
        if t <= trim_start:
            i_start = i
        if t >= trim_stop:
            i_stop = i
            break
    benchpos = core.G3TimestreamMap()
    oldpos = frame.pop(bench_key, None)
    for k in oldpos.keys():
        benchpos[k] = oldpos[k][i_start:i_stop + 1]
    frame[bench_key] = benchpos
    return
