import numpy, os
from spt3g import core
from spt3g.gcp import ACUStatus, ACUState, TrackerStatus, TrackerState, CalFile

@core.indexmod
def CalibrateFrame(f):
    '''Apply gain / offset / units from GCP cal file'''
    
    if f.type != core.G3FrameType.GcpSlow:
        return

    cf = CalFile.CalFileReader()
    cd = cf.readCalFile()
    if 'array' in f.keys():
        for map in f['array'].keys():
            for reg in f['array'][map].keys():
                try:
                    

@core.indexmod
def UnpackACUData(f):
    '''Extracts ACU status information to ACUStatus key in frame'''

    if f.type != core.G3FrameType.GcpSlow:
        return

    a = ACUStatus()
    a.time = f['antenna0']['frame']['utc']
    a.az_pos = f['antenna0']['acu']['az_pos'].value
    a.el_pos = f['antenna0']['acu']['el_pos'].value
    a.az_err = f['antenna0']['acu']['az_err'].value
    a.el_err = f['antenna0']['acu']['el_err'].value

    a.az_rate = f['antenna0']['acu']['az_rate'].value
    a.el_rate = f['antenna0']['acu']['el_rate'].value

    # 'new_*' registers not actually filled by GCP; ignore them

    a.px_checksum_error_count = f['antenna0']['acu']['px_checksum_error_count'].value
    a.px_resync_count = f['antenna0']['acu']['px_resync_count'].value
    a.px_resync_timeout_count = f['antenna0']['acu']['px_resync_timeout_count'].value
    a.px_resyncing = f['antenna0']['acu']['px_resyncing'].value
    a.px_timeout_count = f['antenna0']['acu']['px_timeout_count'].value
    a.restart_count = f['antenna0']['acu']['restart_count'].value

    a.state = ACUState(f['antenna0']['acu']['state'].value)
    a.acu_status = f['antenna0']['acu']['acu_status'].value

    f['ACUStatus'] = a

@core.indexmod
def UnpackTrackerData(f, rewrite_source_current=True):
    '''
    Extracts tracker status information to frame. If rewrite_source_current
    is True (the default), will try to rewrite source names of "current"
    if DecryptFeatureBit() has been run and either "elnod", "calibrator",
    or "noise" is present in the feature bit list to that value.
    '''

    if f.type != core.G3FrameType.GcpSlow:
        return

    # Grab the GCP source name. If it is "current", fill in something more
    # helpful from the feature bits if possible.
    source = f['antenna0']['tracker']['source'].value
    if rewrite_source_current and 'GCPFeatureBits' in f and source == 'current':
        if 'elnod' in f['GCPFeatureBits']:
            source = 'elnod'
        if 'calibrator' in f['GCPFeatureBits']:
            source = 'calibrator'
        if 'noise' in f['GCPFeatureBits']:
            source = 'noise'
    f['SourceName'] = source

    t = TrackerStatus()
    # List comprehensions are due to funny business with G3VectorFrameObject
    t.time = [tm for tm in f['antenna0']['tracker']['utc'][0]]

    # Measured values
    t.az_pos = numpy.asarray(f['antenna0']['tracker']['actual'][0])
    t.el_pos = numpy.asarray(f['antenna0']['tracker']['actual'][1])
    t.az_rate = numpy.asarray(f['antenna0']['tracker']['actual_rates'][0])
    t.el_rate = numpy.asarray(f['antenna0']['tracker']['actual_rates'][1])
    
    # Expected values
    t.az_command = numpy.asarray(f['antenna0']['tracker']['expected'][0])
    t.el_command = numpy.asarray(f['antenna0']['tracker']['expected'][1])
    t.az_rate_command = numpy.asarray(f['antenna0']['tracker']['expected_rates'][0])
    t.el_rate_command = numpy.asarray(f['antenna0']['tracker']['expected_rates'][1])

    # Status params
    if isinstance(f['antenna0']['tracker']['state'][0], core.G3String):
        # If state is all zero (LACKING), for example due to an ACU glitch,
        # the ARC reader may decide that the 8-bit array field is a string.
        # Treat it as one.
        t.state = [TrackerState(0) for s in f['antenna0']['tracker']['inControl'][0]]
    else:
        t.state = [TrackerState(s) for s in f['antenna0']['tracker']['state'][0]]
    t.acu_seq = f['antenna0']['tracker']['acu_seq'][0]
    t.in_control = core.BoolVector(f['antenna0']['tracker']['inControl'][0])
    t.scan_flag = core.BoolVector(f['antenna0']['tracker']['scan_flag'][0])

    f['TrackerStatus'] = t

@core.indexmod
def DecryptFeatureBit(f):
    '''
    Unpacks the GCP feature flags
    '''

    if f.type != core.G3FrameType.GcpSlow:
        return

    flag_array = core.G3VectorString()
    feature_bit = f['array']['frame']['features'].value

    flags = ['analyze', 'source_scan', 'cabin_shutter', 'elnod', 'pol_cal',
             'calibrator', 'every_pixel_on_src', 'skydip', 'optical', 'noise',
             'trail', 'el_scan']

    for i in enumerate(flags):
        if feature_bit & (1 << i[0]):
            flag_array.append(i[1])

    f['GCPFeatureBits'] = flag_array

@core.pipesegment
def ARCExtract(pipe):
    '''Extract GCP registers into native objects'''
    pipe.Add(CalibrateFrame)
    pipe.Add(UnpackACUData)
    pipe.Add(DecryptFeatureBit)
    pipe.Add(UnpackTrackerData)

# Need tool for tilt meter next
