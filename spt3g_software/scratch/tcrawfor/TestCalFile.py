# tools to check calfile-derived calibration information against
# original, hard-coded values.

import numpy
from spt3g import core, gcp, std_processing, dfmux
from spt3g.gcp import CalFile, ARCExtractor, ARCExtractor_tcmod

class FakeFrame(dict):
    def __init__(self, input_type=core.G3FrameType.GcpSlow):
        self.type = input_type
        self['array'] = {}
        self['antenna0'] = {}

class FakeFrameValueHolder(object):
    def __init__(self, input_value):
        self.value = input_value

def make_fake_frame():
    '''
    Make a fake frame with fake antenna0 and array data to check
    calibration process.

    DEPRECATED, DO NOT USE
    '''

    f = FakeFrame()

    mjd0 = 57754.0 # Jan 1, 2017
    s_to_mjd = 1./365./86400.
    time0 = core.G3Time('20170101_000000')
    nfast = 100
    az0 = 82.6 # arbitrary
    el0 = -57.5 # middle of 500d field
    azrate = 2.0 # deg/s
    elrate = 0.01 # just so I can check cal, should really be 0 up to az tilt
    temp0 = -15.6 # C
    encoffs = numpy.array([13.2,-5.6]) # mas
    tilts = numpy.array([153623.,-93960.]) # mas
    tilts_avg = tilts
    linsens_avg = numpy.array([43.5,672.,26.,123.])
    t_air = -72.4 # C
    p_air = 823.6 # mb

    f['antenna0']['frame'] = {}
    f['antenna0']['frame']['utc'] = time0
    f['antenna0']['acu'] = {}
    f['antenna0']['acu']['az_pos'] = FakeFrameValueHolder(az0)
    f['antenna0']['acu']['el_pos'] = FakeFrameValueHolder(az0)
    f['antenna0']['acu']['az_err'] = FakeFrameValueHolder(0.01) # arbitrary
    f['antenna0']['acu']['el_err'] = FakeFrameValueHolder(0.001) # arbitrary
    f['antenna0']['acu']['az_rate'] = FakeFrameValueHolder(azrate)
    f['antenna0']['acu']['el_rate'] = FakeFrameValueHolder(elrate)
    f['antenna0']['acu']['px_checksum_error_count'] = FakeFrameValueHolder(0)
    f['antenna0']['acu']['px_resync_count'] = FakeFrameValueHolder(0)
    f['antenna0']['acu']['px_resync_timeout_count'] = FakeFrameValueHolder(0)
    f['antenna0']['acu']['px_resyncing'] = FakeFrameValueHolder(0)
    f['antenna0']['acu']['px_timeout_count'] = FakeFrameValueHolder(0)
    f['antenna0']['acu']['restart_count'] = FakeFrameValueHolder(0)
    f['antenna0']['acu']['state'] = FakeFrameValueHolder(0)
    f['antenna0']['acu']['acu_status'] = FakeFrameValueHolder(0)
    f['antenna0']['tracker'] = {}
    f['antenna0']['tracker']['source'] = FakeFrameValueHolder('ra0hdec-57.5')
    time_fast = numpy.arange(nfast)/numpy.float(nfast)
    f['antenna0']['tracker']['utc'] = time_fast + mjd0
    f['antenna0']['tracker']['actual'] = numpy.zeros([2,nfast])
    azpos_cal = time_fast*azrate + az0
    f['antenna0']['tracker']['actual'][0] = azpos_cal/57.2958
    f['antenna0']['tracker']['actual'][1] = el0/57.2958
    f['antenna0']['tracker']['actual_rates'] = numpy.zeros([2,nfast])
    f['antenna0']['tracker']['actual_rates'][0] = azrate/(1.0725/3.6e6)
    f['antenna0']['tracker']['actual_rates'][1] = elrate/(0.4625/3.6e6)
    f['antenna0']['tracker']['expected'] = numpy.zeros([2,nfast])
    f['antenna0']['tracker']['expected'][0] = azpos_cal/57.2958
    f['antenna0']['tracker']['expected'][1] = el0/57.2958
    f['antenna0']['tracker']['expected_rates'] = numpy.zeros([2,nfast])
    f['antenna0']['tracker']['expected_rates'][0] = azrate/(1.0725/3.6e6)
    f['antenna0']['tracker']['expected_rates'][1] = elrate/(0.4625/3.6e6)
    f['antenna0']['tracker']['state'] = numpy.repeat("",nfast)
    f['antenna0']['tracker']['inControl'] = numpy.repeat(1,nfast)
    f['antenna0']['tracker']['acu_seq'] = numpy.repeat(-1,nfast)
    f['antenna0']['tracker']['scan_flag'] = numpy.repeat(1,nfast)
    f['array']['frame'] = {}
    f['array']['frame']['features'] = FakeFrameValueHolder(3)
    f['antenna0']['tracker']['utc'] = numpy.repeat(time0,nfast)
    f['antenna0']['scu'] = {}
    f['antenna0']['scu']['temp'] = numpy.repeat(temp0,nfast)
    f['antenna0']['tracker']['encoder_off'] = numpy.array([encoffs[0],encoffs[1]])
    f['antenna0']['tracker']['tilt_xy_avg'] = numpy.zeros([2,nfast])
    f['antenna0']['tracker']['tilt_xy_avg'][0] = tilts_avg[0]/(3600/3.6e6)
    f['antenna0']['tracker']['tilt_xy_avg'][1] = tilts_avg[1]/(3600/3.6e6)
    f['antenna0']['tracker']['horiz_mount'] = numpy.zeros([2,nfast])
    f['antenna0']['tracker']['horiz_mount'][0] = \
        f['antenna0']['tracker']['actual'][0][0]
    f['antenna0']['tracker']['horiz_mount'][1] = \
        f['antenna0']['tracker']['actual'][1][0]
    f['antenna0']['tracker']['horiz_off'] = numpy.zeros([2,nfast])
    f['antenna0']['tracker']['horiz_off'][0] = \
        numpy.asarray(f['antenna0']['tracker']['actual'][0]) - \
        f['antenna0']['tracker']['actual'][0][0]
    f['antenna0']['tracker']['horiz_off'][1] = \
        numpy.asarray(f['antenna0']['tracker']['actual'][1]) - \
        f['antenna0']['tracker']['actual'][1][0]
    f['antenna0']['tracker']['linear_sensor_avg'] = numpy.zeros([4,100])
    for i in numpy.arange(4):
        f['antenna0']['tracker']['linear_sensor_avg'][i] = linsens_avg[i]
    f['array']['weather'] = {}
    f['array']['weather']['airTemperature'] = FakeFrameValueHolder(t_air)
    f['array']['weather']['pressure'] = FakeFrameValueHolder(p_air)
    
    return f

def grab_boards(f, boardlist = []):
    keylist = ['TrackerPointing',
               'ACUStatus',
               'SourceName',
               'array',
               'TrackerStatus',
               'GCPFeatureBits',
               'antenna0']
    try:
        for key in keylist:
            boardlist.append(f[key])
    except:
        pass

#def test_calfile_against_arcex():
#    '''
#    Test the application of calibration params read from cal file
#    against values hard-coded in the original ARCExtract.py.
#    '''
#
#    start_time='22-Jan-2017:20:06:06'
#    stop_time='22-Jan-2017:20:06:06'
#    arcdir='/buffer/arc/'
##    fntemp1 = '/home/tcrawfor/temp_eraseme_dump_1.g3'
##    fntemp2 = '/home/tcrawfor/temp_eraseme_dump_2.g3'
#
#    data1 = []
#    data2 = []
#
#    pipe1 = core.G3Pipeline()
#    pipe1.Add(std_processing.ARCTimerangeReader, start_time=start_time, stop_time=stop_time, basedir=arcdir)
#    pipe1.Add(gcp.ARCExtract)
#    pipe1.Add(grab_boards, boardlist=data1)
##    pipe1.Add(core.G3Writer, filename=fntemp1)
#    pipe1.Run()
#
#    pipe2 = core.G3Pipeline()
#    pipe2.Add(std_processing.ARCTimerangeReader, start_time=start_time, stop_time=stop_time, basedir=arcdir)
#    pipe2.Add(ARCExtractor_tcmod.ARCExtract_tcmod)
#    pipe2.Add(grab_boards, boardlist=data2)
##    pipe2.Add(core.G3Writer, filename=fntemp2)
#    pipe2.Run()
#
##    for frame1 in core.G3File(fntemp1):
###        print frame1['array']['frame']['utc']
###        print numpy.min(frame1['antenna0']['tracker']['actual'][0])
###        print numpy.max(frame1['antenna0']['tracker']['actual'][0])
###        print frame1
##        pass
##
##    for frame2 in core.G3File(fntemp2):
###        print frame2['array']['frame']['utc']
###        print numpy.min(frame2['antenna0']['tracker']['actual'][0])
###        print numpy.max(frame2['antenna0']['tracker']['actual'][0])
###        print frame2
##        pass
#
##    return frame1, frame2
##    return data1, data2
#
#    acu_attrs =  ['acu_status',
#                  'az_err',
#                  'az_pos',
#                  'az_rate',
#                  'el_err',
#                  'el_pos',
#                  'el_rate',
#                  'px_checksum_error_count',
#                  'px_resync_count',
#                  'px_resync_timeout_count',
#                  'px_resyncing',
#                  'px_timeout_count',
#                  'restart_count',
#                  'state',
#                  'time']
#    whacu1 = (numpy.where(numpy.asarray([isinstance(field, gcp.ACUStatus) for field in data1])))[0][0]
#    acu1 = data1[whacu1]
#    whacu2 = (numpy.where(numpy.asarray([isinstance(field, gcp.ACUStatus) for field in data2])))[0][0]
#    acu2 = data2[whacu2]
#    acu_check = []
#    ftol = 1e-4
#    for attr in acu_attrs:
#        attr1 = getattr(acu1,attr)
#        if type(attr1) is float:
#            if attr1 == 0.:
#                acu_check.append(getattr(acu2,attr) == attr1)
#            else:
#                acu_check.append(numpy.abs((getattr(acu2,attr) - attr1)/attr1) < ftol)
#        else:
#            acu_check.append(getattr(acu2,attr) == attr1)
#
#    ts_attrs = ['acu_seq',
#                'az_command',
#                'az_pos',
#                'az_rate',
#                'az_rate_command',
#                'el_command',
#                'el_pos',
#                'el_rate',
#                'el_rate_command',
#                'in_control',
#                'scan_flag',
#                'state',
#                'time']
#
#    whts1 = (numpy.where(numpy.asarray([isinstance(field, gcp.TrackerStatus) for field in data1])))[0][0]
#    ts1 = data1[whts1]
#    whts2 = (numpy.where(numpy.asarray([isinstance(field, gcp.TrackerStatus) for field in data2])))[0][0]
#    ts2 = data2[whts2]
#    ts_check = []
#    for attr in ts_attrs:
#        attr1 = getattr(ts1,attr)
#        attr2 = getattr(ts2,attr)
#        stemp = numpy.size(attr1)
#        if stemp == 1:
#            if type(attr1) is float:
#                if attr1 == 0.:
#                    ts_check.append(attr2 == attr1)
#                else:
#                    ts_check.append(numpy.abs((attr2 - attr1)/attr1) < ftol)
#            else:
#                ts_check.append(attr2 == attr1)
#        else:
#            tempbool = True
#            shtemp = numpy.shape(attr1)
#            if len(shtemp) == 1:
#                for i in numpy.arange(stemp):
#                    if type(attr1[i]) is float:
#                        if attr1[i] == 0.:
#                            if attr2[i] != attr1[i]:
#                                tempbool = False
#                        else:
#                            if numpy.abs((attr2[i] - attr1[i])/attr1[i]) > ftol:
#                                tempbool = False
#                    else:
#                        if attr2[i] != attr1[i]:
#                            tempbool = False
#            else:
#                for i in numpy.arange(shtemp[0]):
#                    for j in numpy.arange(shtemp[1]):
#                        if type(attr1[i,j]) is float:
#                            if attr1[i,j] == 0.:
#                                if attr2[i,j] != attr1[i,j]:
#                                    tempbool = False
#                                else:
#                                    if numpy.abs((attr2[i,j] - attr1[i,j])/attr1[i,j]) > ftol:
#                                        tempbool = False
#                            else:
#                                if attr2[i,j] != attr1[i,j]:
#                                    tempbool = False
#            ts_check.append(tempbool)
#
#    tp_attrs = ['encoder_off_x',
#                'encoder_off_y',
#                'features',
#                'horiz_mount_x',
#                'horiz_mount_y',
#                'horiz_off_x',
#                'horiz_off_y',
#                'linsens_avg_l1',
#                'linsens_avg_l2',
#                'linsens_avg_r1',
#                'linsens_avg_r2',
#                'scu_temp',
#                'telescope_pressure',
#                'telescope_temp',
#                'tilts_x',
#                'tilts_y',
#                'time']
#
#    whtp1 = (numpy.where(numpy.asarray([isinstance(field, gcp.TrackerPointing) for field in data1])))[0][0]
#    tp1 = data1[whtp1]
#    whtp2 = (numpy.where(numpy.asarray([isinstance(field, gcp.TrackerPointing) for field in data2])))[0][0]
#    tp2 = data2[whtp2]
#    tp_check = []
#    for attr in tp_attrs:
#        attr1 = getattr(tp1,attr)
#        attr2 = getattr(tp2,attr)
#        stemp = numpy.size(attr1)
#        if stemp == 1:
#            if type(attr1) is float:
#                if attr1 == 0.:
#                    tp_check.append(attr2 == attr1)
#                else:
#                    tp_check.append(numpy.abs((attr2 - attr1)/attr1) < ftol)
#            else:
#                tp_check.append(attr2 == attr1)
#        else:
#            tempbool = True
#            shtemp = numpy.shape(attr1)
#            if len(shtemp) == 1:
#                for i in numpy.arange(stemp):
#                    if type(attr1[i]) is float:
#                        if attr1[i] == 0.:
#                            if attr2[i] != attr1[i]:
#                                tempbool = False
#                        else:
#                            if numpy.abs((attr2[i] - attr1[i])/attr1[i]) > ftol:
#                                tempbool = False
#                    else:
#                        if attr2[i] != attr1[i]:
#                            tempbool = False
#            else:
#                for i in numpy.arange(shtemp[0]):
#                    for j in numpy.arange(shtemp[1]):
#                        if type(attr1[i,j]) is float:
#                            if attr1[i,j] == 0.:
#                                if attr2[i,j] != attr1[i,j]:
#                                    tempbool = False
#                                else:
#                                    if numpy.abs((attr2[i,j] - attr1[i,j])/attr1[i,j]) > ftol:
#                                        tempbool = False
#                            else:
#                                if attr2[i,j] != attr1[i,j]:
#                                    tempbool = False
#            tp_check.append(tempbool)
#        print notavariable


def test_calfile_against_arcex():
    '''
    Test the application of calibration params read from cal file
    against values hard-coded in the original ARCExtract.py.
    '''

    start_time='22-Jan-2017:20:06:06'
    stop_time='22-Jan-2017:20:06:06'
#    start_time='28-Jan-2017:00:20:00'
#    stop_time='28-Jan-2017:00:20:00'
    arcdir='/buffer/arc/'
#    fntemp1 = '/home/tcrawfor/temp_eraseme_dump_1.g3'
#    fntemp2 = '/home/tcrawfor/temp_eraseme_dump_2.g3'

    data1 = []
    data2 = []

    pipe1 = core.G3Pipeline()
    pipe1.Add(std_processing.ARCTimerangeReader, start_time=start_time, stop_time=stop_time, basedir=arcdir)
    pipe1.Add(gcp.ARCExtract)
    pipe1.Add(grab_boards, boardlist=data1)
#    pipe1.Add(core.G3Writer, filename=fntemp1)
    pipe1.Run()

    pipe2 = core.G3Pipeline()
    pipe2.Add(std_processing.ARCTimerangeReader, start_time=start_time, stop_time=stop_time, basedir=arcdir)
    pipe2.Add(ARCExtractor_tcmod.ARCExtract_tcmod)
    pipe2.Add(grab_boards, boardlist=data2)
#    pipe2.Add(core.G3Writer, filename=fntemp2)
    pipe2.Run()

#    for frame1 in core.G3File(fntemp1):
##        print frame1['array']['frame']['utc']
##        print numpy.min(frame1['antenna0']['tracker']['actual'][0])
##        print numpy.max(frame1['antenna0']['tracker']['actual'][0])
##        print frame1
#        pass
#
#    for frame2 in core.G3File(fntemp2):
##        print frame2['array']['frame']['utc']
##        print numpy.min(frame2['antenna0']['tracker']['actual'][0])
##        print numpy.max(frame2['antenna0']['tracker']['actual'][0])
##        print frame2
#        pass

#    return frame1, frame2
#    return data1, data2

    acu_attrs =  ['acu_status',
                  'az_err',
                  'az_pos',
                  'az_rate',
                  'el_err',
                  'el_pos',
                  'el_rate',
                  'px_checksum_error_count',
                  'px_resync_count',
                  'px_resync_timeout_count',
                  'px_resyncing',
                  'px_timeout_count',
                  'restart_count',
                  'state',
                  'time']
    whacu1 = (numpy.where(numpy.asarray([isinstance(field, gcp.ACUStatus) for field in data1])))[0][0]
    acu1 = data1[whacu1]
    whacu2 = (numpy.where(numpy.asarray([isinstance(field, gcp.ACUStatus) for field in data2])))[0][0]
    acu2 = data2[whacu2]
    acu_check = []
    for attr in acu_attrs:
        attr1 = getattr(acu1,attr)
        attr2 = getattr(acu2,attr)
        if type(attr1) is core.G3Time:
            acu_check.append(getattr(acu2,attr) == attr1)
        elif type(attr1) is core.G3TimeVector:
            acu_check.append(numpy.all(numpy.asarray(getattr(acu2,attr)) == numpy.asarray(attr1)))
        else:
            acu_check.append(numpy.all(numpy.isclose(attr1,attr2)))

    ts_attrs = ['acu_seq',
                'az_command',
                'az_pos',
                'az_rate',
                'az_rate_command',
                'el_command',
                'el_pos',
                'el_rate',
                'el_rate_command',
                'in_control',
                'scan_flag',
                'state',
                'time']

    whts1 = (numpy.where(numpy.asarray([isinstance(field, gcp.TrackerStatus) for field in data1])))[0][0]
    ts1 = data1[whts1]
    whts2 = (numpy.where(numpy.asarray([isinstance(field, gcp.TrackerStatus) for field in data2])))[0][0]
    ts2 = data2[whts2]
    ts_check = []
    for attr in ts_attrs:
        attr1 = getattr(ts1,attr)
        attr2 = getattr(ts2,attr)
        if type(attr1) is core.G3Time:
            ts_check.append(getattr(ts2,attr) == attr1)
        elif type(attr1) is core.G3TimeVector:
            ts_check.append(numpy.all(numpy.asarray(getattr(ts2,attr)) == numpy.asarray(attr1)))
        else:
            ts_check.append(numpy.all(numpy.isclose(attr1,attr2)))

    tp_attrs = ['encoder_off_x',
                'encoder_off_y',
                'features',
                'horiz_mount_x',
                'horiz_mount_y',
                'horiz_off_x',
                'horiz_off_y',
                'linsens_avg_l1',
                'linsens_avg_l2',
                'linsens_avg_r1',
                'linsens_avg_r2',
                'scu_temp',
                'telescope_pressure',
                'telescope_temp',
                'tilts_x',
                'tilts_y',
                'time']

    whtp1 = (numpy.where(numpy.asarray([isinstance(field, gcp.TrackerPointing) for field in data1])))[0][0]
    tp1 = data1[whtp1]
    whtp2 = (numpy.where(numpy.asarray([isinstance(field, gcp.TrackerPointing) for field in data2])))[0][0]
    tp2 = data2[whtp2]
    tp_check = []
    for attr in tp_attrs:
        attr1 = getattr(tp1,attr)
        attr2 = getattr(tp2,attr)
        if type(attr1) is core.G3Time:
            tp_check.append(getattr(tp2,attr) == attr1)
        elif type(attr1) is core.G3TimeVector:
            tp_check.append(numpy.all(numpy.asarray(getattr(tp2,attr)) == numpy.asarray(attr1)))
        else:
            tp_check.append(numpy.all(numpy.isclose(attr1,attr2)))

    for attr, tbool in zip(acu_attrs, acu_check):
        print attr, tbool
        if tbool == False:
            print('Check failed for attribute ACUStatus.' + attr + '.\n')
            print(' Old value: ' + str(getattr(acu1,attr)) + '\n')
            print(' New value: ' + str(getattr(acu2,attr)) + '\n')

    for attr, tbool in zip(ts_attrs, ts_check):
        print attr, tbool
        if tbool == False:
            print('Check failed for attribute TSStatus.' + attr + '.\n')
            print(' Old value: ' + str((getattr(ts1,attr))[0]) + '\n')
            print(' New value: ' + str((getattr(ts2,attr))[0]) + '\n')

    for attr, tbool in zip(tp_attrs, tp_check):
        print attr, tbool
        if tbool == False:
            print('Check failed for attribute TrackerPointing.' + attr + '.\n')
            print(' Old value: ' + str((getattr(tp1,attr))[0]) + '\n')
            print(' New value: ' + str((getattr(tp2,attr))[0]) + '\n')

