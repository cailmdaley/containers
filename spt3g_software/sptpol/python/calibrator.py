from spt3g import core

@core.indexmod
def DiscretizeCalibratorTimestream(frame, Input='Calibrator', Output='CalibratorOn', Level=4e5):
    '''
    Process SPTpol-style calibrator data (ADC values from a hacked-up DfMux
    board) into 3G-style data (1 if calibrator sync on, 0 if off).
    Output timestream will be 1 if the ADC value is above Level and 0 otherwise.
    '''
    if frame.type != core.G3FrameType.Scan:
        return
    if Input not in frame:
        return

    ts = core.G3Timestream(frame[Input])
    ts[ts[:] <= Level] = 0.
    ts[ts[:] > Level] = 1.
    frame[Output] = ts

