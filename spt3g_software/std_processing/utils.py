from spt3g import core, dfmux

class _GetBolosHelper(object):
    def __init__(self):
        self.bolos = []
    def __call__(self, frame):
        if frame.type == core.G3FrameType.Wiring:
            self.bolos = frame['WiringMap'].keys()
            return []

@core.usefulfunc
def get_bolos_in_obs(fn):
    gb = _GetBolosHelper()
    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader,filename = fn)
    pipe.Add(gb)
    pipe.Run()
    return gb.bolos

epoch_3g = core.G3Time('01-Jan-2017:00:00:00')

@core.usefulfunc
def obsid_to_g3time(obsid):
    """
    Convert an observation ID into a G3Time instance

    Arguments
    ---------
    obsid : int
        Observation ID number

    Returns
    -------
    time : G3Time instance
        The (approximate) time corresponding to the beginning of the observation
    """
    return epoch_3g + int(obsid) * core.G3Units.seconds

@core.usefulfunc
def time_to_obsid(*time_in):
    """
    Convert an input time into an observation ID (in seconds since the 3G epoch)

    Arguments
    ---------
    time : G3Time instance, string, ctime or IRIG B code.
        This function accepts the same arguments that one can use to initialize
        a G3Time instance

    Returns
    -------
    obsID : int
        Seconds since the SPT3G epoch
    """

    time_in = core.G3Time(*time_in)
    return int((time_in.time - epoch_3g.time) / core.G3Units.seconds)
