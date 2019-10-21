from spt3g import core
import time

pipe = core.G3Pipeline()

sample_rate = 152.
i = 0


def DummyFrames(frame, n):
    """
    Adds a bunch of dfmux, encoder, and irig timepoint frames
    For now, just store zeroes in all the entries
    """
    global i
    outframes = []
    while i<n:
        # dfmux frame followed by a whwp encoder frame
        dfmux = core.G3Frame(core.G3FrameType.Timepoint)
        dfmux['DfMux'] = 0 # just something to store there for now
        outframes.append(dfmux)
        
        encoder = core.G3Frame(core.G3FrameType.Timepoint)
        encoder['isDfmux'] = core.G3Bool(False)
        encoder['whwp_encoder_clk_cnts'] = 0
        encoder['whwp_encoder_cnts'] = 0
        outframes.append(encoder)

        # IRIG frames once every second
        if i % 152 == 0:
            irig = core.G3Frame(core.G3FrameType.Timepoint)
            irig['isDfMux'] = core.G3Bool(True)
            irig['whwp_irig_rising_edge_time'] = 0
            irig['whwp_irig_info'] = 0
            irig['whwp_irig_synch_clk_cnts'] = 0
            outframes.append(irig)
            
        time.sleep(1./sample_rate)
        i += 1
        
        return outframes
        
    return [] # after we're done, end by returning an empty list to the next module
pipe.Add(DummyFrames, n=10000)


pipe.Add(core.G3Writer, filename='dummy_hwp_data.g3')
pipe.Run(profile=True)
