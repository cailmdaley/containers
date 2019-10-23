from spt3g import core
from spt3g.dfmux import DfMuxWiringMap, DfMuxChannelMapping


def just_make_wiring_map(frame):
    hwm = DfMuxWiringMap()
    for i in range(1729):
        mapping = DfMuxChannelMapping()
        mapping.board_serial = -1
        mapping.board_slot = -1
        mapping.crate_serial = -1
        mapping.board_ip = int(i/48)
        mapping.module = int(i/12) - 1
        mapping.channel = int(i) - 1
        hwm['channel_%s'%i] = mapping
    wiring_frame = core.G3Frame(core.G3FrameType.Wiring)
    wiring_frame['WiringMap'] = hwm
    return wiring_frame

def add_obs_header(frame):
    obs_frame = core.G3Frame(core.G3FrameType.Observation)
    obs_frame['SomeDataOrSomething'] = 3
    return [obs_frame, frame]

pipe = core.G3Pipeline()
pipe.Add(just_make_wiring_map)
pipe.Add(add_obs_header)
pipe.Add(core.Dump)
pipe.Add(core.G3MultiFileWriter, filename = 'test-%03u.g3.gz', size_limit = 2 * 1024**3)
#pipe.Add(core.G3Writer, filename = 'test-%03u.g3.gz')
pipe.Run()
