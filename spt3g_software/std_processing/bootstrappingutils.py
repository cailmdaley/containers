from spt3g import core, dfmux
from spt3g.calibration import BolometerPropertiesMap, BolometerProperties


@core.indexmod
class MakeBoresightBolometerProperties(object):
    '''
    Creates a bolometer properties map with all pointing offsets set to zero. Renames
    the previous bolometer properties map (if present) to 'RealBolometerProperties'
    '''
    def __init__(self):
        core.log_warn("This is deprecated, please move away")
        self.wiring_map = None
        self.sent_off = False
        self.cached_cal_frame = core.G3Frame(core.G3FrameType.Calibration)

    def __call__(self, frame):
        if frame.type == core.G3FrameType.Scan:
            if self.sent_off or not self.wiring_map:
                return

            # No cal data yet, fabricate from scratch
            calframe = self.cached_cal_frame
            bolo_props = BolometerPropertiesMap()
            for k in self.wiring_map.keys():
                wm = self.wiring_map[k]
                bp = BolometerProperties()
                bp.physical_name = k
                bp.band = 0
                bp.pol_angle = 0
                bp.pol_efficiency = 1
                bp.x_offset = 0
                bp.y_offset = 0
                bp.wafer_id = ''
                bp.squid_id = ''
                bolo_props[k] = bp
            calframe['BolometerProperties'] = bolo_props
            self.sent_off = True
            return [calframe, frame]
        elif frame.type == core.G3FrameType.Wiring:
            self.wiring_map = frame['WiringMap']
        elif frame.type == core.G3FrameType.Calibration:
            if 'BolometerProperties' in frame:
                new_bolo_props = BolometerPropertiesMap(frame['BolometerProperties'])
                for bp in new_bolo_props.values():
                    bp.x_offset = 0
                    bp.y_offset = 0
                frame['RealBolometerProperties'] = frame['BolometerProperties']
                del frame['BolometerProperties']
                frame['BolometerProperties'] = new_bolo_props
            elif not self.sent_off:
                # Hrm. Early calibration frame, but no data to put in it yet.
                # Hold off on everything for a bit
                self.cached_cal_frame = frame
                return []
