from spt3g import core
import matplotlib.pyplot as plt
from spt3g.util import bolo_psd_plot


class ViewTimestreamPSD(core.G3Module):
    """
    If bolo_order is specified, individually plots PSDs for each bolo listed.
    Otherwise, plots PSD for every bolometer in scan on a heatmap.
    """
    def __init__(self, in_ts_key="CalTimestreams",
                 bolo_order=None, low_f=None, high_f=None):
        super(ViewTimestreamPSD, self).__init__()
        self.in_ts_key = in_ts_key
        self.bolo_order = bolo_order
        self.low_f = low_f
        self.high_f = high_f
    def Process(self, frame):
        if (frame.type==core.G3FrameType.Scan):
            ts = frame[self.in_ts_key]
            bolo_psd_plot(ts, bolo_order=self.bolo_order,
                          low_f=self.low_f,
                          high_f=self.high_f)
            plt.show()
