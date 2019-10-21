from matplotlib import use
use('agg')
import numpy as np
from spt3g import core, gcp, std_processing, util

core.set_log_level(core.G3LogLevel.LOG_INFO, 'ARCTimerangeReader')

pipe = core.G3Pipeline()
pipe.Add(
    std_processing.ARCTimerangeReader,
    start_time=core.G3Time('20190130_001838'),
    stop_time=core.G3Time('20190130_002441'),
    basedir='/spt_data/arc',
    no_more_data_error=False,
)

fields = {
    'tp0': ['array', 'spectrometer', 'tp', 0, 0],
    'tp1': ['array', 'spectrometer', 'tp', 0, 1],
    'tp2': ['array', 'spectrometer', 'tp', 0, 2],
    'tp3': ['array', 'spectrometer', 'tp', 0, 3],
    'utc': ['array', 'spectrometer', 'utc', 2],
    'az': ['antenna0', 'tracker', 'actual', 0],
    'el': ['antenna0', 'tracker', 'actual', 1],
    'utc2': ['antenna0', 'tracker', 'utc']
}

macc = util.extractdata.MultiAccumulator(fields)
pipe.Add(macc)
pipe.Run()

data = macc.extract_values()
np.savez('vlbi_data.npz', **data)



res = 5 * core.G3Units.arcmin

sl = slice(10, 345)

az = data['az'][::100][sl]
el = data['el'][::100][sl]
tp = data['tp0'][sl]


idx = ((az // res) - (az // res).min()).astype(int)
idy = ((el // res) - (el // res).min()).astype(int)
nx = idx.max() + 1
ny = idy.max() + 1
pix = idx + idy * nx
sz = nx * ny

m = np.bincount(pix, weights=tp, minlength=sz).reshape(ny, nx)

from matplotlib import pyplot as plt
plt.imshow(m)
# plt.colorbar()
plt.savefig('vlbi_map.png', bbox_inches='tight')
