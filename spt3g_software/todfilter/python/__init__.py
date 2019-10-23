from spt3g.core.load_pybindings import load_pybindings
load_pybindings(__name__, __path__)

from . import dftutils
from . import polyutils
from . import downsampler
from . import notchfilter
from . import util
from . import convolve_filter
from . import timeconstant_filter

try:
    import matplotlib.pyplot
    mplavail = True
except:
    # Allow this to fail if no graphics are available
    mplavail = False

if mplavail:
    from .timestream_psd_modules import ViewTimestreamPSD

del mplavail

