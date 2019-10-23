from spt3g import core, coordinateutils
import sys, numpy, tmplfilter

# Usage: filter.py inputmap.h5 mapname output.npz [beam.npz]
#
# If passed a beam file (see processvenustobeammap.py), that will be used as the
# beam kernel. If no beam is passed, this script will use a Gaussian template that
# could produce results differing by as much as 15%. Note that no checks are made
# that the angular resolution of the beam map matches the input map: you have to
# take care of this yourself, at least for now.

mapname = sys.argv[2]
band = int(mapname.split('-')[-1].split('GHz')[0])

mframe = None
for frame in core.G3File(sys.argv[1]):
	if 'Id' in frame and frame['Id'] == mapname:
		mframe = frame
		break
assert(mframe is not None)

m = mframe['T']/core.G3Units.K
w = mframe['Wunpol'].TT

# Conversion to mJy:
# - Use SPT-SZ data release paper table 1 to get 396.3e6 Jy/K/sr
# - Convert sr to arcmin^2 (8.46159496e-8)
# - Convert Jy to mJy (10^3)
# - Multiply by resolution^2 to get per-pixel density
if band == 150:
	factor = 384e6
elif band == 95:
	factor = 384e6 * 28.9/44.4
else:
	assert(False, 'Unknown band!')
tomjy = factor*8.46159496e-8*1000*((m.res/core.G3Units.arcmin)**2)

print('Loaded %s' % sys.argv[1])

def downsample_map(data, N=8):
	from numpy import average, split
	width = data.shape[0]
	height= data.shape[1]
	return average(split(average(split(data, width // N, axis=1), axis=-1), height // N, axis=1), axis=-1)

if False:
	m = downsample_map(m, 2) # Downsample to 0.5-arcminute pixels
	w = downsample_map(w, 2)
	tomjy *= 4
	pixel_width = 2*m.res/core.G3Units.arcmin
else:
	pixel_width = m.res/core.G3Units.arcmin

# Make an annulus of inner radius 2 arcmin and outer radius 5 arcmin
filter = numpy.zeros((45, 45))
r = numpy.hypot(numpy.indices(filter.shape)[0] - filter.shape[0]/2, numpy.indices(filter.shape)[1] - filter.shape[1]/2)
filter[r > 2/pixel_width] = 1
filter[r > 5/pixel_width] = 0

# Calculate annulus to remove 
csigma = w**-.5
filtered, filteredsigma = tmplfilter.tmplfilter(m/w, csigma, filter)

# Remove it and convert to mJy
cm = (m/w - filtered)*tomjy
csigma = numpy.hypot(csigma, filteredsigma)*tomjy

# Rescale to map RMS
csigma *= numpy.std((cm/csigma)[numpy.isfinite(cm/csigma)])

print('Computed raw high-pass-filtered map')

# Deconcolve the instrument beam
def gaussian(height, center_x, center_y, width_x, *args):
	"""Returns a gaussian function with the given parameters"""
	width_x = float(width_x)
	width_y = float(width_x)
	rotation = 0
	if len(args) > 0:
		width_y = float(args[0])
	if len(args) > 1:
		rotation = float(args[1])
	return lambda x,y: height*numpy.exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2 + 2*rotation*(center_x-x)*(center_y-y))/2)/(2*numpy.pi*width_x*width_y)

if int(band) == 150:
	beamsigma = 1.1/2.35482 / pixel_width
elif int(band) == 95:
	beamsigma = 1.7/2.35482 / pixel_width
else:
	assert(False, 'Unknown band')
filterwidth = int(numpy.ceil(3*beamsigma))
psfilter = gaussian(1, filterwidth, filterwidth, beamsigma)(*numpy.indices((2*filterwidth+1,2*filterwidth+1)))

# Use Venus instead (or another raster template) if passed an extra argument
if len(sys.argv) > 4:
	psfilter = numpy.load(sys.argv[4])

psfilter = psfilter - tmplfilter.tmplfilter(psfilter, numpy.ones(psfilter.shape), filter)[0]

psfiltered, psfilteredsigma = tmplfilter.tmplfilter(cm, csigma, psfilter)

print('Computed point-source-filtered map')

time = mframe.get('Observation', None)

numpy.savez(sys.argv[3], map=psfiltered, rms=psfilteredsigma, start=time)
