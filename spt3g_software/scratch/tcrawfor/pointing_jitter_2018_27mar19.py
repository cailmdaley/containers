from spt3g import core, std_processing, dfmux, calibration
from spt3g.mapmaker import mapmakerutils as mm  
from spt3g.util import gauss_fit
from astropy.io import ascii
from scipy import ndimage 
import glob, pickle

field = 'ra0hdec-44.75'
#field = 'ra0hdec-52.25'
#field = 'ra0hdec-59.75'
#field = 'ra0hdec-67.25'

if field == 'ra0hdec-44.75':
    yp1 = 4280
    xp1 = 6391
    source = 'RCW38'
if field == 'ra0hdec-52.25':
    yp1 = 3201
    xp1 = 2128
    source = 'RCW38'
if field == 'ra0hdec-59.75':
    yp1 = 1641
    xp1 = 1986
    source = 'MAT5A'
if field == 'ra0hdec-67.25':
    yp1 = 1245
    xp1 = 5734
    source = 'MAT5A'

dddir1 = '/spt/user/ddutcher/'+field+'/lpf18k_150t_0.5res_online_20190308/'
files1 = glob.glob(dddir1+'*/*.g3')
files1.sort()
nfiles = len(files1)
## !!!
#nfiles = 2
## !!!

xcens1 = np.zeros(nfiles)
ycens1 = np.zeros(nfiles)
xcens2 = np.zeros(nfiles)
ycens2 = np.zeros(nfiles)
ycoll2 = np.zeros(nfiles)
eltilt2 = np.zeros(nfiles)
cutouts1 = np.zeros([nfiles,40,40])
cutouts2 = np.zeros([nfiles,40,40])

for i in np.arange(nfiles):
    file1 = files1[i]
    print(file1)
    f1 = core.G3File(file1)
    for frame in f1:
        if frame.type is core.G3FrameType.Map:
            if '150' in frame['Id']:
                mm.RemoveWeightModule(frame)
                cutouts1[i,:,:] = np.asarray(frame['T'])[yp1-20:yp1+20,xp1-20:xp1+20]
                ptemp1 = np.asarray(gauss_fit.fit2Dgaussian(cutouts1[i,:,:]))
                ycens1[i] = ptemp1[2]
                xcens1[i] = ptemp1[3]

    file2 = file1.replace('online','offline')
    print(file2)
    f2 = core.G3File(file2)
    for frame in f2:
        if frame.type is core.G3FrameType.Calibration:
            ycoll2[i] = frame[source+'PointingModelCorrection']['fixedCollimation'][1]
            eltilt2[i] = frame[source+'PointingModelCorrection']['tilts'][2]
        if frame.type is core.G3FrameType.Map:
            if '150' in frame['Id']:
                mm.RemoveWeightModule(frame)
                cutouts2[i,:,:] = np.asarray(frame['T'])[yp1-20:yp1+20,xp1-20:xp1+20]
                ptemp2 = np.asarray(gauss_fit.fit2Dgaussian(cutouts2[i,:,:]))
                ycens2[i] = ptemp2[2]
                xcens2[i] = ptemp2[3]

outdict = {}
outdict['xcens1'] = xcens1
outdict['ycens1'] = ycens1
outdict['xcens1'] = xcens1
outdict['ycens1'] = ycens1
outdict['ycoll2'] = ycoll2
outdict['eltilt2'] = eltilt2
outdict['cutouts1'] = cutouts1
outdict['cutouts2'] = cutouts2
outfile = '/spt/user/tcrawfor/public/pointing_jitter_one_source_' + field + '.pkl'
pickle.dump(outdict,open(outfile,'w'))
