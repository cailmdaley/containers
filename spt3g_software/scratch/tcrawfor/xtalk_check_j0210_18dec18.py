from spt3g import core, dfmux, mapmaker
from spt3g.util import tctools
import numpy as np
import pickle
import glob
import scipy
from scipy import ndimage

obsid = '46849140'

dir = '/spt/data/bolodata/downsampled/PMNJ0210-5101-pixelraster/'+obsid+'/'
g3files = glob.glob(dir+'0*.g3')
g3files.sort()
nfiles = len(g3files)

poly_order = 4
npts2 = 10

outdict = {}

ndone = 0
for g3file in g3files:
    f1 = core.G3File(g3file)
    nfdone = 0
    for frame in f1:
        if frame.type is core.G3FrameType.Scan:
            if 'Turnaround' not in frame:
                if ndone == 0 and nfdone == 0:
                    names = frame['RawTimestreams_I'].keys()
                    nbolo = len(names)
                    for name in names:
                        tempdict = {}
                        tempdict['maxsn'] = -100.
                        tempdict['xtalk_coeff'] = np.zeros(nbolo)
                        outdict[name] = tempdict
                npts = len(frame['RawTimestreams_I'][names[0]])
                xtemp = np.arange(npts)
                bdtemp = np.zeros([nbolo,npts])
                for q in np.arange(nbolo):
                    tstemp = np.asarray(frame['RawTimestreams_I'][names[q]])
                    tstemp *= -1.
                    coeffs = np.polyfit(xtemp,tstemp,poly_order)
                    for j in np.arange(poly_order+1):
                        tstemp -= coeffs[j]*xtemp**(poly_order-j)
                    bdtemp[q,:] = tstemp
                for q in np.arange(nbolo):
                    tstemp = bdtemp[q,:]
                    tstemp = ndimage.gaussian_filter1d(tstemp,4.)
                    thisimax = np.argmax(np.abs(tstemp))
                    thissn = np.max(np.abs(tstemp))/np.std(tstemp)
                    name = names[q]
                    if thissn > 5. and thissn > outdict[name]['maxsn']:
                        try:
                            tstsq = np.sum(bdtemp[q,thisimax-npts2:thisimax+npts2]**2)
                            xttemp = np.asarray([np.sum(bdtemp[q,thisimax-npts2:thisimax+npts2]*bdtemp[qq,thisimax-npts2:thisimax+npts2])/np.sqrt(tstsq*np.sum(bdtemp[qq,thisimax-npts2:thisimax+npts2]**2)) for qq in np.arange(nbolo)])
                            outdict[name]['xtalk_coeff'] = xttemp
#                            print(notavariable)
                        except:
                            pass
