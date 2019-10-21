from spt3g import core, dfmux, std_processing
from spt3g.util import tctools
import numpy as np
import pickle
from scipy import ndimage
import glob

cfile = '/spt/data/bolodata/fullrate/calibrator/18576429/nominal_online_cal.g3'
f2 = core.G3File(cfile)
bp = f2.next()['NominalBolometerProperties']
names = np.asarray(bp.keys())
# !!!
#names = names[10:12]
# !!!
nbolo = len(names)
bands = np.asarray([bp[name].band/10. for name in names])

files_all=glob.glob('/spt/data/bolodata/fullrate/calibrator/1762*/0000.g3')
files_all.sort()
files = files_all[0:6]
calfreqs = [6.,9.,15.,30.,52.,65.]

npts_assumed = 27922
npts_chunk = 4096
nchunk = 5
npts_psd = 512
freqs = np.arange(npts_psd)/np.float(npts_psd)*152.6/2.
whcf = []
for cf in calfreqs:
    whcf.append(np.argmin(np.abs(freqs-cf)))
whcf = np.asarray(whcf)

crdict = {}
dcrdict = {}
psddict = {}
for name in names:
    crdict[name] = np.zeros(6)
    dcrdict[name] = np.zeros(6)
    psddict[name] = np.zeros([6,npts_psd])

ifile = 0
for file1 in files:
    f1 = core.G3File(file1)
    frame = f1.next()
    frame = f1.next()
    frame = f1.next()
    for name in names:
        try:
            bdtemp = frame['RawTimestreams_I'][name]
            bdtemp = bdtemp[1000:]
            crtemp = np.zeros(nchunk)
            noisetemp = np.zeros(nchunk)
            for j in np.arange(nchunk):
                bdt2 = bdtemp[npts_chunk*j:npts_chunk*(j+1)]
                bdt2 -= np.mean(bdt2)
                psdtemp = tctools.quick_pspec(bdt2,npts_psd=npts_psd)
                psd = psdtemp['psd']
                crtemp[j] = np.sqrt(np.sum(psd[whcf[ifile]-2:whcf[ifile]+3]**2))
                noisetemp[j] = np.sqrt(np.sum(psd[whcf[ifile]-10:whcf[ifile]-5]**2))
                psddict[name][ifile,:] += psd**2
            psddict[name][ifile,:] = np.sqrt(psddict[name][ifile,:]/np.float(nchunk))
            ctemp = np.mean(crtemp**2)
            ntemp = np.mean(noisetemp**2)
            if ctemp > ntemp:
                crdict[name][ifile] = np.sqrt(ctemp-ntemp)
            else:
                crdict[name][ifile] = 0.
            dcrdict[name][ifile] = np.std(crtemp)
        except:
            dcrdict[name][ifile] = 1e6
    ifile += 1

tau = {}
dtau = {}
for name in names:
    rtemp = crdict[name][4]/crdict[name][0]
    tau[name] = np.sqrt((rtemp**2-1.)/4./np.pi**2/(calfreqs[0]**2-calfreqs[4]**2*rtemp**2))
    dtfrac = np.sqrt((dcrdict[name][4]/crdict[name][4])**2+(dcrdict[name][0]/crdict[name][0])**2)
    dtau[name] = tau[name]*dtfrac

taus = np.asarray([tau[name] for name in names])
dtaus = np.asarray([dtau[name] for name in names])


