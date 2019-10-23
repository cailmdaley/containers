## gather a bunch of stats on cal and elnod s/n
#
#import glob
#from spt3g.util import tctools
#from spt3g import core, std_processing, dfmux, gcp
#
#cdir = '/spt/analysis/production/calibration/'
#calfiles = glob.glob(cdir+'calibrator/*.g3')
##calfiles = calfiles[0:10]
#cdict = {}
#for calfile in calfiles:
#    dtemp = {}
#    strtemp = filter(None,calfile.split('/'))
#    calfile_short = strtemp[len(strtemp)-1]
#    thisobsid = (filter(None,calfile_short.split('.g3')))[0]
#    thistime = tctools.obsid_to_g3time(thisobsid)
#    dtemp['g3time'] = thistime
#    dtemp['mjd'] = thistime.mjd
#    try:
#        f1 = core.G3File(calfile)
#        cframe = f1.next()
#        calsns = np.asarray((cframe['CalibratorResponseSN']).values())
#    except:
#        core.log_warn('No good frames in '+calfile+'.\n')
#        continue
#    calsns[np.where(np.isfinite(calsns) == False)] = 0.
#    wh5 = np.where(calsns > 5.)
#    dtemp['n5'] = len(wh5[0])
#    wh20 = np.where(calsns > 20.)
#    dtemp['n20'] = len(wh20[0])
#    cdict[np.int(thisobsid)] = dtemp
#
#cmjds = np.asarray([cdict[key]['mjd'] for key in cdict.keys()])
#cn5 = np.asarray([cdict[key]['n5'] for key in cdict.keys()])
#cn20 = np.asarray([cdict[key]['n20'] for key in cdict.keys()])
#scmjd = np.argsort(cmjds)
#cmjds = cmjds[scmjd]
#cn5 = cn5[scmjd]
#cn20 = cn20[scmjd]

edir = '/spt/analysis/production/calibration/'
elnodfiles = glob.glob(cdir+'elnod/*.g3')
#elnodfiles = elnodfiles[0:2]
edict = {}
for elnodfile in elnodfiles:
    dtemp = {}
    strtemp = filter(None,elnodfile.split('/'))
    elnodfile_short = strtemp[len(strtemp)-1]
    thisobsid = (filter(None,elnodfile_short.split('.g3')))[0]
    thistime = tctools.obsid_to_g3time(thisobsid)
    dtemp['g3time'] = thistime
    dtemp['mjd'] = thistime.mjd
    try:
        f1 = core.G3File(elnodfile)
        eframe = f1.next()
        elnodsns = -np.asarray((eframe['ElnodSNSlopes']).values())
    except:
        core.log_warn('No good frames in '+elnodfile+'.\n')
        continue
    elnodsns[np.where(np.isfinite(elnodsns) == False)] = 0.
    wh5 = np.where(elnodsns > 100.)
    dtemp['n5'] = len(wh5[0])
    wh20 = np.where(elnodsns > 20.)
    dtemp['n20'] = len(wh20[0])
    edict[np.int(thisobsid)] = dtemp

emjds = np.asarray([edict[key]['mjd'] for key in edict.keys()])
en5 = np.asarray([edict[key]['n5'] for key in edict.keys()])
en20 = np.asarray([edict[key]['n20'] for key in edict.keys()])
semjd = np.argsort(emjds)
emjds = emjds[semjd]
en5 = en5[semjd]
en20 = en20[semjd]

