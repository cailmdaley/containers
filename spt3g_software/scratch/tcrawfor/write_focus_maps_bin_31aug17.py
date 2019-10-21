from spt3g import core, std_processing, dfmux, calibration, xtalk, todfilter, mapmaker
import glob

files_all = glob.glob('/spt/user/production/calibration/0537-441-pixelraster/*.g3')
files_all.sort()
files = np.asarray(files_all[len(files_all)-7:])
nfiles = len(files)

for file1 in files:
    f1 = core.G3File(file1)
    frame90=f1.next()
    frame150=f1.next()
    frame220=f1.next()
    t150=frame150['T']
    w150=frame150['Wunpol']
    nt150=np.asarray(frame150['T'])
    nw150=np.asarray(w150.TT)
    nt150[np.where(nw150 > 0.)] /= nw150[np.where(nw150 > 0.)]
    file150=file1.replace('/spt/user/production/calibration/0537-441-pixelraster/','/spt/user/tcrawfor/public/focusmaps_v2/')
    file150=file150.replace('.g3','')
    nt150.tofile(file150)





