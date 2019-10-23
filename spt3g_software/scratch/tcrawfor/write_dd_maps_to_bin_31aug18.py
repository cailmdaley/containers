from spt3g import core, std_processing
from spt3g.mapmaker import mapmakerutils as mm
from spt3g.mapspectra import basicmaputils
import glob

files_all = glob.glob('/spt/user/ddutcher/coadds/lr_wafCM_poly19mhpf300_*.g3')
files_all.sort()
nfiles = len(files_all)
print files_all

for file1 in files_all:
    print file1
    f1 = core.G3File(file1)
    frame90=f1.next()
    frame150=f1.next()
#    frame220=f1.next()
    mm.RemoveWeightModule(frame150)
    w150=frame150['Wpol']
    nt150=np.asarray(frame150['T'])
    qtemp,utemp = basicmaputils.flatten_pol(frame150['Q'],frame150['U'])
    nq150=np.asarray(qtemp)
    nu150=np.asarray(utemp)
    nw150=np.asarray(w150.TT)
    file150=file1.replace('/spt/user/ddutcher/coadds/','/spt/user/tcrawfor/')
    file150=file150.replace('.g3','')
    tfile150=file150.replace('lr_','tlr_')
    nt150.tofile(tfile150)
    qfile150=file150.replace('lr_','qlr_')
    nq150.tofile(qfile150)
    ufile150=file150.replace('lr_','ulr_')
    nu150.tofile(ufile150)
    wfile150=file150.replace('lr','lrwts')
    nw150.tofile(wfile150)




