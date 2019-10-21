from spt3g import core, std_processing
from spt3g.mapmaker import mapmakerutils as mm
from spt3g.mapspectra import basicmaputils
import glob

file1 = '/spt/user/ddutcher/coadds/sequential/364_maps.g3'

f1 = core.G3File(file1)
frame150l=f1.next()
frame150r=f1.next()
mm.RemoveWeightModule(frame150l)
mm.RemoveWeightModule(frame150r)
w150l=frame150l['Wpol']
nt150l=np.asarray(frame150l['T'])
nw150l=np.asarray(w150l.TT)
file150l=file1.replace('/spt/user/ddutcher/coadds/','/spt/user/tcrawfor/')
file150l=file150l.replace('.g3','')
tfile150l=file150l.replace('maps','t_map_150l')
nt150l.tofile(tfile150l)
wfile150l=tfile150l.replace('t_map','t_weight')
nw150l.tofile(wfile150l)
w150r=frame150r['Wpol']
nt150r=np.asarray(frame150r['T'])
nw150r=np.asarray(w150r.TT)
file150r=file1.replace('/spt/user/ddutcher/coadds/','/spt/user/tcrawfor/')
file150r=file150r.replace('.g3','')
tfile150r=file150r.replace('maps','t_map_150r')
nt150r.tofile(tfile150r)
wfile150r=tfile150r.replace('t_map','t_weight')
nw150r.tofile(wfile150r)




