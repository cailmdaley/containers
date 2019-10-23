from spt3g import core, dfmux, std_processing
from spt3g.util import tctools
import pickle

zpdata = pickle.load(open('/home/panz/s_n_20_6hz.pkl'))
#tcdata = pickle.load(open('/spt/user/tcrawfor/public/taus_for_bender.pkl'))
tcdata = pickle.load(open('/spt/user/tcrawfor/public/taus_from_zp_obs.pkl'))

znames = zpdata.keys()
ntau = len(znames)

tnames = tcdata.keys()
cfile = '/spt/data/bolodata/fullrate/calibrator/18576429/nominal_online_cal.g3'
f2 = core.G3File(cfile)
bp = f2.next()['NominalBolometerProperties']
tpnames = np.asarray([(bp[name].physical_name).replace('_','.') for name in tnames])

ztaus = np.zeros(ntau)
ttaus = np.zeros(ntau)
for i in np.arange(ntau):
    zn = znames[i]
    ztaus[i] = zpdata[zn][0]
    whn = np.where(tpnames == str(zn))
    if len(whn[0]) > 0:
        ttaus[i] = tcdata[tnames[whn[0][0]]]['tau']
