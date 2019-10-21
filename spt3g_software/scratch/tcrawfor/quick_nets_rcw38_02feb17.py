from spt3g import core, dfmux
from spt3g.scratch.tcrawfor import tctools
import numpy as np
import pickle

file2 = '/spt_data/bolodata/fullrate/RCW38/2653414/nominal_online_cal.g3'
f2 = core.G3File(file2)
bp = f2.next()['NominalBolometerProperties']

whitenoise_w = pickle.load(open('whitenoise_w_02feb17.pkl'))
names = pickle.load(open('bolo_names_02feb17.pkl'))
wndict = {}
for i in np.arange(len(names)):
    wndict[names[i]] = whitenoise_w[i]

netdict = {}
wunpol = []
file1 = '/home/nwhitehorn/rcw38-secondlight-2.g3'
f1 = core.G3File(file1)
for frame in f1:
    if frame.type is core.G3FrameType.Map:
        if 'Wunpol' in frame and frame['Id'] != 'bs':
            wunpol.append(frame['Wunpol'])
        if frame['Id'] in names:
            net, w_per_k, whitelevel = \
                tctools.quick_net_from_rcw38(frame['T'],wndict[frame['Id']], 
                                             data_is_wn=True, invert_map=True,
                                             band=np.int(bp[frame['Id']].band/10.))
            netdict[frame['Id']] = net


bands = np.asarray([np.int(bp[name].band/10.) for name in names])
nets = np.asarray([netdict[name] for name in names])


