'''
weight_stats.py

Module to compute statistics about maps
'''

import numpy as np
import argparse
import pickle as pkl
from spt3g import core, mapmaker

parser = argparse.ArgumentParser()
parser.add_argument('input_files')
parser.add_argument('-b','--band', default = '150GHz')
parser.add_argument('-o','--output', default='output.pkl')
pargs = parser.parse_args()

m = core.G3File(pargs.input_files)

right = core.G3Frame(core.G3FrameType.Map)
left = core.G3Frame(core.G3FrameType.Map)
for fr in m:
    if fr.type == core.G3FrameType.Map:
        if 'Right-'+pargs.band in fr['Id']:
            for k in fr.keys():
                if k not in right:
                    right[k] = fr[k]
                else:
                    val = right[k]
                    del right[k]
                    right[k]=val+fr[k]
        if 'Left-'+pargs.band in fr['Id']:
            for k in fr.keys():
                if k not in left:
                    left[k] = fr[k]
                else:
                    val = left[k]
                    del left[k]
                    left[k]=val+fr[k]
    if 'CalibratorResponseSN' in fr:
        calsn = np.array(fr['CalibratorResponseSN'].values())
        calsn = calsn[np.where(np.isfinite(calsn))]

combined = core.G3Frame(core.G3FrameType.Map)
difference = core.G3Frame(core.G3FrameType.Map)

combined['T'] = left['T']+right['T']
combined['Q'] = left['Q']+right['Q']
combined['U'] = left['U']+right['U']
combined['Wpol'] = left['Wpol']+right['Wpol']

lt,lq,lu = mapmaker.mapmakerutils.remove_weight(left['T'], left['Q'], left['U'], left['Wpol'])
rt,rq,ru = mapmaker.mapmakerutils.remove_weight(right['T'], right['Q'], right['U'], right['Wpol'])
t,q,u = mapmaker.mapmakerutils.remove_weight(combined['T'], combined['Q'], combined['U'], combined['Wpol'])

difference['T'] = (lt - rt)/2.
difference['Q'] = (lq - rq)/2.
difference['U'] = (lu - ru)/2.
difference['Wpol'] = left['Wpol']+right['Wpol']

results = dict()
w = np.asarray(combined['Wpol'].TT)
w = w[np.where(w != 0.0)]
m = np.asarray(combined['T']/combined['Wpol'].TT)

results['nbolos'] = len(np.where(calsn>20)[0])
results['sum_w'] = np.sum(w)
results['std_w'] = np.std(w)
results['med_w'] = np.median(w)
results['tmap_rms'] = np.nanstd(np.asarray(t))
results['pmap_rms'] = np.mean([np.nanstd(q),
                               np.nanstd(u)])
results['tnoise_rms'] = np.nanstd(difference['T'])
results['pnoise_rms'] = np.mean([np.nanstd(difference['Q']), 
                                np.nanstd(difference['U'])])

pkl.dump(results, open(pargs.output,'wb'))