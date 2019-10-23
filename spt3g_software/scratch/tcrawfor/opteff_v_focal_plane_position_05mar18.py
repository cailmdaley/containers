from spt3g import core, dfmux, std_processing, todfilter
from spt3g.util import tctools
import numpy as np
import pickle
from scipy import ndimage

obsid_str = '35123058'

file2 = '/spt/data/bolodata/fullrate/RCW38-pixelraster/'+obsid_str+'/nominal_online_cal.g3'
f2 = core.G3File(file2)
bp = f2.next()['NominalBolometerProperties']
outdict = pickle.load(open('/big_scratch/tcrawfor/script_output/'+obsid_str+'.pkl','r'))

bands = np.asarray([int(bp[name].band/10.) for name in outdict['opteff'].keys()])
xoffs = np.asarray([bp[name].x_offset for name in outdict['opteff'].keys()])
yoffs = np.asarray([bp[name].y_offset for name in outdict['opteff'].keys()])
opteffs = np.asarray([outdict['opteff'][name] for name in outdict['opteff'].keys()])

wh90 = np.where(bands == 90)
wh150 = np.where(bands == 150)
wh220 = np.where(bands == 220)

xoffs90 = xoffs[wh90]
yoffs90 = yoffs[wh90]
opteffs90 = opteffs[wh90]
xoffs150 = xoffs[wh150]
yoffs150 = yoffs[wh150]
opteffs150 = opteffs[wh150]
xoffs220 = xoffs[wh220]
yoffs220 = yoffs[wh220]
opteffs220 = opteffs[wh220]

oe_low = np.arange(6)*0.05
oe_hi = np.arange(6)*0.05 + 0.05
oe_hi[5] = 1.
cols = ['k','b','c','g','y','r']






