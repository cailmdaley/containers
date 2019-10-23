from spt3g import core, dfmux, std_processing, todfilter
from spt3g.util import tctools
import numpy as np
import pickle
from scipy import ndimage

outdict = pickle.load(open('/spt/user/tcrawfor/public/cals_rcw38_mat5a_08mar19.pkl'))
#outdict = pickle.load(open('/spt/user/tcrawfor/public/cals_rcw38_mat5a_2018_08mar19.pkl'))
names = [name for name in outdict.keys() if 'opteff' in outdict[name]['RCW38']]
#cframe = core.G3File('/spt/user/production/calibration/calframe/RCW38/68192210.g3').next()
#bp = cframe['BolometerProperties']
bpframe = core.G3File('/spt/user/production/calibration/boloproperties/60000000.g3').next()
#bpframe = core.G3File('/spt/user/production/calibration/boloproperties/40000000.g3').next()
bp = bpframe['BolometerProperties']
wafnames = np.asarray([bp[name].wafer_id for name in names])

bands = np.zeros(len(names))
xoffs = np.zeros(len(names))
yoffs = np.zeros(len(names))
opteffs = np.zeros(len(names))
for i in np.arange(len(names)):
    name = names[i]
    bands[i] = bp[name].band/10.
    xoffs[i] = bp[name].x_offset/core.G3Units.arcmin
    yoffs[i] = bp[name].y_offset/core.G3Units.arcmin
    if 'y' in bp[name].physical_name:
        xoffs[i] += 1.
        yoffs[i] += 1.
    opteffs[i] = outdict[name]['RCW38']['opteff']

wh90 = np.where(bands == 90.)
wh150 = np.where(bands == 150.)
wh220 = np.where(bands == 220.)

xoffs90 = xoffs[wh90]
yoffs90 = yoffs[wh90]
opteffs90 = opteffs[wh90]
xoffs150 = xoffs[wh150]
yoffs150 = yoffs[wh150]
opteffs150 = opteffs[wh150]
xoffs220 = xoffs[wh220]
yoffs220 = yoffs[wh220]
opteffs220 = opteffs[wh220]

cols = ['k','b','c','g','y','r']

#do90 = True
do90 = False
#do150 = False
do150 = True
do220 = False
#do220 = True
frac2plot = 0.6
nbins = len(cols)
binsize = frac2plot/np.float(nbins)
oe_low_frac = 1. + (np.arange(nbins)-np.float(nbins)/2.)*binsize
oe_hi_frac = oe_low_frac + binsize
if do90:
    xtemp = xoffs90
    ytemp = yoffs90
    oetemp = opteffs90*2.
    bandstr = '90'
if do150:
    xtemp = xoffs150
    ytemp = yoffs150
    oetemp = opteffs150*2.
    bandstr = '150'
if do220:
    xtemp = xoffs220
    ytemp = yoffs220
    oetemp = opteffs220*2.
    bandstr = '220'
oe_low = oe_low_frac*np.nanmedian(oetemp)
oe_hi = oe_hi_frac*np.nanmedian(oetemp)
oe_low[0] = 0.
oe_hi[5] = 1.0

for i in np.arange(6):
#for i in np.arange(6)[::-1]:
    whtemp = np.where(np.logical_and(oetemp > oe_low[i], oetemp < oe_hi[i]))
#    plot(xtemp[whtemp],-ytemp[whtemp],'o',color=cols[i])
    plot(xtemp[whtemp],-ytemp[whtemp],'o',color=cols[i],markersize=4)
#    plot(xtemp[whtemp],-ytemp[whtemp],'x',color=cols[i])
    if i < nbins/2:
        plot(7000.,0.,'o',color=cols[i],label="median - " + str(np.int(np.round((1.-(oe_low_frac[i]+oe_hi_frac[i])/2.)*100.))) + "%")
    else:
        plot(7000.,0.,'o',color=cols[i],label="median + " + str(np.int(np.round(((oe_low_frac[i]+oe_hi_frac[i])/2. - 1.)*100.))) + "%")
xlim(-88,94)
ylim(-64,64)
xlabel('X offset from boresight [arcmin]')
ylabel('Y offset from boresight [arcmin]')
title(bandstr + ' GHz')
tctools.annotate_focal_plane()
legend(loc=4,fontsize='small')
#savefig('opteff_v_focal_plane_position_'+bandstr+'_10mar19.png')
