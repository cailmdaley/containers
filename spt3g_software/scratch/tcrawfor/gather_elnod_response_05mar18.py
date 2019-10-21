from spt3g import core, std_processing, dfmux, calibration
import pickle, glob

bp = core.G3File('/spt/user/production/calibration/boloproperties/32512242.g3').next()['BolometerProperties']
efiles = glob.glob('/spt/user/production/calibration/elnod/36*.g3')
efiles.sort()
elnod_response = {}
for name in bp.keys():
    elnod_response[name] = []
for efile in efiles:
    eframe = core.G3File(efile).next()
    enr = eframe['ElnodSlopes']
    if np.abs(np.median(np.asarray(enr.values()))) > 50.:
        for name in enr.keys():
            try:
                elnod_response[name].append(enr[name])
            except:
                pass

med_enr = {}
for name in elnod_response.keys():
    med_enr[name] = -np.median(np.asarray(elnod_response[name]))
        

bands = np.asarray([int(bp[name].band/10.) for name in med_enr.keys()])
xoffs = np.asarray([bp[name].x_offset for name in med_enr.keys()])
yoffs = np.asarray([bp[name].y_offset for name in med_enr.keys()])
enrs = np.asarray([med_enr[name] for name in med_enr.keys()])

wh90 = np.where(bands == 90)
wh150 = np.where(bands == 150)
wh220 = np.where(bands == 220)

xoffs90 = xoffs[wh90]
yoffs90 = yoffs[wh90]
enrs90 = enrs[wh90]
xoffs150 = xoffs[wh150]
yoffs150 = yoffs[wh150]
enrs150 = enrs[wh150]
xoffs220 = xoffs[wh220]
yoffs220 = yoffs[wh220]
enrs220 = enrs[wh220]

oe_low = np.arange(6)*0.05
oe_hi = np.arange(6)*0.05 + 0.05
oe_hi[5] = 1.
cols = ['k','b','c','g','y','r']


