from spt3g import core, dfmux, std_processing, todfilter, calibration
from spt3g.util import tctools
import numpy as np
import pickle
from scipy import ndimage

kcmb_conversion_factors = {
  'RCW38': {
    90.0*core.G3Units.GHz: 4.0549662e-07*core.G3Units.K,
    150.0*core.G3Units.GHz: 2.5601153e-07*core.G3Units.K,
    220.0*core.G3Units.GHz: 2.8025804e-07*core.G3Units.K,
  },
  'MAT5A': {
    90.0*core.G3Units.GHz: 2.5738063e-07*core.G3Units.K, # center (608, 555)
    150.0*core.G3Units.GHz: 1.7319235e-07*core.G3Units.K,
    220.0*core.G3Units.GHz: 2.145164e-07*core.G3Units.K,
  },
}

obsid_str = '68216210'

bpframe = core.G3File('/spt/user/production/calibration/boloproperties/60000000.g3').next()
bp = bpframe['BolometerProperties']
wafnames = np.asarray([bp[name].wafer_id for name in bp.keys()])
uwn = np.unique(wafnames)
bandstr = ['90','150','220']
calframe = core.G3File('/spt/user/production/calibration/calframe/ra0hdec-52.25/'+obsid_str+'.g3').next()
names1 = calframe['RCW38FluxCalibration'].keys()
#converter = dfmux.unittransforms.ConvertTimestreamUnits(Input='RawTimestreams_I')
converter = dfmux.ConvertTimestreamUnits(Input='RawTimestreams_I',Output='TimestreamsWatts',Units=core.G3TimestreamUnits.Power)
converter(calframe)

file1 = '/spt/data/bolodata/downsampled/ra0hdec-52.25/' + obsid_str + '/0000.g3'
f3 = core.G3File(file1)
npts_psd = 1024
noisedict1 = {}
noisedict2 = {}
nbdict = {}
for bs in bandstr:
    noisedict1[bs] = {}
    noisedict2[bs] = {}
    nbdict[bs] = {}
    for wn in uwn:
        noisedict1[bs][wn] = np.zeros(npts_psd)
        noisedict2[bs][wn] = np.zeros(npts_psd)
        nbdict[bs][wn] = np.zeros(npts_psd)
samplerate = 76.3
nfnoise = 5
#nfnoise = 2
i = 0
while i < nfnoise:
    frame = f3.next()
    if frame.type is core.G3FrameType.Wiring:
        converter(frame)
    if frame.type is core.G3FrameType.Scan:
        if 'Turnaround' not in frame:
            print(i)
            converter(frame)
            names2 = frame['RawTimestreams_I'].keys()
            names = np.intersect1d(names1,names2)
            for name in names:
                if calframe['CalibratorResponseSN'][name] > 20.:
                    bpt = bp[name]
                    if 'x' in bpt.physical_name:
                        partner = bpt.physical_name.replace('x','y')
                        pnames = [name2 for name2 in frame['RawTimestreams_I'].keys() if bp[name2].physical_name == partner]
                        if len(pnames) > 0:
                            pname = pnames[0]
                            thisbandstr = str(np.int(bpt.band/10.))
                            thiswn = bpt.wafer_id
                            if calframe['CalibratorResponseSN'][pname] > 20. and pname in calframe['RCW38FluxCalibration']:
                                calx = calframe['CalibratorResponse'][name]*calframe['RCW38FluxCalibration'][name]*calframe['RCW38IntegralFlux'][name]/kcmb_conversion_factors['RCW38'][bpt.band]*calframe['RCW38SkyTransmission'][thisbandstr]
                                caly = calframe['CalibratorResponse'][pname]*calframe['RCW38FluxCalibration'][pname]*calframe['RCW38IntegralFlux'][pname]/kcmb_conversion_factors['RCW38'][bpt.band]*calframe['RCW38SkyTransmission'][thisbandstr]
                                bdx = np.asarray(frame['TimestreamsWatts'][name])/calx
                                bdx -= np.mean(bdx)
                                bdy = np.asarray(frame['TimestreamsWatts'][pname])/caly
                                bdy -= np.mean(bdy)
                                bdxsm = ndimage.gaussian_filter1d(bdx,32.)
                                bdysm = ndimage.gaussian_filter1d(bdy,32.)
                                rctemp = np.sum(bdxsm*bdysm)/np.sum(bdxsm*bdxsm)
                                bdx *= np.sqrt(rctemp)
                                bdy /= np.sqrt(rctemp)
                                psdd1 = tctools.quick_pspec((bdx+bdy)/2./core.G3Units.K,rate=samplerate,npts_psd=npts_psd)
                                psdd2 = tctools.quick_pspec((bdx-bdy)/2./core.G3Units.K,rate=samplerate,npts_psd=npts_psd)
                                if np.max(np.abs(psdd1['psd'])) > 0. and np.max(np.abs(psdd2['psd'])) > 0.:
                                    noisedict1[thisbandstr][thiswn] += psdd1['psd']**2/2.
                                    noisedict2[thisbandstr][thiswn] += psdd2['psd']**2/2.
                                    nbdict[thisbandstr][thiswn] += 1.
            i+=1

for bs in bandstr:
    for wn in uwn:
        noisedict1[bs][wn] /= nbdict[bs][wn]
        noisedict2[bs][wn] /= nbdict[bs][wn]
        noisedict1[bs][wn] = np.sqrt(noisedict1[bs][wn])
        noisedict2[bs][wn] = np.sqrt(noisedict2[bs][wn])

