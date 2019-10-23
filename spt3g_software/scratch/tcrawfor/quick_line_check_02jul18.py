from spt3g import core, dfmux, std_processing, util
import glob, pickle

npts = 7000
freqs = np.arange(npts/2)/(npts/2.)*76.3/2.

cf = core.G3File('/spt/user/production/calibration/calframe/ra0hdec-44.75/47084030.g3').next()
bp = cf['BolometerProperties']
cnames = np.intersect1d(cf['CalibratorResponseSN'].keys(),bp.keys())
goodnames = np.asarray([name for name in cnames if cf['CalibratorResponseSN'][name] > 20. and bp[name].band == 1500.0])
allwafs = np.asarray([bp[name].wafer_id for name in goodnames])
uwafs = np.unique(allwafs)
uwafs.sort()

leftdict = {}
rightdict = {}
for name in goodnames:
    leftdict[name] = np.zeros(npts/2)
    rightdict[name] = np.zeros(npts/2)

f1 = core.G3File('/spt/data/bolodata/downsampled/ra0hdec-44.75/47084030/0000.g3')
for i in np.arange(20):
    frame = f1.next()
    if frame.type is core.G3FrameType.Scan:   
        if 'Turnaround' not in frame:
            daz = np.mean(np.diff(frame['OnlineBoresightAz']))
            if daz < 0.:
                for name in goodnames:
                    dtemp = frame['RawTimestreams_I'][name][0:npts]
                    dtemp -= np.mean(dtemp)
                    dtemp *= np.hanning(npts)
                    dtempf = np.abs(np.fft.fft(dtemp))
                    leftdict[name] += dtempf[0:npts/2]**2
