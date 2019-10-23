import numpy
import scipy
from spt3g import core, std_processing, gcp, dfmux

def grab_bdata(frame, cmap_dict = {}, data1 = [], data2 = [], data3 = []):

    if frame.type != core.G3FrameType.Timepoint:
        return
    try:
        bdata = []
        data1.append(frame['CalibratorOn'])
        data2.append(frame['EventHeader'])
        for key in cmap_dict:
            inds = cmap_dict[key]
            bdata.append(frame['DfMux'][inds[0]][inds[1]][inds[2]])
        data3.append(bdata)
    except:
        pass

def get_one_file(bolos2get, gfile='/poleanalysis/sptdaq/20170129_rawdata/20170129_202437.g3'):

    f1 = core.G3File(gfile)
    wframe = f1.next()
    wmap = wframe['WiringMap']
    bnames = np.asarray(wmap.keys())
    cframe = f1.next()
    bp = cframe['NominalBolometerProperties']

    cmap_dict = {}
    for key in bolos2get:
        fkey = filter(None,key.split('.'))
        wafername = fkey[1]
        key2 = wafername + '/' + key
        wmk = wmap[key2]
        ilist = [wmk.board_serial, wmk.module, wmk.channel*2]
        cmap_dict[key2] = ilist
    names = cmap_dict.keys()

    data1 = []
    data2 = []
    data3 = []

    pipe1 = core.G3Pipeline()
    pipe1.add(grab_bdata, frame, cmap_dict = cmap_dict, data1 = data1, data2 = data2, data3 = data3)
    pipe.Add(dfmux.ConvertTimestreamUnits, Input=args.input_ts, Output='BoloMapTimestreams', Units=core.G3TimestreamUnits.Watts)


bdata = np.zeros([len(data3[0]),len(data3)])
for i in np.arange(len(data3)):
    bdata[:,i] = data3[i]

npts_psd = 1024
win = np.hanning(npts_psd)
freqs = np.arange(npts_psd/2)/np.float(npts_psd/2)*76.5
psds = np.zeros([len(names),npts_psd])
for i in np.arange(len(names)):
    thisbdata = bdata[i,:] - np.mean(bdata[i,:])
    for j in np.arange(len(thisbdata)/1024):
        psds[i,:] += (np.abs(np.fft.fft(thisbdata[npts_psd*j:npts_psd*(j+1)])))**2
psds = np.sqrt(psds)

whwhite = np.where(np.logical_and(freqs > 10., freqs < 20.))
whitenoise_cts = np.zeros(len(names))
for i in np.arange(len(names)):
    psdtemp = psds[i,:]
    whitenoise_cts[i] = np.sqrt(np.mean(psdtemp[whwhite]**2))
