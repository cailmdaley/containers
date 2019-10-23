#import numpy
#import scipy
#import pickle
#from spt3g import core, std_processing, gcp
#from spt3g.scratch.tcrawfor import tctools
#
#np = numpy
#
#def grab_data(frame, cmap_dict, data1 = [], data2 = [], data3 = []):
#
#    if frame.type != core.G3FrameType.Scan:
#        return
#    try:
#        bdata = []
#        data1.append(frame['BoresightRa'])
#        data2.append(frame['BoresightDec'])
#        for key in cmap_dict:
#            bdata.append(frame['RawTimestreams_I'][key])
#        data3.append(bdata)
#    except:
#        pass
#
#
#f2 = core.G3File('/home/tcrawfor/spt_code/spt3g_software/calibration/scripts/elnod_output_2668507.g3')
#eframe = f2.next()
#names = eframe['ElnodSlopes'].keys()
#elnod_s = np.array([eframe['ElnodSlopes'][key] for key in names])
#elnod_ss = np.array([eframe['ElnodSigmaSlopes'][key] for key in names])
#elnod_ss[(np.where(elnod_ss == 0.))[0]] = 1e12
#n2get = 10
#stemp = np.argsort(elnod_sn)
#rstemp = stemp[::-1]
#bolos2get = (np.asarray(names))[rstemp[0:10]]
#
#file1 = '/spt_data/bolodata/fullrate/0537-441/2668644/0000.g3'
#f1 = core.G3File(file1)
#oframe = f1.next()
#wframe = f1.next()
#wmap = wframe['WiringMap']
#bnames = np.asarray(wmap.keys())
##cframe = f1.next()
##bp = cframe['NominalBolometerProperties']
#
#cmap_dict = {}
#for key in bolos2get:
#    key2 = key
#    wmk = wmap[key2]
#    ilist = [wmk.board_serial, wmk.module, wmk.channel*2]
#    cmap_dict[key2] = ilist
#
#data1 = []
#data2 = []
#data3 = []
#
#for frame in f1:
#    grab_data(frame, cmap_dict, data1 = data1, data2 = data2, data3 = data3)

ra = tctools.list_to_array(data1)/core.G3Units.deg
dec = tctools.list_to_array(data2)/core.G3Units.deg
nframes = len(data1)
nbolos = len(bolos2get)
npts_tot = len(ra)
bdata = np.zeros([nbolos,npts_tot])
npts = 0
for i in np.arange(nframes):
    tnpts = len(data3[i][0])
    for j in np.arange(nbolos):
        bdata[j,npts:npts+tnpts] = data3[i][j]
    npts += tnpts


