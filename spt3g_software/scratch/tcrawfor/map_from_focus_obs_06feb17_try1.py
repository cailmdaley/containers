#import numpy as np
#import scipy
#from scipy import ndimage
#import pickle
#import argparse as ap
#from spt3g import core, std_processing, gcp, dfmux
#import os
#import glob
#
## offsets
#f1 = core.G3File('/home/nwhitehorn/saturnbpm2.g3')
#bpframe = f1.next()
#bp = bpframe['BolometerProperties']
## "R.A. detector is R.A. boresight *plus* offset." (NW, 05Feb17, 9:40pm)
#
## other stuff
#f3 = core.G3File('/spt_data/bolodata/fullrate/elnod/3003099/nominal_online_cal.g3')
#cframe = f3.next()
#bp2 = cframe['NominalBolometerProperties']
#
## get decent bolos
#gnames = pickle.load(open('gnames2_06feb17.pkl'))
#
#start_mjd = 57789.778981481
#stop_mjd = 57789.83125
#
#files_all = glob.glob('/buffer/bolodata/20170205*.g3')
#files = []
#filemjds = []
#for file1 in files_all:
#    file2 = filter(None,file1.split('/'))[2]
#    dstring = filter(None,file2.split('.'))[0]
#    thismjd = core.G3Time(dstring).mjd
#    if thismjd >= start_mjd and thismjd <= stop_mjd:
#        files.append(file1)
#        filemjds.append(thismjd)
#
#files.sort()
#filemjds.sort()
    
map_weighted = np.zeros([720,720])
nhits = np.zeros([720,720])
mapshape = np.shape(map_weighted)
reso_arcmin = 0.25

def grab_boards(f, boardlist = []):
    key = 'antenna0'
    try:
        boardlist.append(f[key])
    except:
        pass

arcdir='/spt_data/arc/'

for file1, fmjd in zip(files, filemjds):
    data4 = []
    start_time = core.G3Time()
    start_time.mjd = fmjd
    stop_time = core.G3Time()
    stop_time.mjd = fmjd + 49./86400.
    pipe1 = core.G3Pipeline()
    pipe1.Add(std_processing.ARCTimerangeReader, start_time=start_time, stop_time=stop_time, basedir=arcdir)
    pipe1.Add(gcp.ARCExtract)
    pipe1.Add(grab_boards, boardlist=data4)
    pipe1.Run()
    ra_src = data4[0]['tracker']['equat_geoc'][0][0]/core.G3Units.hours*15.
    dec_src = data4[0]['tracker']['equat_geoc'][1][0]/core.G3Units.deg
    az = np.zeros(len(data4)*100)
    el = np.zeros(len(data4)*100)
    lst = np.zeros(len(data4)*100)
    elmjds = np.zeros(len(data4)*100)
    for i in np.arange(len(data4)):
        az[100*i:100*(i+1)] = data4[i]['tracker']['actual'][0]
        el[100*i:100*(i+1)] = data4[i]['tracker']['actual'][1]
        lst[100*i:100*(i+1)] = data4[i]['tracker']['lst'][0]
        for j in np.arange(100):
            elmjds[100*i+j] = (data4[i]['tracker']['utc'][0][j]).mjd

    lst_hr = lst/core.G3Units.h*60. # I think lst is wrong in cal file
    lst_deg = lst_hr*15.
    ra = np.mod(az/core.G3Units.deg+lst_deg,360.)
    dec = -el/core.G3Units.deg
    azrate = np.diff(az/core.G3Units.deg)*100. # actual archived rate is hella-noisy
    azrate2 = np.zeros(len(az))
    azrate2[0:len(azrate)] = azrate
    whg = np.where(np.logical_and(np.abs(np.abs(azrate2)-0.4) < 0.1, 
                                  np.abs(ra-84.7) < 1.7))

    if len(whg[0]) > 1000:
        print(file1)
        mjds = []
        bdata = []
        meanoff = dec_src - np.mean(dec)
        if np.abs(meanoff) > 0.9:
            continue
        gnames2 = []
        for name in gnames:
            try:
                if np.abs(bp[name].y_offset/core.G3Units.deg - meanoff) < 10./60.:
                    gnames2.append(name)
            except KeyError:
                pass
        if len(gnames2) > 0:
            f10 = core.G3File(file1)
            for fr in f10:
                if fr.type == core.G3FrameType.Wiring:
                    wmap = fr['WiringMap']
                if fr.type == core.G3FrameType.Timepoint:
                    mjds.append(fr['EventHeader'].mjd)
                    bdtemp = []
                    for name in gnames2:
                        wmk = wmap[name]
                        inds = [wmk.board_serial, wmk.module, wmk.channel*2]
                        bdtemp.append(fr['DfMux'][inds[0]][inds[1]][inds[2]])
                    bdata.append(bdtemp)

            bdata2 = np.zeros([len(gnames2),len(bdata)])
            for i in np.arange(len(bdata)):
                bdtemp = bdata[i]
                bdata2[:,i] = bdtemp
            xtemp = (np.arange(len(bdata))).astype(float)
            xtemp -= np.mean(xtemp)
            for i in np.arange(len(gnames2)):
                bdtemp = bdata2[i,:]
                ptemp = scipy.polyfit(xtemp,bdtemp,5)
                yfit = np.zeros(len(xtemp))
                for p in np.arange(6):
                    yfit += ptemp[p]*xtemp**(5-p)
                bdata2[i,:] = bdtemp - yfit

            mjds = np.asarray(mjds)
            ra_interp = np.interp(mjds-elmjds[0], elmjds-elmjds[0], ra)
            dec_interp = np.interp(mjds-elmjds[0], elmjds-elmjds[0], dec)
            xp_bs = (np.round((ra_interp-ra_src)*np.cos(dec_interp*np.pi/180.)*60./reso_arcmin)).astype(int) + mapshape[1]/2
            yp_bs = (np.round((dec_interp-dec_src)*60./reso_arcmin)).astype(int) + mapshape[0]/2
            xp_all = np.zeros([len(gnames2),len(xp_bs)],dtype=int) 
            yp_all = np.zeros([len(gnames2),len(xp_bs)],dtype=int)
            for i in np.arange(len(gnames2)):
                xp_all[i,:] = xp_bs + bp[gnames2[i]].x_offset/core.G3Units.arcmin/reso_arcmin
                yp_all[i,:] = yp_bs + bp[gnames2[i]].y_offset/core.G3Units.arcmin/reso_arcmin
            for i in np.arange(len(gnames2)):
                for j in np.arange(len(xp_bs)):
                    try:
                        map_weighted[yp_all[i,j],xp_all[i,j]] += bdata2[i,j]
                        nhits[yp_all[i,j],xp_all[i,j]] += 1.
                    except:
                        pass
