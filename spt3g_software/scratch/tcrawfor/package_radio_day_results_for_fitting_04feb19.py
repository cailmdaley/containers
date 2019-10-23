import numpy as np
from spt3g import core, std_processing, gcp, coordinateutils
from spt3g.util import tctools
import glob, os, pickle
from astropy.io import ascii

def grab_boards(f, boardlist = []):
#    keylist = ['TrackerPointing',
#               'ACUStatus',
#               'SourceName',
#               'array',
#               'TrackerStatus',
#               'GCPFeatureBits',
    key = 'antenna0'
#    key = 'TrackerPointing'
    try:
        boardlist.append(f[key])
    except:
        pass

in_datestamp = '04feb19_newboloprops_nopm'
out_datestamp = '04feb19_newboloprops_nopm'
userdir = os.path.expanduser('~/')
tcfile = userdir + '/fit_radio_day_maps_'+in_datestamp+'.pkl'
dtc = pickle.load(open(tcfile,'rb'))

f = open('/home/tcrawfor/spt_code/gcp/config/ephem/source.cat')
srclines = f.readlines()
f.close()
srclines = srclines[1978:]

obsdict = {}

for i in np.arange(len(dtc['obsids'])):
#for i in np.arange(len(dtc['obsids'])) + 19:

    obsid1 = dtc['obsids'][i]
    src1 = dtc['srcnames'][i]
    time1 = std_processing.obsid_to_g3time(obsid1)

    print(src1)
    print(time1.Summary())

    thisdict = {}
    src1_up = src1.upper()
    thisdict['source'] = src1_up
    for line in srclines:
        if src1_up in line and 'vlbi' not in line:
            thisdict['ra_src'] = tctools.ten_tc(line[19:32],h2deg=True)
            thisdict['dec_src'] = tctools.ten_tc(line[36:50])

    thisdict['had_good_fit'] = False
    xofftemp = []
    yofftemp = []
    dofftemp = []
    dyofftemp = []
    ratemp = []
    dectemp = []
    if dtc['ampls'][0,i] > 0. and dtc['ampls'][0,i] < 1. and dtc['x_off_arcsec'][0,i] > -100. and dtc['x_off_arcsec'][0,i] < 1000. and dtc['y_off_arcsec'][0,i] > 400. and dtc['y_off_arcsec'][0,i] < 1000. and dtc['fwhms'][0,i] > 1.4 and dtc['fwhms'][0,i] < 2.5:
        thisdict['had_good_fit'] = True
        xofftemp.append(dtc['x_off_arcsec'][0,i])
        yofftemp.append(dtc['y_off_arcsec'][0,i])
        dofftemp.append(dtc['pos_unc_arcsec'][0,i])
        ratemp.append(dtc['ra_src_meas'][0,i])
        dectemp.append(dtc['dec_src_meas'][0,i])
        if dtc['ampls'][1,i] > 0. and dtc['ampls'][1,i] < 1. and np.abs(dtc['x_off_arcsec'][1,i]-dtc['x_off_arcsec'][0,i]) < 15. and np.abs(dtc['y_off_arcsec'][1,i]-dtc['y_off_arcsec'][0,i]) < 15. and dtc['fwhms'][1,i] > 1.1 and dtc['fwhms'][1,i] < 2.:
            xofftemp.append(dtc['x_off_arcsec'][1,i])
            yofftemp.append(dtc['y_off_arcsec'][1,i])
            dofftemp.append(dtc['pos_unc_arcsec'][1,i])
            ratemp.append(dtc['ra_src_meas'][1,i])
            dectemp.append(dtc['dec_src_meas'][1,i])
        if dtc['ampls'][2,i] > 0. and dtc['ampls'][2,i] < 1. and np.abs(dtc['x_off_arcsec'][2,i]-dtc['x_off_arcsec'][0,i]) < 15. and np.abs(dtc['y_off_arcsec'][2,i]-dtc['y_off_arcsec'][0,i]) < 15. and dtc['fwhms'][2,i] > 0.9 and dtc['fwhms'][2,i] < 1.4:
            xofftemp.append(dtc['x_off_arcsec'][2,i])
            yofftemp.append(dtc['y_off_arcsec'][2,i])
            dofftemp.append(dtc['pos_unc_arcsec'][2,i])
            ratemp.append(dtc['ra_src_meas'][2,i])
            dectemp.append(dtc['dec_src_meas'][2,i])
        wtemp = 1./np.asarray(dofftemp)**2
        thisdict['xoff'] = np.sum(np.asarray(xofftemp)*wtemp)/np.sum(wtemp)/60.
        thisdict['dxoff'] = np.sqrt(1./np.sum(wtemp))/60.
        thisdict['yoff'] = np.sum(np.asarray(yofftemp)*wtemp)/np.sum(wtemp)/60.
        thisdict['dyoff'] = thisdict['dxoff']
        thisdict['ra_meas'] = np.sum(np.asarray(ratemp)*wtemp)/np.sum(wtemp)
        thisdict['dec_meas'] = np.sum(np.asarray(dectemp)*wtemp)/np.sum(wtemp)

        bdfiles = glob.glob('/spt_data/bolodata/fullrate/'+src1_up+'/'+str(obsid1)+'/0*.g3')
        bdfiles.sort()
        if (len(bdfiles)//2)*2 == len(bdfiles):
            bdfiles = bdfiles[len(bdfiles)//2-1:len(bdfiles)//2+1]
        else:
            bdfiles = [bdfiles[len(bdfiles)//2]]
        ratemp2 = []
        dectemp2 = []
        mjdtemp = []
        framestarts = []
        nptot = 0
        for bdf in bdfiles:
            f1 = core.G3File(bdf)
            for frame in f1:
                if frame.type is core.G3FrameType.Scan:
                    if 'Turnaround' not in frame:
                        framestarts.append(nptot)
                        thisra = frame['OnlineBoresightRa']
                        nptot += len(thisra)
                        ratemp2.append(thisra)
                        dectemp2.append(frame['OnlineBoresightDec'])
                        mjds = np.arange(len(thisra))/np.float(len(thisra)-1)*(thisra.stop.mjd-thisra.start.mjd) + thisra.start.mjd
                        mjdtemp.append(mjds)
                        thisdict['collimation'] = np.asarray(frame['OnlinePointingModel']['fixedCollimation'])/core.G3Units.deg
                        thisdict['flexure'] = np.asarray(frame['OnlinePointingModel']['flexure'])/core.G3Units.deg
        nratemp2 = tctools.list_to_array(ratemp2)/core.G3Units.deg
        ndectemp2 = tctools.list_to_array(dectemp2)/core.G3Units.deg
        nmjds = tctools.list_to_array(mjdtemp)
        dxdectemp2 = np.abs((nratemp2-thisdict['ra_src'])*np.cos(ndectemp2*core.G3Units.deg))
        ddectemp2 = np.abs(ndectemp2-thisdict['dec_src'])
        disttemp = np.sqrt(dxdectemp2**2 + ddectemp2**2)
        indtemp = np.argmin(disttemp)
        mjd = nmjds[indtemp]
        framestarts = np.asarray(framestarts)
        whind= np.where(framestarts < indtemp)
        if len(whind[0]) == 0:
            goodframe = len(framestarts) - 1
        else:
            goodframe = np.max(whind[0])
        indinframe = indtemp - framestarts[goodframe]
        framecounter = 0
        for bdf in bdfiles:
            f1 = core.G3File(bdf)
            for frame in f1:
                if frame.type is core.G3FrameType.Scan:
                    if 'Turnaround' not in frame:
#                        print(framecounter)
                        if framecounter == goodframe:
#                            print("getting trans")
                            coordinateutils.coordsysmodules.FillCoordTransRotations(frame,transform_store_key='AzElRot',end_coord_sys = core.MapCoordReference.Local) 
                            trans1 = frame['OnlineRaDecRotation']
                            trans2 = frame['AzElRot']
                            tempazcheck = frame['RawBoresightAz']
                            tempelcheck = frame['RawBoresightEl']
                            tempazcheck2 = frame['OnlineBoresightAz']
                            tempelcheck2 = frame['OnlineBoresightEl']
                            break
                        framecounter += 1
        data1 = []
        pipe1 = core.G3Pipeline()
        pipe1.Add(std_processing.ARCTimerangeReader, start_time=tempazcheck.start, stop_time=tempazcheck.stop, basedir='/spt_data/arc/')
        pipe1.Add(grab_boards, boardlist=data1)
        pipe1.Run()
        aztemp = []
        eltemp = []
        for tframe in data1:
            aztemp.append(tframe['tracker']['horiz_topo'][0])
            eltemp.append(tframe['tracker']['horiz_topo'][1])
        naztemp = tctools.list_to_array(aztemp)/3.6e6
        neltemp = tctools.list_to_array(eltemp)/3.6e6
        naztemp2 = np.interp(np.arange(len(tempazcheck))/np.float(len(tempazcheck)-1),np.arange(len(naztemp))/np.float(len(naztemp)-1),naztemp)
        neltemp2 = np.interp(np.arange(len(tempelcheck))/np.float(len(tempelcheck)-1),np.arange(len(neltemp))/np.float(len(neltemp)-1),neltemp)
        thisdict['az_src'] = naztemp2[indinframe]
        thisdict['el_src'] = neltemp2[indinframe]
        thisdict['mjd'] = mjd
    else:
        thisdict['xoff'] = np.nan
        thisdict['yoff'] = np.nan
        thisdict['dxoff'] = np.nan
        thisdict['dyoff'] = np.nan
        thisdict['az_src'] = np.nan
        thisdict['el_src'] = np.nan
        thisdict['ra_meas'] = np.nan
        thisdict['dec_meas'] = np.nan
        thisdict['mjd'] = np.nan

    thisdict['ra_off'] = (thisdict['ra_meas'] - thisdict['ra_src'])*60.
    thisdict['dec_off'] = (thisdict['dec_meas'] - thisdict['dec_src'])*60.
    thisdict['dra'] = thisdict['dxoff']/np.cos(thisdict['dec_src']*np.pi/180.)
    thisdict['ddec'] = thisdict['dyoff']
    thisdict['tilts'] = np.asarray(frame['OnlinePointingModel']['tilts']/core.G3Units.deg)
    thisdict['flexure'] = np.asarray(frame['OnlinePointingModel']['flexure']/core.G3Units.deg)
    thisdict['fixedCollimation'] = np.asarray(frame['OnlinePointingModel']['fixedCollimation']/core.G3Units.deg)

    obsdict[time1.Summary()] = thisdict

delim = ', '
f = open('pointing_model_fit_data/summer2019/radio_day_summary_for_fitting_'+out_datestamp+'.csv','w')
f.write('#Written by package_radio_day_results_for_fitting_04feb19.py \n')
f.write('#Using fit results in ' + tcfile + '\n')
f.write('#tilts 0.006933888889, -0.007933888889, -0.022062777778 \n') # note these values will be overridden by per-observation tilts (just here for formatting)
f.write('#flexure radio, ' + '%10.8f' % thisdict['flexure'][0] + ', ' + '%10.8f' % thisdict['flexure'][1] + '\n')
f.write('#collimate_fixed radio, ' + '%10.8f' % thisdict['fixedCollimation'][0] + ', ' + '%10.8f' % thisdict['fixedCollimation'][1] + '\n')
f.write('# Name, MJD, R.A. [deg], Dec. [deg], R.A. offset [arcmin], Dec. offset [arcmin], dR.A. [arcmin], dDec. [arcmin], Source az [deg], Source el [deg], Tilt0 [deg], Tilt1 [deg], Tilt2 [deg]\n')

timestrs = list(obsdict.keys())
timestrs.sort()

for timestr1 in timestrs:
    od = obsdict[timestr1]
    if od['had_good_fit']:
        mjdstr = '%10.4f' % od['mjd']
        rastr = '%08.4f' % od['ra_src']
        decstr = '%08.4f' % od['dec_src']
        raoffstr = '%06.4f' % od['ra_off']
        decoffstr = '%06.4f' % od['dec_off']
        drastr = '%6.4f' % od['dra']
        ddecstr = '%6.4f' % od['ddec']
        srcazstr = '%08.4f' % od['az_src']
        srcelstr = '%07.4f' % od['el_src']
        tilt0str = '%9.6f' % od['tilts'][0]
        tilt1str = '%9.6f' % od['tilts'][1]
        tilt2str = '%9.6f' % od['tilts'][2]
        
        f.write(od['source']+delim+mjdstr+delim+rastr+delim+decstr+delim+raoffstr+delim+decoffstr+delim+drastr+delim+ddecstr+delim+srcazstr+delim+srcelstr+delim+tilt0str+delim+tilt1str+delim+tilt2str+'\n')

        print(od['source']+delim+mjdstr+delim+rastr+delim+decstr+delim+raoffstr+delim+decoffstr+delim+drastr+delim+ddecstr+delim+srcazstr+delim+srcelstr+delim+tilt0str+delim+tilt1str+delim+tilt2str+'\n')

f.close()
