import numpy as np
from spt3g import core, std_processing, gcp, coordinateutils
from spt3g.util import tctools
import glob, os, pickle
from astropy.io import ascii

in_datestamp = '01feb19'
out_datestamp = '01feb19'
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
    if dtc['ampls'][0,i] > 0. and dtc['ampls'][0,i] < 1. and dtc['x_off_arcsec'][0,i] > -200. and dtc['x_off_arcsec'][0,i] < 200. and dtc['y_off_arcsec'][0,i] > -200. and dtc['y_off_arcsec'][0,i] < 200. and dtc['fwhms'][0,i] > 1.4 and dtc['fwhms'][0,i] < 2.5:
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
        framestarts = []
        nptot = 0
        for bdf in bdfiles:
            f1 = core.G3File(bdf)
            for frame in f1:
                if frame.type is core.G3FrameType.Scan:
                    if 'Turnaround' not in frame:
                        framestarts.append(nptot)
                        nptot += len(frame['OnlineBoresightRa'])
                        ratemp2.append(frame['OnlineBoresightRa'])
                        dectemp2.append(frame['OnlineBoresightDec'])
                        thisdict['collimation'] = np.asarray(frame['OnlinePointingModel']['fixedCollimation'])/core.G3Units.deg
                        thisdict['flexure'] = np.asarray(frame['OnlinePointingModel']['flexure'])/core.G3Units.deg
                        stime1 = core.G3Time(frame['RawTimestreams_I'][list(frame['RawTimestreams_I'].keys())[0]].stop.time)
        nratemp2 = tctools.list_to_array(ratemp2)/core.G3Units.deg
        ndectemp2 = tctools.list_to_array(dectemp2)/core.G3Units.deg
        dxdectemp2 = np.abs((nratemp2-thisdict['ra_src'])*np.cos(ndectemp2*core.G3Units.deg))
        ddectemp2 = np.abs(ndectemp2-thisdict['dec_src'])
        disttemp = np.sqrt(dxdectemp2**2 + ddectemp2**2)
        indtemp = np.argmin(disttemp)
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
                        print(framecounter)
                        if framecounter == goodframe:
                            print("getting trans")
                            coordinateutils.coordsysmodules.FillCoordTransRotations(frame,transform_store_key='AzElRot',end_coord_sys = core.MapCoordReference.Local) 
                            trans1 = frame['OnlineRaDecRotation']
                            trans2 = frame['AzElRot']
                            tempazcheck = frame['RawBoresightAz']
                            tempelcheck = frame['RawBoresightEl']
                        framecounter += 1

        radec_q_0 = coordinateutils.ang_to_quat(thisdict['ra_src']*core.G3Units.deg,thisdict['dec_src']*core.G3Units.deg)
        xy_q_0 = trans1**-1 * radec_q_0 * trans1 
        aznel_q_0 = trans2 * xy_q_0 * trans2**-1
        az_0, nel_0 = coordinateutils.quat_to_ang(aznel_q_0)
        thisdict['az_src'] = az_0[indinframe]/core.G3Units.deg
        thisdict['el_src'] = -nel_0[indinframe]/core.G3Units.deg
        if thisdict['az_src'] < 0:
            thisdict['az_src'] += 360.
        if thisdict['el_src'] < 0:
            thisdict['el_src'] += 360.

        radec_q_1 = coordinateutils.ang_to_quat(thisdict['ra_meas']*core.G3Units.deg,thisdict['dec_meas']*core.G3Units.deg)
        xy_q_1 = trans1**-1 * radec_q_1 * trans1 
        aznel_q_1 = trans2 * xy_q_1 * trans2**-1
        az_1, nel_1 = coordinateutils.quat_to_ang(aznel_q_1)
        thisdict['az_meas'] = az_1[indinframe]/core.G3Units.deg
        thisdict['el_meas'] = -nel_1[indinframe]/core.G3Units.deg
        if thisdict['az_meas'] < 0:
            thisdict['az_meas'] += 360.
        if thisdict['el_meas'] < 0:
            thisdict['el_meas'] += 360.

#        print(notavariable)

    else:
        thisdict['xoff'] = np.nan
        thisdict['yoff'] = np.nan
        thisdict['dxoff'] = np.nan
        thisdict['dyoff'] = np.nan
        thisdict['az_src'] = np.nan
        thisdict['el_src'] = np.nan
        thisdict['az_meas'] = np.nan
        thisdict['el_meas'] = np.nan

    thisdict['dra'] = thisdict['dxoff']/np.cos(thisdict['dec_src']*np.pi/180.)
    thisdict['ddec'] = thisdict['dyoff']
    thisdict['mjd'] = (time1.mjd+stime1.mjd)/2.
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
f.write('# Name, R.A. [deg], dec. [deg], MJD, Source az [deg], Source el [deg], Measured Az [deg], dAz [arcmin], Measured El [deg], dEl [arcmin], tilt0 [deg], tilt1 [deg], tilt2 [deg]\n')

timestrs = list(obsdict.keys())
timestrs.sort()

for timestr1 in timestrs:
    od = obsdict[timestr1]
    if od['had_good_fit']:
        rastr = '%08.4f' % od['ra_src']
        decstr = '%08.4f' % od['dec_src']
        mjdstr = '%10.4f' % od['mjd']
        srcazstr = '%08.4f' % od['az_src']
        srcelstr = '%07.4f' % od['el_src']
        measazstr = '%08.4f' % od['az_meas']
        measelstr = '%07.4f' % od['el_meas']
        dazstr = '%6.3f' % od['dra']
        delstr = '%6.3f' % od['ddec']
        tilt0str = '%9.6f' % od['tilts'][0]
        tilt1str = '%9.6f' % od['tilts'][1]
        tilt2str = '%9.6f' % od['tilts'][2]
        
        f.write(od['source']+delim+rastr+delim+decstr+delim+mjdstr+delim+srcazstr+delim+srcelstr+delim+measazstr+delim+dazstr+delim+measelstr+delim+delstr+delim+tilt0str+delim+tilt1str+delim+tilt2str+delim+'\n')

f.close()
