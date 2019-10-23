import numpy as np
import os
import copy
import glob
from spt3g import core, coordinateutils
from spt3g.pointing import pointing_tools as pt

'''
   This module contains functions for obtaining a list of pointing parameters 
   to pass to the offline and/or online pointing model.
'''

def grabAzTiltParams(input_files, tilts_dir='/poleanalysis/sptdaq/azTilts'):
    '''
    Extract tilt parameters from the nearest azimuth tilt obsID to this
    observations's obsID.

    Arguments
    ---------
    input_files : array-like
        A list or array of absolute paths (strings) to G3 files.
    tils_dir : string, optional
        the absilute path to a directory containing azimuth tilt
        observation data.

    Returns
    -------
    a2, a3 : tuple, 2 elements
        A 2-element tuple comprised of a2 (float) and a3 (float), which
        are tilt values mapped to pointing model parameters
    '''
    a2 = None
    a3 = None
    for fname in input_files:
        for frame in core.G3File(fname):
            #Extract the obsID of the frame.
            if frame.type == core.G3FrameType.Observation:
                obsID = frame['ObservationID']

            #Now track down which tilt obsID should be matched to this.
            dirs = glob.glob(os.path.join(tilts_dir,'*'))
            obslist = [int(f.split('/')[-1]) for f in dirs]
            obslist.sort()
            thisID = np.where(obslist >= obsID)[0][0]
            tilt_obsID1 = obslist[thisID]
            tilt_obsID2 = obslist[thisID+1]

            #Which tilt_obsID is closest to obsID?
            tilt_obs = np.array([tilt_obsID1, tilt_obsID2])
            diffs = obsID - tilt_obs

            index = np.where(diffs == np.min(diffs))[0]

            tilt_obsID = tilt_obs[index]

            #Now load tilt information from the relevant g3 file.
            info = [f for f in core.G3File(os.path.join(tilts_dir,str(tilt_obsID[0]), '0001.g3'))][0]
            tiltHA = info['tiltHA']/core.G3Units.degrees
            tiltLat = info['tiltLat']/core.G3Units.degrees

            #Now map the tilt values to pointing model parameters.
            a2 = 0.88*tiltHA
            a3 = -0.88*tiltLat
            break
        
        #Break the for loop once we have a2,a3
        if a2 != None:
            break

    return a2, a3


#Online Az model given raw Az/El and pointing parameters) 
def CorrectedAz(az, el, a2, a3, a4, a5, az0, DET=0,
                flags=['az_tilts', 'el_tilts',
                       'flexure', 'collimation',
                       'refraction']):
    """
    Return an azimuth value that has been corrected by the pointing
    model.

    Arguments
    ---------
    az : float
        An azimuth position in units of degrees.
    el : float
        An elevation position in units of degrees.
    a2 : float
        Description
    a3 : float
        Description
    a4 : float
        Description
    a5 : float
        Description
    az0 : float
        Description
    DET : number, optional
        Description
    flags : array-like
        An array of flags describing which modifications to make to the
        input azimuth value. The accepted flags are:
        'az_tilts', 'el_tilts', 'collimation', 'thermolin', 'flexure',
        'collimation', 'refraction'
        The default behavior is to use the flags:
        'az_tilts', 'el_tilts', 'flexure', 'collimation', 'refraction'

    Returns
    -------
    az + d_az : float
        The corrected azimuth position as calculated by the azimuth model.
    """
    d_az = -az0

    if 'az_tilts' in flags:
        d_az += (a2* np.cos(az/core.G3Units.rad) +\
                   a3 *np.sin(az/core.G3Units.rad))*np.tan(el/core.G3Units.rad)

    if 'el_tilts' in flags:
        d_az += a4*np.tan(el/core.G3Units.rad)

    if 'collimation' in flags:
        d_az += - a5/np.cos(el/core.G3Units.rad)


    if 'thermolin' in flags:
        d_az += -DET*np.tan(el/core.G3Units.rad)
    
    return az + d_az


#Online El model given raw Az/El and pointing parameters, and refraction)
def CorrectedEl(az, el, a0, a1, a2, a3, a6, refraction,
                DEL = 0.0,
                flags=['az_tilts', 'el_tilts', 'flexure',
                       'collimation', 'refraction']):
    """
    Return an elevation value that has been corrected by the pointing
    model.

    Arguments
    ---------
    az : float
        An azimuth position in units of degrees.
    el : float
        An elevation position in units of degrees.
    a0 : float
        Description
    a1 : float
        Description
    a2 : float
        Description
    a3 : float
        Description
    a6 : float
        Description
    refraction : float
    DET : number, optional
        Description
    flags : array-like
        An array of flags describing which modifications to make to the
        input azimuth value. The accepted flags are:
        'az_tilts', 'el_tilts', 'collimation', 'thermolin', 'flexure',
        'collimation', 'refraction'
        The default behavior is to use the flags:
        'az_tilts', 'el_tilts', 'flexure', 'collimation', 'refraction'

    Returns
    -------
    el + d_el : float
        The corrected elevation position as calculated by the elevation
        model.
    """

    d_el = 0.

    if 'flexure' in flags:
        d_el += a0*np.sin(el/core.G3Units.rad) + a1*np.cos(el/core.G3Units.rad)

    if 'az_tilts' in flags:
        d_el += -(a2*np.sin(az/core.G3Units.rad) - a3*np.cos(az/core.G3Units.rad))

    if 'collimation' in flags:
        d_el += -a6

    if 'refraction' in flags:
        d_el += -refraction

    if 'thermolin' in flags:
        d_el += -DEL

    return el + d_el

@core.indexmod
class CorrectBadRefraction(object):
    '''
    Until Jan 26 2018, we were using the wrong refraction register
    to correct elevation.  This module goes back to the arcfile(s),
    and takes the median refraction correction for the observation.
    This is only valid for shorter observations that cover a small 
    elevation range.

    The temperature was also wrong during this time, so we use that to
    detect bad data.
    '''
    def __init__(self, bad_el_key = 'OnlineBoresightEl', dec = False,
                 obstime = 30 * core.G3Units.minutes, 
                 arcdir = '/spt_data/arc'):
        '''
        bad_el_key: str ('OnlineBoresightEl')
            The Elevation-like key to be fixed.  The fixed key is replaced with
            fixed data.
        dec: bool (False)
            If `bad_el_key` points to a declination key, set to True.
            This module assumes that el = -dec, so expect small errors.
        obstime: float (30 * core.G3Units.minutes)
            Approximate length of the observation.  Longer is better,
            but this only controls how much arcfile data we read to get
            the true refraction correction.  If this matters, you need to 
            do all of this properly, rather than the approximations we're using.
        arcdir: str '/spt_data/arc'
            Directory containing arcfiles
        '''
        # put the import here so our dependence is minimal
        from spt3g.util.extractdata import extract_keys
        self.extract_keys = extract_keys
        self.bad_el_key = bad_el_key
        self.dec = dec
        self.median_refraction = None
        self.obstime = obstime
        self.arcdir = arcdir
        
    def __call__(self, fr):
        def t_from_arc(path):
            return core.G3Time(os.path.basename(path).split('.')[0])

        if 'ObservationStart' in fr:
            self.start = fr['ObservationStart']
            self.stop = self.start + self.obstime
            return
        if fr.type != core.G3FrameType.Scan:
            return
        T = np.median(fr['TrackerPointing'].telescope_temp) / core.G3Units.K
        old_refrac = np.median(fr['TrackerPointing'].refraction)
        if T > 0 or abs(old_refrac / core.G3Units.deg - 0.012463) > 1e-3:
            # Magic number for something close to the value of the thing
            # we accidentally stored in TrackerPointing.refraction.
            # The second check hopefully makes this robust to weird
            # temperature readings in otherwise good data.
            # We probably have good data that doesn't need correcting
            return
            
        # get the data from the arcfile(s), if we haven't yet
        if self.median_refraction is None:
            from spt3g.std_processing import FindARCFiles
            arcfiles = [FindARCFiles.ARCForTime(self.start, self.arcdir)]
            # remove this line
            # grab arcfiles until we get to the nominal end of the ob
            while t_from_arc(arcfiles[-1]).time + 1000 * core.G3Units.s < \
                  self.stop.time:
                arcfiles.append(FindARCFiles.NextARCFile(arcfiles[-1]))
                keys = {'refraction': ['antenna0', 'tracker', 'refraction', 2],
                        'time': ['antenna0', 'tracker', 'utc']}
                arcdata = self.extract_keys(arcfiles, keys)
                # Now, grab all the times that we're observing the source
                inds = np.where((arcdata['time'] < self.stop) & 
                                (arcdata['time'] > self.start))[0]
                self.median_refraction = np.median(arcdata['refraction'][inds])
            
        core.log_warn('Fixing bad refraction correction.  If this data is from Jan 26, 2018 or later, something else has probably gone wrong.')
        if self.dec:
            raw_el = -1 * fr.pop(self.bad_el_key, None) + old_refrac 
            fr[self.bad_el_key] = -1 * (raw_el - self.median_refraction)
        else:
            raw_el = fr.pop(self.bad_el_key, None) + old_refrac
            fr[self.bad_el_key] = raw_el - self.median_refraction

@core.indexmod
class CorrectBoresightPointing(object):
    '''
    Makes two timesteams (<output>Az, <output>El) that apply offsets
    from a given pointing model to environmental data to the specified
    raw encoder timestreams. Interpolates over pointing dropouts by default.

    The usual model choices are 'OnlinePointingModel' and 'OfflinePointingModel'
    and can reference dictionaries of parameters in any frame type.

    flags: list of pointing model segments to turn on.  Default is everything
           in online model (no thermolin corrections).

    ***WARNING***
    CalcAzElToRaDecOffsetFromBoresight
    depends on this module working on coordinates that are slightly offset from the boresight.
    If that functionality changes you will need to modify that or all of the pointing will
    be wrong.

    '''
    def __init__(self,
                 raw_az_key = 'RawBoresightAz', raw_el_key = 'RawBoresightEl',
                 output = 'OnlineBoresight', model = 'OnlinePointingModel',
                 flags = ['az_tilts', 'el_tilts', 'flexure', 'collimation',
                         'refraction']):
        # XXX: contrary to documentation, does no interpolation
        self.raw_az_key = raw_az_key
        self.raw_el_key = raw_el_key
        self.output = output
        self.model_key = model
        self.flags = flags

        self.model = None # Will be filled later

    def __call__(self, frame):
        # Try to cache the model, wherever it appears
        if self.model_key in frame:
            self.model = frame[self.model_key]

        # Otherwise, ignore non-scan frames
        if frame.type != core.G3FrameType.Scan:
            return

        # Now get model params
        a0, a1 = self.model['flexure'][0:2]
        a2, a3, a4 = self.model['tilts'][0:3]
        a5, a6 = self.model['fixedCollimation'][0:2]
        az0 = frame['TrackerPointing'].encoder_off_x[0]
        refraction = np.median(frame['TrackerPointing'].refraction)

        p = {'a0':a0, 'a1':a1, 'a2':a2, 'a3':a3, 'a4':a4,
             'a5':a5, 'a6':a6, 'az0':az0, 'refraction':refraction}

        # And apply corrections
        if 'thermolin' in self.flags:
            l1 =  np.median(frame['TrackerPointing'].linsens_avg_l1)
            l2 =  np.median(frame['TrackerPointing'].linsens_avg_l2)
            r1 =  np.median(frame['TrackerPointing'].linsens_avg_r1)
            r2 =  np.median(frame['TrackerPointing'].linsens_avg_r2)

            lin_data = get_lin_sens(l1, l2, r1, r2)
            p_thermolin = thermo2pointing(frame['TrackerPointing'].scu_temp, 
                                          lin_data)

            p['DET'] = p_thermolin['DET']
            p['DEL'] = p_thermolin['DEL']
        else:
            p['DET'] = 0.0
            p['DEL'] = 0.0

        corrected_az = CorrectedAz(frame[self.raw_az_key],
				   frame[self.raw_el_key],
                                   p['a2'], p['a3'], p['a4'], 
                                   p['a5'], p['az0'], p['DET'], 
                                   flags=self.flags)
        corrected_el = CorrectedEl(frame[self.raw_az_key],
				   frame[self.raw_el_key],
                                   p['a0'], p['a1'], p['a2'], 
                                   p['a3'], p['a6'], p['refraction'], p['DEL'],
                                   flags=self.flags)

        frame[self.output + 'Az'] = corrected_az
        frame[self.output + 'El'] = corrected_el

def TestThermolin(frame):
    '''
    Gather DET and DEL corretions for a given observation.
    '''
    if frame.type == core.G3FrameType.Scan:
        p = extractOfflinePointingParameters(frame)
        l1 =  np.median(frame['TrackerPointing'].linsens_avg_l1)
        l2 =  np.median(frame['TrackerPointing'].linsens_avg_l2)
        r1 =  np.median(frame['TrackerPointing'].linsens_avg_r1)
        r2 =  np.median(frame['TrackerPointing'].linsens_avg_r2)

        lin_data = get_lin_sens(l1, l2, r1, r2)
        p_thermolin = thermo2pointing(frame['TrackerPointing'].scu_temp, 
                                  lin_data)

        p['DET'] = p_thermolin['DET']
        p['DEL'] = p_thermolin['DEL']

        frame['DET'] = core.G3Double(p['DET'])
        frame['DEL'] = core.G3Double(p['DEL'])
        


#----------------------------------------------------------------------------------------
#Below is code used to calculate collimation corrections based on scu temperature sensor
#and yoke arm metrology readings.
#This may or may not get replaced with the thermolin code RK wrote for use with EHT.
#----------------------------------------------------------------------------------------
def get_lin_sens(l1, l2, r1, r2):
    """
    Translated from the IDL get_lin_sens.pro written by RK.

    The purpose of this function is to return the linear sensor and thermometry sensor data
    from a given time window.  Linearly interpolate over dropouts in the 
    sensor data.
    
    INPUTS
        date: Array of dates.

    OUTPUTS
        S: a dictionary with the following substructures:
            'utc': The 100 Hz UTC.
            'l1': The 100 Hz L1 length, in mm.
            'l2': The 100 Hz L2 length, in mm.
            'r1': The 100 Hz R1 length, in mm.
            'r2': The 100 Hz R2 length, in mm.
            'del': The 100 Hz elevation correction, in arcseconds.
            'daz': The 100 Hz azimuth correction, in arcseconds.
            'det': The 100 Hz elevation tilt correction, in arcseconds.
            'temp': The thermometry data, which is an array with [nthermos, nsamples].

     Translated: October 2012, JWH.
     Originally Written: April 2008, RK.
     Modifications: Take linear sensor data in place of the 'date' input. 7 Dec 2012, SH
    """

    #Yoke dimensions in mm.
    Rs = 1652.
    Rh = 3556.
    Ry = 6782.

    #Calculate corrections in arcsec.
    DEL = (1./(2.*Rs))*(l2 - l1 + r2 - r1)*(3600.*180./np.pi)
    DAZ = (Rh/(Ry*Rs))*(l1 - l2 - r1 + r2)*(3600.*180./np.pi)
    DET = (1./(2.*Ry))*(r1 + r2 - l1 - l2)*(3600.*180./np.pi)

    #Subtract medians calculated from RCW38 observations.
    DAZ -= 0.0 #38.7
    DEL -= 0.0 #27.6
    DET -= 0.0 #18.6

    #Fill the output dictionary.
    s = {#'utc':utc, 
         'l1':l1, 'l2':l2, 'r1':r1, 'r2':r2,
         'del':DEL, 'daz':DAZ, 'det':DET}

    return s


def thermo2pointing(scu_temp, this_lin, 
                    thermometry_config_file='thermometer_pointing_coefficients',
                    nointerp=True):
    """
    The purpose of this function is to provide pointing corrections DET and DEL, given an
    input array of structure thermometry and/or linear sensor data.  The model is just linear 
    in the thermometry + linear sensors.  The coefficients for the model are stored in an 
    external common txt file.

    INPUTS:
        scu_temp_in: the array of thermometry + linear sensor data.  It has
                     dimensions of (63, nsamples) where there are 60 thermometers
                     and 3 linear sensors.  The thermometry should be raw (degrees C).

    OUTPUTS:
        s - a dictionary with the following fields:
            DET: the DET (elevation axis tilt) correction, in arcseconds.
            DEL: the DEL (plain old elevation) correction, in arcseconds.


    EXCEPTIONS
        ValueError if the config file doesn't match the size of the scu_temp register data. 

    Translated to python from the original IDL written by RK, Jan 2009.
    Translated by JWH October 2012.
    """
    scu_temp = np.median(np.reshape(scu_temp, (int(len(scu_temp)/60), 60)), axis=0)
    scu_temp = scu_temp/core.G3Units.K - 273.15 # convert to C

    scu_temp = np.hstack([scu_temp, this_lin['daz'], this_lin['del'], this_lin['det']])

    scu_temp = np.array([scu_temp]).T

    # Read in a config file that contains the coefficients for going from
    # structure temperatures to pointing offsets DET and DEL.
    index, det_coeff, del_coeff, neighbors1, neighbors2 = readThermoConfig()

    nsamples = scu_temp.shape[1]
    nthermo = scu_temp.shape[0]
    if nthermo != len(det_coeff)-1 or nthermo != len(del_coeff)-1:
        raise ValueError('# of thermometers in scu_temp does not match # of coefficients')

    #Interpolate over dropouts
    npts = nsamples
    thermo_zero = -200.0
    if nointerp==True and nsamples > 1:
        for i in range(nthermo-3):
            whnodrop1 = np.nonzero((scu_temp[i] != thermo_zero) &
                                   (scu_temp[i] != 0.0) &
                                   (scu_temp[i] > -150.) &
                                   (scu_temp[i] < 40.))[0]
            nnodrop = len(whnodrop1)
            if nnodrop < npts/2.:
                continue
            thisdata = scu_temp[i]
            thisdata = pt.interp_over_dropouts(thisdata, whnodrop=whnodrop1)
            scu_temp[i] = thisdata

    #The thermometer indexed in IDL by i=40 begain to have problems in 2011.
    #The cause of these problems are unknown, but basically the temperatures
    #that are recorded are crazy.  This can screw up the pointing, since the 
    #offline pointing model depends on the telescope temperatures.  So if it's
    #2011 or later, and the i=40 thermometer looks crazy, let's replace its data
    #with that of a nearby thermometer, i=42.
    if (np.abs(np.median(scu_temp[40]) - np.median(scu_temp[42])) > 10.):
        scu_temp[40] = scu_temp[42]

    # Look through each thermometer, checking to see if there are any
    # which are still equal to the "thermometer zero" value, typically
    # -200.0 C, which were not interpolated over because there's no good
    # data to use for the interpolation.  For these thermometers we want
    # to replace their output with that of their "neighbors", where 
    # neighbor is defined as a thermometer that historically had a 
    # similar temperature.
    scu_tempo = scu_temp.copy()
    for i in range(nthermo-3):
        this_temp = np.array(scu_temp[i]).reshape(len(scu_temp[i]))
        wh_zero = np.nonzero(this_temp == thermo_zero)[0]
        n_zero = len(wh_zero)
        if n_zero > 0:
            #This thermometer is returning the "zero" value.  Replace its
            #data from that from its closest possible neighbor with similar data.
            this_n1 = neighbors1[i]
            this_temp_n = scu_temp[int(this_n1)]
            wh_zero_n = np.nonzero(this_temp_n == thermo_zero)[0]
            n_zero_n = len(wh_zero_n)  
            if n_zero_n == 0:
                scu_temp[i] = this_temp_n
            else:
                this_n2 = neighbors2[i]
                this_temp_n = np.array(scu_temp[this_n2])
                wh_zero_n = np.nonzero(this_temp_n == thermo_zero)[0]
                n_zero_n = len(wh_zero_n)
                if n_zero_n == 0:
                    scu_temp[i] = this_temp_n
                else:
                    if (i<=25) or (i >= 40):
                        pass
                    wh_not_zero = np.nonzero(scu_tempo != thermo_zero)[0]
                    n_not_zero = len(wh_not_zero)
                    wh_not_zero = 0.
                    if n_not_zero > 0:
                        scu_temp[i] = ( np.zeros(len(scu_temp[i])) + 
                                        np.median(scu_tempo[scu_tempo != thermo_zero]) )

    #Convert temps to Kelvin.
    scu_temp[0:60,0] += 273.15

    #Restructure the coefficients
    det_dc = det_coeff[-1]
    del_dc = del_coeff[-1]
    det_coeff = np.matrix(det_coeff[:-1])
    del_coeff = np.matrix(del_coeff[:-1])
    scu_temp = np.matrix(scu_temp)

    #Subtract off the median (calculated from 2014 RCW38 obs) so we have mean zero corrections.
    DET = np.array(det_coeff*scu_temp + det_dc) - 0.0 #-20.2
    DEL = np.array(del_coeff*scu_temp + del_dc) - 0.0 #-7.0


    #Return corrections in units of degrees.
    return {'DET':DET[0][0]/3600., 'DEL':DEL[0][0]/3600.}


def readThermoConfig(config_file='thermometer_pointing_coefficients.txt'):
    '''
    Read in thermometer pointing coefficients calculated originally by RK for SPT-SZ.
    It's probably a good idea to update these coefficients...
    '''
    d = open(os.path.join(os.path.dirname(os.path.realpath(__file__)), config_file), 'r').read().split('\n')[:-1]
    
    index = []
    det_coeff = []
    del_coeff = []
    neighbors1 = []
    neighbors2 = []

    for i in range(len(d)):
        index.append(int(float(d[i].split(' ')[0])))
        det_coeff.append(float(d[i].split(' ')[1]))
        del_coeff.append(float(d[i].split(' ')[2]))
        neighbors1.append(int(float(d[i].split(' ')[3])))
        neighbors2.append(int(float(d[i].split(' ')[4])))

    return index, det_coeff, del_coeff, neighbors1, neighbors2

@core.usefulfunc
def OfflinePointingParamsAtTime(t, config_files):
    '''
    Read SPTpol-style configuration text files and interpolate the values
    into an appropriately-formatted G3MapVectorDouble that can be placed into a
    frame for use by CorrectBoresightPointing.

    t should be either a G3Time or a string that the G3Time constructor can
    interpret.
    '''

    params = {}
    for fname in config_files:
        with open(fname, 'r') as f:
            # Can't figure out a way to parse these with numpy.loadtxt :(
            for line in f.readlines():
                if len(line[0].strip()) == 0:
                    continue
                if line[0] == '#' and line[1] == '#':
                    headers = line[2:].split()
                elif line.split()[0] == '#mjd':
                    headers = line[1:].split()
                elif line[0] == '#':
                    continue
                elif line.startswith('VALIDITY'):
                    continue
                else:
                    fields = {headers[i]: float(j) for i, j in
                              enumerate(line.split())}
                    if fields['mjd'] not in params:
                        params[fields['mjd']] = {}
                    params[fields['mjd']].update(fields)

    keys = {k for v in params.values() for k in v if k != 'mjd'}

    if isinstance(t, str):
        t = core.G3Time(t)
    desiredmjd = t.mjd
    p_at_t = {}

    for k in keys:
        mjd = []
        vals = []
        for datum in params.values():
            if k not in datum:
                continue
            mjd.append(datum['mjd'])
            vals.append(datum[k])
        p_at_t[k] = np.interp([desiredmjd], mjd, vals)[0]

    out = core.G3MapVectorDouble()
    if 'a0' in p_at_t:
        out['flexure'] = core.G3VectorDouble([p_at_t['a0'], p_at_t['a1']])
    if 'a2' in p_at_t:
        out['tilts'] = core.G3VectorDouble([p_at_t['a2'], p_at_t['a3'], p_at_t['a4']])
    if 'a5' in p_at_t:
        out['fixedCollimation'] = core.G3VectorDouble([p_at_t['a5'], p_at_t['a6']])
    if 'az0' in p_at_t:
        out['az0'] = core.G3VectorDouble([p_at_t['az0']])

    return out

@core.indexmod
class ApplyPointingCorrection(object):
    '''
    Apply an additive correction (usually derived from HII-region
    very-fast-point obs) to the online pointing model.
    '''
    def __init__(self,                                           
                 in_model='OnlinePointingModel',
                 out_model='OfflinePointingModel',
                 correction='PointingModelCorrection',
                 source=['FieldScan','RCW38','MAT5A']):
        self.in_model = in_model
        self.out_model = out_model
        self.correction = correction
        self.source = source
        self.dpmodel = None
    
    def __call__(self, frame):
        if frame.type == core.G3FrameType.Calibration:
            correction = None
            if isinstance(self.source, list):
                for src in self.source:
                    if src + self.correction in frame:
                        correction = src + self.correction
                        break
            elif self.source + self.correction in frame:
                correction = self.source + self.correction
            if correction is None:
                return
            core.log_info("Found pointing model correction {}".format(correction))
            self.dpmodel = frame[correction]

        if frame.type == core.G3FrameType.Scan:
            model_orig = frame[self.in_model]
            model_corr = core.G3MapVectorDouble()
            if self.dpmodel is None:
                for key in model_orig.keys():
                    model_corr[key] = model_orig[key]
            else:
                for key in model_orig.keys():
                    try:
                        model_corr[key] = model_orig[key] + self.dpmodel[key]
                    except:
                        # core.log_warn("key " + key + " not found in PointingModelCorrection field")
                        model_corr[key] = model_orig[key]
            frame[self.out_model] = model_corr
#
#
#def OfflinePointingFromAdditiveCorrection(frame,
#                                          in_model='OnlinePointingModel',
#                                          out_model='OfflinePointingModel',
#                                          correction='PointingModelCorrection',
#                                          source='RCW38'):
#    if frame.type == core.G3FrameType.Scan:
#        model_orig = frame[in_model]
#        model_corr = core.G3MapVectorDouble()
#        if source+correction not in frame:
#            core.log_warn("PointingModelCorrection field not found in calibration frame; offline model will be identical to online mode.")
#            for key in model_orig.keys():
#                model_corr[key] = model_orig[key]
#        else:
#            dpmodel = frame[source+correction]
#            for key in model_orig.keys():
#                try:
#                    model_corr[key] = model_orig[key] + dpmodel[key]
#                except:
#                    core.log_warn("key " + key + " not found in both OnlinePointingModel and PointingModelCorrection fields")
#        frame[out_model] = model_corr
