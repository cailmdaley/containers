from spt3g import core
from spt3g.sources import source_utils
import numpy as np
import os, glob
import astropy
import time

def ang_diam_from_gcp_ephem(source,mjd,gcp_dir=None):
    '''
    get distance from GCP ephem file & calculate angular diameter from
    that.  agrees with JPL Horizons to <0.01 arcsec.
    '''

    units = core.G3Units
    radii = {"Mercury":2440.*units.km,
             "Venus":6052.*units.km,
             "Mars":3397.*units.km,
             "Jupiter":71492.*units.km,
             "Saturn":60268.*units.km,
             "Uranus":25559.*units.km,
             "Neptune":24766.*units.km}

    goodsources = radii.keys()
    isgood = source.capitalize() in goodsources
    if not(isgood):
        core.log_warn('Source '+source+' not in list of planet radii.')
        return -1

    if gcp_dir is None:
        gcp_dir = os.getenv('GCP_DIR')
        if gcp_dir is None:
            core.log_warn("Can't find environment variable $GCP_DIR, please specify GCP location.")
            return -2
#    ephemfile = gcp_dir+'/control/ephem/'+source.lower()+'.ephem'
    ephemfile = gcp_dir+'/config/ephem/'+source.lower()+'.ephem'
#    print(ephemfile)
    data = ascii.read(ephemfile)
    mjds = np.zeros(len(data))
    distances = np.zeros(len(data))
    for i in np.arange(len(data)):
        dtemp = data[i]
        mjds[i] = dtemp[0]
        distances[i] = dtemp[3]
    thisdist = np.interp(mjd,mjds,distances)*units.au
    thisangdiam_radians = 2.*radii[source.capitalize()]/thisdist

    return thisangdiam_radians

def temperature_solid_angle_product(source,mjd,band,gcp_dir=None,TCMB=True):
    '''
    get integrated temperature-solid angle product (in
    Kelvins-steradians) from angular diameter (calculated above) and
    thermodynamic temperature. temps come from Planck (1612.07151) for
    Mars and beyond and the abstract of doi:10.1016/0019-1035(78)90082-9
    for Mercury and Venus. latter source is only measured at 1mm (270
    GHz), so frequency-independent temperature assumed. all values
    assume time independence. unset TCMB keyword to get answer in
    Rayleigh-Jeans temp. if someone wants to dig up the JCMT FLUXES
    code, that would be an improvement on this.
    '''

    angdiam_radians = ang_diam_from_gcp_ephem(source,mjd,gcp_dir=gcp_dir)
    if angdiam_radians < 0.:
        return -1

    band_ghz = band/core.G3Units.Hz/1e9
    planck_bands_ghz = np.asarray([100.,143.,217.,353.,545.,857.])
    npb = len(planck_bands_ghz)
    temps = {"Mercury":np.zeros(npb)+320.,
             "Venus":np.zeros(npb)+276.,
             "Mars":np.asarray([194.3,198.4,201.9,209.9,209.2,213.5]),
             "Jupiter":np.asarray([172.3,173.6,174.7,166.3,136.5,160.3]),
             "Saturn":np.asarray([145.7,147.0,144.9,141.5,102.4,115.5]),
             "Uranus":np.asarray([120.5,108.4,98.5,86.2,73.9,66.2]),
             "Neptune":np.asarray([117.4,106.4,97.4,82.6,72.3,65.3])}
    thesetemps = temps[source.capitalize()]
    thistemp = np.interp(band_ghz,planck_bands_ghz,thesetemps)
    if TCMB:
        thistemp *= source_utils.trj_to_tcmb(band_ghz)

#    print("Using "+str(thistemp)+"K as thermodynamic temperature for source "+source.capitalize()+" at frequency "+str(band)+" GHz.")

    result = thistemp*np.pi*(angdiam_radians/2.)**2

    return result

def sun_distance(source,mjd,gcp_dir=None,input_ra=None,input_dec=None):
    '''
    get distance from sun (in degrees) using info from 
    GCP ephem file. agrees with JPL Horizons to < 5 arcsec.
    '''

    units = core.G3Units

    if gcp_dir is None:
        gcp_dir = os.getenv('GCP_DIR')
        if gcp_dir is None:
            core.log_warn("Can't find environment variable $GCP_DIR, please specify GCP location.")
            return -2

    if input_ra is None or input_dec is None:
        ephemfile = gcp_dir+'/config/ephem/'+source.lower()+'.ephem'
        data = ascii.read(ephemfile)
        srcmjds = np.zeros(len(data))
        srcra = np.zeros(len(data))
        srcdec = np.zeros(len(data))
        for i in np.arange(len(data)):
            dtemp = data[i]
            srcmjds[i] = dtemp[0]
            srcra[i] = astropy.coordinates.Angle(dtemp[1],unit='hour').deg
            srcdec[i] = astropy.coordinates.Angle(dtemp[2],unit='deg').deg
        this_srcra = np.interp(mjd,srcmjds,srcra)
        this_srcdec = np.interp(mjd,srcmjds,srcdec)
    else:
        this_srcra = input_ra
        this_srcdec = input_dec

    sunephemfile = gcp_dir+'/config/ephem/sun.ephem'
    sundata = ascii.read(sunephemfile)
    sunmjds = np.zeros(len(sundata))
    sunra = np.zeros(len(sundata))
    sundec = np.zeros(len(sundata))
    for i in np.arange(len(sundata)):
        dtemp = sundata[i]
        sunmjds[i] = dtemp[0]
        sunra[i] = astropy.coordinates.Angle(dtemp[1],unit='hour').deg
        sundec[i] = astropy.coordinates.Angle(dtemp[2],unit='deg').deg
    this_sunra = np.interp(mjd,sunmjds,sunra)
    this_sundec = np.interp(mjd,sunmjds,sundec)
    k1 = astropy.coordinates.SkyCoord(this_srcra, this_srcdec, frame='icrs', unit='deg')
    k2 = astropy.coordinates.SkyCoord(this_sunra, this_sundec, frame='icrs', unit='deg')
    thisdist = astropy.coordinates.SkyCoord.separation(k1,k2).deg

    return thisdist

def show(source,gcp_dir=None):
    '''
    meant to mimic "show" command in GCP viewer. returns
    RA & dec of source at current time.
    '''

    units = core.G3Units

    if gcp_dir is None:
        gcp_dir = os.getenv('GCP_DIR')
        if gcp_dir is None:
            core.log_warn("Can't find environment variable $GCP_DIR, please specify GCP location (with kwarg gcp_dir).")
            return -2

    utcnow = time.gmtime()
    timenow = core.G3Time(str(utcnow.tm_year)+str('%02i' % utcnow.tm_mon)+str('%02i' % utcnow.tm_mday)+'_'+str('%02i' % utcnow.tm_hour)+str('%02i' % utcnow.tm_min)+str('%02i' % utcnow.tm_sec))
    mjd = timenow.mjd

    sun_and_planets = ['sun','mercury','venus','mars','jupiter','saturn','uranus','neptune','pluto']
    isplanet = False
    if source.lower() in sun_and_planets:
        isplanet = True
        ephemfile = gcp_dir+'/config/ephem/'+source.lower()+'.ephem'
        data = ascii.read(ephemfile)
        srcmjds = np.zeros(len(data))
        srcra = np.zeros(len(data))
        srcdec = np.zeros(len(data))
        srcdist = np.zeros(len(data))
        for i in np.arange(len(data)):
            dtemp = data[i]
            srcmjds[i] = dtemp[0]
            srcra[i] = astropy.coordinates.Angle(dtemp[1],unit='hour').deg
            srcdec[i] = astropy.coordinates.Angle(dtemp[2],unit='deg').deg
            srcdist[i] = dtemp[3]
        this_srcra = np.interp(mjd,srcmjds,srcra)
        this_srcdec = np.interp(mjd,srcmjds,srcdec)
        this_srcdist = np.interp(mjd,srcmjds,srcdist)
    else:
        ephemfile = gcp_dir+'/config/ephem/source.cat'
        f = open(ephemfile)
        nextline = '#foo'
        thissrc = 'nosourceyouhaveeverheardof'
        while thissrc.lower() != source.lower() and nextline != '':
            newline = nextline
            nextline = f.readline()
            if newline[0] == '#':
                continue
            strs = (' '.join(newline.split())).split(' ')
            if strs[0].upper() == 'J2000':
                thissrc = strs[1]
                rastr = strs[2]
                decstr = strs[3]
        this_srcra = astropy.coordinates.Angle(dtemp[1],unit='hour').deg
        this_srcdec = astropy.coordinates.Angle(dtemp[2],unit='deg').deg
    
    k = astropy.coordinates.SkyCoord(this_srcra, this_srcdec, frame='icrs', unit='deg')
#    if isplanet:
#        k = astropy.coordinates.SkyCoord(this_srcra*astropy.units.deg, this_srcdec*astropy.units.deg, distance=this_srcdist*astropy.units.au, frame='icrs')
#    else:
#        k = astropy.coordinates.SkyCoord(this_srcra, this_srcdec, frame='icrs', unit='deg')
    loc = astropy.coordinates.EarthLocation(lon=314.211*astropy.units.deg, lat=-89.993*astropy.units.deg, height=2997.29*astropy.units.m)
    aptime = astropy.time.Time(mjd,format='mjd')
    kt = k.transform_to(astropy.coordinates.AltAz(obstime=aptime,location=loc))

    print("Source: " + source.upper() + "  (" + timenow.Summary()[0:20] + ")")
    print(" AZ: " + str(kt.az.to_string(unit='hour',sep=':')) + "  EL: " + str(kt.alt.to_string(sep=':')))
    print(" RA: " + dtemp[1] + " DEC: " + dtemp[2])
    if kt.alt.deg > 0.:
        sundist = sun_distance('',mjd,input_ra=this_srcra,input_dec=this_srcdec,gcp_dir=gcp_dir)
        if source.lower() == 'sun':
            print(" Currently above horizon.\n")
        else:
            print(" Currently above horizon, currently " + "%6.2f" % sundist + " degrees from sun.\n")
    else:
        print(" Currently below horizon.")

    return
