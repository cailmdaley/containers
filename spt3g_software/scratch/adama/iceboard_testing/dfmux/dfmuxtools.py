# dfmuxtools.py
#
# Some tools for making maps and analyzing timestreams.
# Warning: some of these tools might be deprecated,
# only minimally useful, or poorly-written. This was 
# written in the course of learning the spt3g_software
# framework, so surely there are many non-optimalities.
#
# Adam Anderson
# adama@fnal.gov

from spt3g import core, gcp, dfmux, todfilter, calibration, mapmaker
from spt3g.mapmaker.mapmakerutils import MakeMap, ExtractTheMaps
import numpy as np
import os, calendar, time, re
import datetime as dt
import scipy.signal

import astropy.coordinates as astrocoord
from astropy.coordinates import SkyCoord
from astropy.coordinates import AltAz, ICRS
from astropy.time import Time
import astropy.units


# set up astropy data for the South Pole site to allow for unit
# conversions between Az/El and Ra/Dec
AltAz.location = astrocoord.EarthLocation.from_geodetic(lon=(-44. + 39./60.), lat=(-89 + 59./60. + 27.8/60./60.), height=2843.0)
AltAz.pressure = astropy.units.Quantity(0.69, astropy.units.bar)
AltAz.temperature = astropy.units.Quantity(-27.0, astropy.units.deg_C)
AltAz.humidity = 0.7
AltAz.obswl = astropy.units.Quantity(2.0, astropy.units.mm)

# variables needed for running 'generateIDF'. Annoyingly, python does not allow a function
# defined within another function to modify variables within the enclosing function. 
# Functions only can have scope over local and global variables. Thus, we make these
# variables global, even though they are only ever useful in the scope of 'generateIDF'.
arcfile_times = []
arcfile_ra = []
arcfile_dec = []
arcfile_scan_flag = []
cal_frame_added = False


def loadpointings(arcfiles):
    """
    Load the pointing information from specific ARC files.
    
    Parameters
    ----------
    arcfiles : list of full paths to the ARC files to load

    Returns
    -------
    times : list of times corresponding to each pointing
    ra : azimuth of each pointing
    dec : elevation of each pointing
    """
    times = []
    ra = []
    dec = []
    scan_flag = []

    def extractpointings(frame):
        if frame.type == core.G3FrameType.Timepoint:
            frame_az = frame['antenna0']['tracker']['actual'][0][0] * 180. / np.pi
            frame_el = frame['antenna0']['tracker']['actual'][1][0] * 180. / np.pi
            frame_datetime = dt.datetime.utcfromtimestamp(frame['antenna0']['tracker']['utc'][0][0].time / core.G3Units.second)
            frame_time = frame['antenna0']['tracker']['utc'][0][0].time / core.G3Units.second

            AltAz.obstime = Time(frame_datetime)
            coord = SkyCoord(frame=AltAz, az=astropy.units.Quantity(frame_az, astropy.units.degree), alt=astropy.units.Quantity(frame_el, astropy.units.degree))

            times.append(frame_time)
            ra.append(coord.icrs.ra.deg)
            dec.append(coord.icrs.dec.deg)
            scan_flag.append(frame['antenna0']['tracker']['scan_flag'][0][0])

            print frame['antenna0']['tracker']['utc'][0][0]
            #print "az = " + str(frame_az)
            #print "el = " + str(frame_el)
            #print "ra = " + str(coord.icrs.ra.deg)
            #print "dec = " + str(coord.icrs.dec.deg)

    for arcfile in arcfiles:
        print "Reading arcfile: " + arcfile
        pipe = core.G3Pipeline()
        pipe.Add(gcp.ARCFileReader(arcfile))
        pipe.Add(extractpointings)
        pipe.Run()

    times = np.array(times)
    ra = np.array(ra)
    dec = np.array(dec)
    scan_flag = np.array(scan_flag)
    
    return times, ra, dec, scan_flag


jframe = 0
def process_observation(starttime, stoptime, time_offset, raw_files, arcfile_path, outfile_path, scantype='GCP'):
    """
    Injects pointing information from ARC files into a 3g file containing raw
    time samples from the dfmux board. Also adds fake weights and calibration
    information which are needed in order for the mapmaker to run and collates
    raw "Timepoint" frames into "Scan" frames.
    
    Parameters
    ----------
    starttime : datetime object indicating the start time of the scan (UTC)
    stoptime : datetime object indicating the stop time of the scan (UTC)
    time_offset: offset in seconds to add to each data frame (see comment below)
    raw_files : list of raw data files to process
    arcfile_path : path to arcfiles
    outfile_path ; path for output file to create
    scantype : 'GCP' defines scans using the GCP scan flags, 'fixed' creates scans
        of a fixed length (hardcoded to be 10000)

    Returns
    -------
    [none]

    NB: The "TEST" timing in the 3G board is a piece of shit as it 
    produces timestamps with incorrect dates. The 'time_offset' is a number of
    seconds to add to all data frames so that the correct time is obtained.
    Please note that there is a function 'calc_time_offset' in this python
    module that will help you determine the correct number of seconds.
    Set this to '0' if you are sure that the timing is accurate. The IRIG-B
    timing should be accurate, but I have seen it a day off, so it is worth
    verifying.
    """
    global arcfile_times, arcfile_ra, arcfile_dec, arcfile_scan_flag, cal_frame_added
    cal_frame_added = False


    def check_time(frame):
        """
        Check whether time of this frame is in the range indended for
        processing.
        """
        if frame.type == core.G3FrameType.Timepoint:
            frame_timestamp = dt.datetime.utcfromtimestamp(frame['EventHeader'].time / core.G3Units.second)
            if frame_timestamp > starttime and frame_timestamp < stoptime:
                return frame
            else:
                return []
        else:
            print frame
        return frame

    def add_time_offset(frame, offset):
        """
        A piece of shit used for testing only!!! This adds a time offset to the
        frame times in order make the "TEST" time from the dfmux boards match
        some recent ARC files so that we can successfully extract pointing info.
        
        Parameters
        ----------
        offset : offset to add to each frame, in seconds
        """
        if frame.type == core.G3FrameType.Timepoint:
            t = frame['EventHeader'].time
            del frame['EventHeader']
            frame['EventHeader'] = core.G3Time(t + offset * core.G3Units.second)

    def inject_pointing(frame):
        """
        Inject pointing info into the frame by interpolating from the pointing info
        loaded from the current ARC file.
        """
        if frame.type == core.G3FrameType.Timepoint:
            frame_time = frame['EventHeader'].time / core.G3Units.second
            global arcfile_times, arcfile_ra, arcfile_dec, arcfile_scan_flag
            
            ra_interp = np.interp(frame_time, arcfile_times, arcfile_ra)
            dec_interp = np.interp(frame_time, arcfile_times, arcfile_dec)

            frame['BoresightRa'] = ra_interp * core.G3Units.deg
            frame['BoresightDec'] = dec_interp * core.G3Units.deg
            frame['scan_flag'] = arcfile_scan_flag[np.min(np.argmin(np.abs(arcfile_times - frame['EventHeader'].time / core.G3Units.second)))]

    def insert_cal(frame):
        """
        Insert a calibration frame at the beginning of each timestream. This is only
        useful insofar as it is needed in order to run the mapmaker.
        """
        global cal_frame_added
        if cal_frame_added == False:
            frame_to_add = core.G3Frame(core.G3FrameType.Calibration)
            frame_to_add['BolometerProperties'] = calibration.BolometerPropertiesMap()
            for boloname in frame['WiringMap'].keys():
                frame_to_add['BolometerProperties'][boloname] = calibration.BolometerProperties()
                frame_to_add['BolometerProperties'][boloname].x_offset = 0.0;
                frame_to_add['BolometerProperties'][boloname].y_offset = 0.0;
                frame_to_add['BolometerProperties'][boloname].band = 0.0;
                frame_to_add['BolometerProperties'][boloname].pol_angle = 0.0;
                frame_to_add['BolometerProperties'][boloname].pol_efficiency = 0.0;
            return_frame = [frame_to_add, frame]
            cal_frame_added = True
            return return_frame

    def print_frame(frame):
        """
        Print every 1000th frame.
        """
        global jframe
        jframe += 1
        if jframe % 1000 == 0 :
            print frame

    def add_weights(frame):
        """
        Add bolometer weights into each frame. This is only useful insofar as it is
        needed in order to run the mapmaker.
        """
        if frame.type == core.G3FrameType.Scan:
            frame['RawTimestreams_I_weights'] = core.G3MapDouble()
            frame['RawTimestreams_Q_weights'] = core.G3MapDouble()
            for boloname in frame['RawTimestreams_I'].keys():
                frame['RawTimestreams_I_weights'][boloname] = 1.0

            for boloname in frame['RawTimestreams_Q'].keys():
                frame['RawTimestreams_Q_weights'][boloname] = 1.0

    arcfilenames = get_arcfiles(starttime, stoptime, arcfile_path)
    arcfile_times, arcfile_ra, arcfile_dec, arcfile_scan_flag = loadpointings(arcfilenames)

    outfile = outfile_path + '/dfmux_idf_' + datetime2timestamp(starttime) + '_to_' + datetime2timestamp(stoptime) + '.g3'

    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=raw_files)
    pipe.Add(add_time_offset, offset=time_offset)
    pipe.Add(check_time)
    pipe.Add(insert_cal)
    pipe.Add(inject_pointing)
    pipe.Add(print_frame)
    if scantype == 'fixed':
        pipe.Add(dfmux.FixedLengthScans, N=10000)
    elif scantype == 'GCP':
        pipe.Add(dfmux.GCPFlaggedScans)
    else:
        print "ERROR: Invalid scan type selected!"
    pipe.Add(dfmux.DfMuxCollator)
    pipe.Add(add_weights)
    pipe.Add(core.G3Writer, filename=outfile)
    pipe.Run()


def make_map(infile):
    """
    Makes a map from a scan file.

    Parameters
    ----------
    infile : full path to the file containing the raw time samples from which
        to build map.

    Returns
    -------
    T_map : temperature map
    """
    # need to pass mapmaker inputs as arguments to this function!!
    out_map_parameters = coordinateutils.FlatSkyMap(x_len = 500, 
                                             y_len = 150,
                                             res = 0.1 * core.G3Units.deg, #remember to always use units when specifying things
                                             alpha_center = 15.0 * core.G3Units.deg,
                                             delta_center = -25.0* core.G3Units.deg,
                                             coord_ref = core.MapCoordReference.Equatorial,
                                             proj= mapmaker.MapProjection.Proj0) #map projections are stored in enums
    
    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=infile)
    pipe.Add(MakeMap, map_in=out_map_parameters,
             #poly_order = 4,
             map_id="test_out", ts_in_key="RawTimestreams_I",
             timestream_weight_key="RawTimestreams_I_weights",
             do_weight = False)
    #prints the contents of a frame
    def print_fr(frame):
        print frame
        if frame.type == core.G3FrameType.Scan:
            for k in frame['DetectorAlphaPointing'].keys():
                for p in frame['DetectorAlphaPointing'][k]:
                    print p
            for p in frame['BoresightRa']:
                print p
            for k in frame['DetectorDeltaPointing'].keys():
                for p in frame['DetectorDeltaPointing'][k]:
                    print p
            for p in frame['BoresightDec']:
                print p
    #pipe.Add(print_fr)

    #We use this object to get map data out
    map_extractor = ExtractTheMaps()
    pipe.Add(map_extractor)
    pipe.Run()

    T_map = map_extractor.maps["test_out"]["T"]
    Q_map = map_extractor.maps["test_out"]["Q"]
    U_map = map_extractor.maps["test_out"]["U"]

    return T_map, Q_map, U_map


class scan_counter(object):
    """
    Counts the number of scans in a file and their length.
    """
    def __init__(self, Output='nscans'):
        self.n_scans = 0
        self.scan_lengths = []
    def __call__(self, frame):
        if frame.type == core.G3FrameType.Scan:
            self.n_scans += 1
            self.scan_lengths.append(len(frame['RawTimestreams_I']))


class scan_extractor(object):
    """
    In a file that contains fixed-length scans, extracts
    the scans into an array.
    """
    def __init__(self, nscans, downsample_factor, Output="scan_extractor"):
        self.nscans = nscans
        self.downsample_factor = downsample_factor
        self.jscan = 0
        self.scans = dict()
    def __call__(self, frame):
        if frame.type == core.G3FrameType.Scan:
            for key in frame['RawTimestreams_I'].keys():
                n_samples = len(frame['RawTimestreams_I'][key]) / self.downsample_factor
                if key not in self.scans.keys():
                    self.scans[key] = np.zeros([self.nscans, n_samples])
                timestream = np.asarray(frame['RawTimestreams_I'][key])
                if self.scans[key].shape[1] == n_samples:
                    timestream = scipy.signal.resample(timestream, n_samples)
                    timestream = polyfilter(timestream, order=3)
                    self.scans[key][self.jscan, :] = timestream
            self.jscan += 1


class bolo_maps(object):
    """
    Extract maps from the scans of an observation.
    """
    def __init__(self, nscans, downsample_factor, Output="bolo_maps"):
        self.maps_left = dict()
        self.maps_right = dict()
        self.jscan_right = 0
        self.jscan_left = 0
        self.downsample_factor = downsample_factor
        self.mapwidth = None
        self.nscans = nscans
        self.mid_ra = None
        self.first_frame = True
    def __call__(self, frame):
        if frame.type == core.G3FrameType.Scan:
            if self.first_frame == True:
                self.first_frame = False
                self.mapwidth = round(1.1 * len(frame['BoresightRa']) / self.downsample_factor)

            # we align each scan to the median RA of the first scan,
            # which is hopefully close to the median RA of the field
            scan_ra = np.asarray(frame['BoresightRa'])
            if self.mid_ra == None:
                self.mid_ra = scan_ra[round(len(scan_ra)/2)]
            # determine whether scan is left-going or right-going
            if scan_ra[0] < scan_ra[-1]:
                scan_direction = +1
                maps = self.maps_right
                jscan = self.jscan_right
                print 'Left-going scan %d' % jscan
            else:
                scan_direction = -1
                maps = self.maps_left
                jscan = self.jscan_left
                print 'Right-going scan %d' % jscan
            scan_ra = scan_ra[::scan_direction]
            n_downsampled = round(len(scan_ra)/self.downsample_factor)
            scan_ra = scipy.signal.resample(scan_ra, n_downsampled)
            mid_ra_index = np.argmin(np.abs(self.mid_ra - scan_ra)) # need to unwrap RA in general

            if self.mapwidth == None:
                self.mapwidth = len(scan_ra)

            # align scans in center of map
            ncenter = self.mapwidth / 2
            min_map_ind = ncenter - mid_ra_index
            max_map_ind = ncenter - mid_ra_index + len(scan_ra)
            
            # copy the data out of the frame if the scan fits in the map array
            if min_map_ind >=0 and max_map_ind <= self.mapwidth:
                for key in frame['RawTimestreams_I'].keys():
                    if key not in maps.keys():
                        maps[key] = np.zeros([self.nscans/2 + 1, self.mapwidth]) + np.nan

                    timestream = np.asarray(frame['RawTimestreams_I'][key])
                    timestream = timestream[::scan_direction]
                    timestream = scipy.signal.resample(timestream, n_downsampled) 
                    timestream = polyfilter(timestream, order=3)

                    maps[key][jscan, min_map_ind:max_map_ind] = timestream
            else:
                print 'Map index error:', min_map_ind, max_map_ind, self.mapwidth, len(scan_ra)
            if scan_direction == 1:
                self.jscan_right += 1
            else:
                self.jscan_left += 1
        return frame


def make_bolo_maps(infile, downsample_factor):
    """
    Make single-bolometer maps. This is essentially just a 3G version of
    Jason Gallicchio's "skyarc" grid maps. Nothing fancy: we just assign 
    the timestreams to a matrix with 1 row per scan and rows aligned crudely
    using the mean RA of each scan.

    Parameters
    ----------
    infile : full path of file to read containing scans
    downsample_factor : downsampling factor

    Returns
    -------
    maps : python dictionary of maps, indexed by channel name
    """
    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=infile)
    count = scan_counter()
    pipe.Add(count)
    pipe.Run()
    n_scans = count.n_scans
    print 'File has %d scans'%n_scans

    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=infile)
    pipe.Add(core.Dump)
    bmaps = bolo_maps(n_scans, downsample_factor)
    pipe.Add(bmaps)
    pipe.Run()

    return bmaps.maps_left, bmaps.maps_right


def extract_fixed_scans(infile, downsample_factor):
    """
    Extract all the scans from an IDF with fixed-length scans as a dictionary
    of numpy arrays. This is useful for making "fake" maps of data taken when
    the telescope is docked.

    Parameters
    ----------
    infile : full path of file to read containing scans
    downsample_factor : factor by which to downsample scans

    Returns
    -------
    data : python dictionary of numpy arrays, containing 1 scan per row
    """
    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=infile)
    count = scan_counter()
    pipe.Add(count)
    pipe.Run()
    n_scans = count.n_scans

    # actually extract the scans
    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=infile)
    extractor = scan_extractor(n_scans, downsample_factor)
    pipe.Add(core.Dump)
    pipe.Add(extractor)
    pipe.Run()
    data = extractor.scans

    return data


def get_arcfiles(starttime, stoptime, arcfile_path):
    """
    Get a list of SPTpol arcfiles in a specific path which cover
    the interval between starttime and stoptime. This is dumb 
    and simply figures this out by the timestamps in the filename.

    Parameters
    ----------
    starttime : python datetime object representing the start time
    stoptime : python datetime object representing the stop time
    arcfile_path : path to search for files

    Returns
    -------
    arcfiles : a list of names of arcfiles
    """
    arc_filenames = np.sort(os.listdir(arcfile_path))
    timerange_strings = np.sort([fname.strip('.dat') for fname in arc_filenames])
    
    arcfiles = []
    for jfile in range(len(arc_filenames)):
        if starttime < timestamp2datetime(timerange_strings[jfile]) and \
           stoptime > timestamp2datetime(timerange_strings[jfile]):
            arcfiles.append(arcfile_path + '/' + arc_filenames[jfile])
        if jfile != len(arc_filenames)-1 and \
           starttime > timestamp2datetime(timerange_strings[jfile]) and \
           starttime < timestamp2datetime(timerange_strings[jfile+1]):
            arcfiles.append(arcfile_path + '/' + arc_filenames[jfile])

    return arcfiles


def timestamp2datetime(timestr):
    """
    Utility function for converting file timestamps of the form:
    'YYYYMMDD_hhmmss' to datetime objects.
    
    Parameters
    ----------
    timestr : a string containing the timestamp to convert

    Returns
    -------
    dtime : a datetime object corresponding to the time of the timestamp
    """
    dtime = dt.datetime(int(timestr[0:4]), int(timestr[4:6]), int(timestr[6:8]),
                        int(timestr[9:11]), int(timestr[11:13]), int(timestr[13:15]))
    return dtime


def datetime2timestamp(dtime):
    """
    Utility function for converting datetime objects into
    timestamp strings of the form 'YYYYMMDD_hhmmss'

    Parameters
    ----------
    dtime : a datetime object to convert

    Returns
    -------
    timestr : a string containing the timestamp to convert
    """
    timestr = '%d'%dtime.year + '%02d'%dtime.month + '%02d'%dtime.day + \
              '_%02d'%dtime.hour + '%02d'%dtime.minute + '%02d'%dtime.second
    return timestr


def polyfilter(timestream, order=1):
    """
    Quick and dirty polyfiltering of the timestream using some
    numpy functions.

    Parameters
    ----------
    timestream : the data timestream to filter
    order : the order of polynomial to fit

    Returns
    -------
    filtered_timestream : the filtered timestream
    """
    t = np.linspace(0., 1., num=len(timestream))
    coeffs = np.polyfit(t, timestream, order)
    p_fit = np.poly1d(coeffs)
    filtered_timestream = timestream - p_fit(t)
    return filtered_timestream


class timestream_extractor(object):
    def __init__(self, bolo_names, observation_length, Output='timestream_extractor'):
        self.observation_length = observation_length
        self.timestreams = dict()
        self.bolo_names = bolo_names
    def __call__(self, frame):
        if frame.type == core.G3FrameType.Scan:
            for bolo in self.bolo_names:
                if bolo not in self.timestreams.keys():
                    self.timestreams[bolo] = np.zeros(observation_length)
                nonzero_ind = np.nonzero(self.timestreams[bolo])[-1]
                if len(nonzero_ind) != 0:
                    start_ind = nonzero_ind[-1]
                else:
                    start_ind = 0
                self.timestreams[bolo][start_ind:(start_ind+len(self.timestreams[bolo]))] = self.timestreams[bolo]
        
def extract_timestreams(filename, bolos):
    """
    Returns the raw timestreams of some bolometers as one long numpy array.

    Parameters
    ---------
    filename : name of the file to read
    bolos : list of bolometers for which to get timestream info

    Returns
    -------
    timestreams : dict of {bolo: timestream}
    """
    # first do a loop over the data file to figure out the length
    # of the array to save
    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=filename)
    counter = scan_counter()
    pipe.Add(scan_counter)
    pipe.Run()
    obs_length = np.sum(np.asarray(scan_counter.scan_lengths))



class first_frame(object):
    """
    A patently stupid class that allows us to read in just the first
    timepoint frame in a G3 file.

    Attributes
    ----------
    first_frame : bool
        Indicates whether the first frame has been read.
    """
    def __init__(self, Output='first_frame'):
        self.first_frame = True
    def __call__(self, frame):
        if self.first_frame == False:
            return []
        if frame.type == core.G3FrameType.Timepoint and self.first_frame == True:
            self.first_frame == False

class store_frame_value(object):
    """
    Let's try to be generic here and write a class that will store
    arbitrary values from frames in a g3 file for manipulation.

    Attributes
    ----------
    stored_data : dictionary of numpy arrays, indexed by variable name
        The data to store for retrieval later.
    access_functions : list of functions to access data from the frames (see below)
    jentry : the number of entry that we are on in reading the data
    frame_type : type of frame from which to retrieve data
    """
    def __init__(self, data_size, frame_type, access_functions, Output='store_frame_value'):
        self.jentry = 0
        self.data = dict()
        self.frame_type = frame_type
        self.access_functions = access_functions
        for func in self.access_functions:
            self.data[func.func_name] = np.zeros(data_size)
    def __call__(self, frame):
        if frame.type == self.frame_type and self.jentry < len(self.data[self.data.keys()[0]]):
            for func in self.access_functions:
                self.data[func.func_name][self.jentry] = func(frame)


def calc_time_offset(filename, starttime):
    """
    Some fecal matter to account for the fact that the timestamp on each g3 frame
    is incorrect when data is taken in TEST mode. There is some large global
    offset that needs to be added to the data--possibly a different one each time
    the iceboard is power cycled. The time offset appears to be stable between 
    power cycles, although this is not necessarily guaranteed. This function will
    calculate the timing offset, given the true start time of the data file. At
    minimum, you should be able to obtain the true start time of the file from 
    the filename. Like fecal matter, this appears to be disgusting but unavoidable.
    
    Parameters
    ----------
    filename : full path to the file with which to compute offsets
    starttime : a Python datetime object representing the true starting time of
        the first frame in the file

    Returns
    -------
    offset : the offset is seconds
    """
    def frametime(frame):
        return frame['EventHeader'].time / core.G3Units.second
    data_storer = store_frame_value(1, core.G3FrameType.Timepoint, [frametime])

    first_frame_checker = first_frame()

    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=filename)
    pipe.Add(first_frame_checker)
    pipe.Add(data_storer)
    pipe.Run()

    epoch_starttime = time.mktime(starttime.timetuple()) - time.altzone
    frame_time = data_storer.data['frametime'][0]
    offset = epoch_starttime - frame_time

    return offset
