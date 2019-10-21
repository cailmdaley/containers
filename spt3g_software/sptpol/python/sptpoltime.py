"""
Perform date/time functions that we want in the SPT code. The big attraction here
is the SptDatetime object, which does string parsing and addition/subtraction.

Functions:
   date_toolkit     : For backward compatibility. Mirrors IDL's "date_toolkit" function.
   datetime_to_jul  : Convert a datetime object into a Julian day.
   datetime_to_mjd  : Convert a datetime object into a modified Julian day.
   now              : Current (local) time.
   utcnow           : Current (UTC) time.
   midpoint         : The time halfway between the two input times.
   toSptDatetime    : Convert arrays of strings or floats to SptDatetime objects.
   mjdToSptDatetime : As toSptDatetime, but assumes the input is a MJD float.
   doubleStringToSptDatetime : Takes a string with two dates and returns a tuple of SptDatetimes.
   filenameToTime   : Extracts a date tag from a filename or list of filenames.

Classes:
   SptDatetime      : Parse date/time strings in different formats, output date/time strings
                in different formats, and find the sum and difference of dates.

Non-Core Dependencies:
   ephem
   NumPy (not used frequently)
"""

__metaclass__ = type  #Use "new-style" classes
__author__    = ["Stephen Hoover"]
__email__     = ["hoover@kicp.uchicago.edu"]
__version__   = "1.05"
__date__      = "2013-10-29" #Date of last modification

from datetime import datetime
from datetime import timedelta
import re
import numpy as np
from spt3g.util.genericutils import class_name 

sec_per_day  = 86400.

# It speeds us up a bit to have these pre-made.
filename_date_format = re.compile(r'\d{8}_\d{6}')
one_us = timedelta(microseconds=1)

#################################
def date_toolkit(input, output_format='arc', add=None, units=None):
    """
    This method replicates the behavior of the IDL "date_toolkit" function
    by using utils.SptDatetime objects. Valid "units" are those accepeted
    by the datetime.timedelta class: 'microseconds', 'milliseconds', 'seconds',
    'minutes', 'hours', 'days', or 'weeks'.
    """

    # Check if the input is string-like.
    try:
        stringlike = bool(input+"")
    except TypeError:
        stringlike = False

    # First check if we got an iterable of dates as input. If so, process all of them
    # and return the result as an array. If the input is a list with 6 members,
    # the last 5 of which are between 0 and 60, then assume it's a single "array" input.
    try:
        if not (stringlike or len(input)==6 and 
                [ isinstance(single_input,(int,float,long)) and (0<=single_input<=60)
                                        for single_input in input[1:]]):
                                            return [ date_toolkit(single_input, output_format, add, units)
                                                    for single_input in input]
    except TypeError:
        pass

    # If the input was a single date, use the SptDatetime object to
    # handle any time addition/subtraction and format conversion.
    date_obj = SptDatetime(input)
    if add is not None and units is not None:
        date_obj += timedelta(**{units:add})
    return date_obj.formatDate(output_format)

#################################
def toSptDatetime(input, asarray=False):
    """
    Converts the input to SptDatetime object(s). This function can handle
    both single time inputs and arrays of time inputs.
    
    INPUTS
        input: Either a single valid input to the SptDatetime constructor,
            or an iterable of such inputs.
        
        asarray [False]: (bool) If True, the output will be a NumPy array.
            Automatically set if the input is an array.
            
    OUTPUT
        A single SptDatetime object (if the input is a single time), 
        a NumPy array of SptDatetime objects (if the input is an array, or
        if asarray==True), or a list of SptDatetime objects.
    """
    try:
        return SptDatetime(input)
    except TypeError:
        time_list = [SptDatetime(_time) for _time in input]
        if asarray or isinstance(input, np.ndarray): time_list = np.array(time_list, dtype=object)
        return time_list
    
#################################
def mjdToSptDatetime(_input, asarray=False):
    """
    Assumes the input is a float which represents a date in MJD format.
    Converts the input to SptDatetime object(s). This function can handle
    both single time inputs and arrays of time inputs.
    
    INPUTS
        _input: Either a single valid MJD, or an iterable of MJD dates.
        
        asarray [False]: (bool) If True, the output will be a NumPy array.
            Automatically set if the input is an array.
            
    OUTPUT
        A single SptDatetime object (if the input is a single time), 
        a NumPy array of SptDatetime objects (if the input is an array, or
        if asarray==True), or a list of SptDatetime objects.
    """
    import ephem

    try:
        return SptDatetime(_input)
    except TypeError:
        # The ephem module wants inputs in DJD, so start by converting from MJD to DJD.
        _input = _input - 15019.5
        
        # Now convert everything to datetime objects.
        _input = map(lambda x: ephem.Date(x).datetime(), _input)
        
        # Finally, convert them to SptDatetime objects.
        time_list = [SptDatetime(_time) for _time in _input]
        
        if asarray or isinstance(_input, np.ndarray): time_list = np.array(time_list, dtype=object)
        return time_list

#################################
def doubleStringToSptDatetime(input_string, *extra_string):
    """
    Takes a string which specifies two SptDatetime objects, splits it up, and converts to 
    a pair of SptDatetime objects.
    
    INPUTS
        input_string: (string) We will first attempt to divide the string by splitting it
            at any whitespace. If this fails to produce a first and last string which are
            each parsable by SptDatareader, then we will divide the string exactly in half.
            
        *extra_string: Optionally, provide the two date strings as two separate arguments.
            
    OUTPUT
        A 2-tuple of SptDatetime objects.
    """
    if len(extra_string)>0:
        input_string = input_string+extra_string[0]
    string_pair = input_string.split()
    try:
        return SptDatetime(string_pair[0]), SptDatetime(string_pair[-1])
    except Exception:
        return SptDatetime(input_string[:len(input_string)/2].strip()), SptDatetime(input_string[len(input_string)/2:].strip())
    

# Old versoins of ephem.Date ignore microseconds when converting a datetime object to a float.
# (Going from a float to a datetime always included microseconds, so this is only relevant in the two functions below.)
# Rather than test for this each time, do it once when this module is loaded and store the value.
# Internally, ephem.Date stores a datetime object as a float Dublin Julian Day, whose Epoch is 12h Dec 31, 1899
# for the new pyephem, this would be a half-second:  0.5/(24*60*60)  and the condition would evaluate to False

#this code is now below to help with not needing pyephem to run things that don't matter

#################################
def datetime_to_jul(datetime_obj):
    """
    Returns the Julian Date corresponding to an input datetime object.
    This function uses the ephem.Date class to do the heavy lifting.
    
    If you want to avoid excessive roundoff, use datetime_to_mjd() instead.
    JD Epoch is 12h Jan 1, 4713 BC
    """
    import ephem
    THIS_IS_OLD_PYEPHEM_THAT_IGNORES_MICROSECONDS = (ephem.Date(datetime(1899, 12, 31, 12, 0, 0, 500000)) == 0.0)

    # Note that ephem.Date stores a datetime object as a float Dublin Julian Day.
    if THIS_IS_OLD_PYEPHEM_THAT_IGNORES_MICROSECONDS:
        # Old versoins of ephem.Date truncate microseconds, so add that back in by hand.
        return float(ephem.Date(datetime_obj)) + 2415020.0 + (datetime_obj.microsecond*1e-6/sec_per_day)
    else:
        return float(ephem.Date(datetime_obj)) + 2415020.0
        
#################################
def datetime_to_mjd(datetime_obj):
    """
    Returns the modified Julian Date corresponding to an input datetime object.
    This function uses the ephem.Date class to do the heavy lifting.
    
    It's useful to have this around as well as datetime_to_jul, since the Julian
    date is such a large number - digitization error reduces accuracy on 
    date that have been converted to Julian dates.
    MDJ Epoch is 0h Nov 17, 1858
    """

    import ephem
    THIS_IS_OLD_PYEPHEM_THAT_IGNORES_MICROSECONDS = (ephem.Date(datetime(1899, 12, 31, 12, 0, 0, 500000)) == 0.0)
    # Note that ephem.Date stores a datetime object as a float Dublin Julian Day.
    if THIS_IS_OLD_PYEPHEM_THAT_IGNORES_MICROSECONDS:
        # Old versoins of ephem.Date truncate microseconds, so add that back in by hand.
        return float(ephem.Date(datetime_obj)) + 15019.5 + (datetime_obj.microsecond*1e-6/sec_per_day)
    else:
        return float(ephem.Date(datetime_obj)) + 15019.5

#################################
def now(tz=None):
    """
    Current time in timezone tz, found by calling datetime.now(tz=tz).
    
    OUTPUT
       SptDatetime object with the current time.
    """
    return SptDatetime(datetime.now(tz=tz))
    
#################################
def utcnow():
    """
    Current UTC time, found by calling datetime.utcnow().

    OUTPUT
       SptDatetime object with the current UTC time.
    """
    return SptDatetime(datetime.utcnow())

#################################
def midpoint(start_time, stop_time):
    """
    Returns the time midway between the two times entered.
    The times are named "start time" and "stop time", but
    they don't need to be in chronological order.
    
    INPUTS
        start_time: (SptDatetime) The first time.
        
        stop_time: (SptDatetime) The second time.
        
    OUTPUT
        An SptDatetime which is midway between the two input times.
    """
    try:
        return start_time + (stop_time - start_time)/2
    except TypeError:
        start_time, stop_time = SptDatetime(start_time), SptDatetime(stop_time)
        return start_time + (stop_time - start_time)/2

#################################
def filenameToTime(filename, identical_formats=False, assume_filename_format=True, 
                   fill_not_a_dates=False, quiet=False):
    """
    Converts input file name(s) which contain "file" format timestamps to
    SptDatetime object(s).
    
    INPUTS
        filename: (string or iterable of strings) Either a filename, or an iterable
            which contains many filenames.
            
        identical_formats [False]: (bool) If "filename" is a list, setting this to
            True lets us assume that the date/time string starts at the same
            character in each string. This saves some computation time on 
            especially long lists. (How much? Seems like it could be a lot.)
            Note that the start and end points of the date string are referenced
            to the end of the filename, so, as long as the string is in the same
            place relative to the end (eg *YYYYMMDD_HHMMSS.fits), you may set
            identical_formats=True.
    
        assume_filename_format [True]: (bool) If we know that the input timestamps are
            in the filename format (YYYYMMDD_HHMMSS), then we can speed up the conversion
            to an SptDatetime object by quite a bit. We assume that anyway, with our
            use of filename_date_format, so there's no reason not to use this option.
            
        fill_not_a_dates [False]: (bool) If True, and the input is a list, put None
            in the output list for any file which can't be turned into a date.
            If this argument is True, the output list will have one entry for each
            member of the input list. If this argument is False, that might not be true.
            
        quiet [False]: (bool) If True, suppress screen output. 
            
    OUTPUT
        If the input was a single string, then the corresponding SptDatetime object
        is output. Otherwise, a list of corresponding SptDatetime objects is output.
    """
    
    # First, try to treat the input as a single string and search for the
    # file format regex inside. If that fails, treat it as a list.
    try:
        output_time = SptDatetime( filename_date_format.findall(filename)[-1] )
    except TypeError:
        output_time = []
        if identical_formats and len(filename)!=0:
            # If we're assuming that the files have the date string in the same place, analyze
            # the first filename to see what that place is. Index from the end of the string - 
            # this is much more robust to changes in filename, since the date string is usually
            # at the end of the filename.
            first_datestring = filename_date_format.findall(filename[0])[-1]
            datestring_start = filename[0].rfind(first_datestring) - len(filename[0])
            datestring_end = datestring_start + len(first_datestring)
            # For quicker processing, assume that the time is in the filename date format.
            if assume_filename_format and datestring_end-datestring_start == 15:
                for name in filename:
                    datestring = name[datestring_start:datestring_end]
                    try:
                        output_time += [SptDatetime(datetime(int(datestring[:4]), int(datestring[4:6]), int(datestring[6:8]),
                                                             int(datestring[9:11]), int(datestring[11:13]), int(datestring[13:15])))]
                    except (IndexError,ValueError):
                        if not quiet: print("Filename %s doesn't contain a date." % name)
                        if fill_not_a_dates: output_time += [None]
                        continue
                
            else:
                for name in filename:
                    try:
                        output_time += [SptDatetime(name[datestring_start:datestring_end])]
                    except TypeError:
                        if not quiet: print("Filename %s doesn't contain a date." % name)
                        if fill_not_a_dates: output_time += [None]
                        continue
                    
        else:
            for name in filename:
                try:
                    output_time += [SptDatetime( filename_date_format.findall(name)[-1] )]
                except IndexError:
                    if not quiet: print("Filename %s doesn't contain a date." % name)
                    if fill_not_a_dates: output_time += [None]
                    continue

        # If the input was an array, then the output should be too.
        if isinstance(filename, np.ndarray):
            output_time = np.asarray(output_time, dtype=object)
    
    return output_time
        
#################################
#################################

# Useful conversion factors
_jul_conversions = {'djd':-2415020.0, 'jul':0.0, 'mjd':-2400000.5, 'uni':-2440587.5}

class SptDatetime(datetime):
    """
    An extension of datetime useful for dates which must be input and output in
    formats common to SPT analysis code. Use the timedelta object to add and
    subtract times from this object.

    This object stores times to the nearest 10 microseconds (10^-5 of a second).

    Add additional formats by including their format strings in the
    "_date_formats" member. If the format is related to the Julian Date by
    an additive constant, put the format name and additive constant in
    "_jul_conversions".
    """

    # Define the date formats for input and output. Note that the Julian Date
    # and associated formats (MJD, DJD) are not covered by the datetime format
    # strings, so I'll handle those separately.
    # Use the first three letters (lowercase) to define each format.
    _date_formats = {'arc':'%d-%b-%Y:%H:%M:%S', 'fil':'%Y%m%d_%H%M%S',
                     'fin':'%Y%m%d_%H%M%S.log', 'fis':'%Y%m%d', 'fih':'%Y-%m-%d',
                     'log':'%y%m%d %H:%M:%S', 'arr':'[%Y, %m, %d, %H, %M, %S]',
                     'dmy':'%d %B %Y', 'mdy':'%m/%d/%y', 'mdf':'%m/%d/%Y',
                     'fan':'%A, %d %B %Y', 'ful':'%A, %d %B %Y, %H:%M:%S',
                     'eph':'%Y/%m/%d %H:%M:%S', 'ard':'%d-%b-%Y',
                     'arm':'%d-%b-%Y:%H:%M:%S.%f',
                     'sky':'%m/%d/%y %H:%M', 'dfm':'%Y-%m-%d-%H:%M:%S'}
    
    # Constants for converting from JD to related date formats.
    # The ephem.Date class uses DJD, so establish a separate dictionary for that.
    # Finally, make a MJD dictionary as well. 
    # 'uni' is the "Unix" epoch, i.e., time since midnight on 1 January 1970.
    _djd_conversions = dict([ (key, num-_jul_conversions['djd']) for key, num in _jul_conversions.items()])
    _mjd_conversions = dict([ (key, num-_jul_conversions['mjd']) for key, num in _jul_conversions.items()])

    # A list of all valid output formats. Not actually used in the code.
    valid_output_formats = list(_date_formats.keys()) + list(_jul_conversions.keys())
    
    # "datetime" is an immutable type, so it's initialized with
    # a "__new__" function rather than an "__init__" function.
    # Therefore SptDatetime must overload the __new__ 
    # function as its constructor.
    def __new__(cls, _input):
        """ Initializes the object from a single input. The input may be any one of
            the following:
              'DD-MMM-YYYY:HH:MM:SS'                         : Archive format
              'DD-MMM-YYYY:HH:MM:SS.UUUUUU'                  : Archive format with microseconds
              'DD-MMM-YYYY'                                  : Archive format without time of day
              'YYYYMMDD_HHMMSS' or 'YYYYMMDD_HHMMSS.log' or 'YYYYMMDD' or 'YYYY-MM-DD' : Filename format
              [Y,M,D,H,M,S]                                  : Date array 
              'YYMMDD HH:MM:SS'                              : Logentry timestamp
              'YYYY/MM/DD HH:MM:SS'                          : PyEphem date string format
              X.X                                            : Julian date (if greater than 6500*365)
              X.X                                            : Modified Julian day (otherwise)
              'now'                                          : The current time, in UTC.
              datetime.datetime                              : An existing datetime.datetime object.
              'DD/MM/YY' or 'DD/MM/YYYY'
              'Nameday, DD Monthname YYYY' or 'DD Monthname YYYY'
              'MM/DD/YY HH:MM'                               : Format used in Skynet logs
            Invalid input formats will generate a TypeError exception.
            """

        # Check input and create the datetime object. Cases: 'now', existing datetime object,
        # a number (treated as JD or MJD, depending on value), or something in the
        # _date_formats dictionary.
        import ephem
        new_date = None
        if isinstance(_input, (float,int,long,np.integer,np.floating)):
            # If we get a single number as input, then we'll assume it's a
            # JD or MJD formatted date. I'll use the "ephem.Date" class
            # to convert them. Note that ephem.Date takes inputs
            # in DJD format, so convert to that first!
            input_type = 'mjd' # Start by assuming MJD. Next check if it could be JD.
            if _input > 6500.*365.: # This is true if we get a JD (epoch 4713 BC!)
                input_type = 'jul'
            new_date = ephem.Date(_input-SptDatetime._djd_conversions[input_type]).datetime()
        elif isinstance(_input, datetime): # Note that this is also true for other SptDatetime objects!
            new_date = _input
        elif _input is None or _input=="None":
            return None
        elif str(_input).lower()=='now':
            new_date = datetime.utcnow()
        else:
            for _format in SptDatetime._date_formats.itervalues():
                try:
                    new_date = datetime.strptime(str(_input), _format)
                    break
                except (ValueError, TypeError):
                    pass

        # Make sure that we got an input that we can understand.
        if new_date is None:
            # Give one last hurrah - try giving the input to the datetime constructor.
            # We'll use this if we're unpickling.
            try:
                new_date=datetime(_input)
            except Exception:
                raise TypeError("Sorry, I can't parse "+'"'+str(_input)+'".')

        # Create the new object. Round off the last digit of the microseconds; we
        # aren't that accurate, and the occasional '9999' is annoying me.
        if new_date.microsecond>=999995:
            new_date += (int(1e6)-new_date.microsecond)*one_us
        new_date = super(SptDatetime, cls).__new__(cls, new_date.year, new_date.month,
                                                   new_date.day, new_date.hour,
                                                   new_date.minute, new_date.second,
                                                   (new_date.microsecond+5)/10*10,
                                                   new_date.tzinfo)
        return new_date

    
    #################################
    def formatDate(self, outformat=None):
        """ Outputs the date formatted in the desired manner. Valid output formats are:
              arc[hive]: 'DD-MMM-YYYY:HH:MM:SS'
              arm (Archive format with microseconds) : 'DD-MMM-YYYY:HH:MM:SS.UUUUUU'
              fil[e]: 'YYYYMMDD_HHMMSS'
              fin (or "file_name"): 'YYYYMMDD_HHMMSS.log'
              fis (or "file_short"): 'YYYYMMDD'
              fih (or "file_hyphen"): 'YYYY-MM-DD'
              arr[ay]: [Y,M,D,H,M,S]
              log[entry]: 'YYMMDD HH:MM:SS'
              eph[em]: 'YYYY/MM/DD HH:MM:SS'
              jul[ian]: Julian date
              mjd: Modified Julian day
              djd: Dublin Julian day
              dmy: Day Month Year - DD Monthname YYYY
              mdy: Month/Day/Year - MM/DD/YY
              mdf: Month/Day/Full year - MM/DD/YYYY
              fan[cy]: Nameday, DD Monthname YYYY
              ful[l]: Full time and date output: Nameday, DD Monthname YYYY, HH:MM:SS
              uni[x]: Unix time, seconds since midnight 1 Jan 1970.
              sky[net]: 'MM/DD/YY HH:MM', the format used in Skynet logs.
        The format selection is case-independent and only relies on the
        first three characters. Unrecognized output formats will
        raise a ValueError exception.
        
        Any fractional seconds are ignored in output modes that don't output
        anything smaller than seconds.
        """
        
        # Set the default output format. Archive format.
        # Display microseconds only if present.
        if outformat is None:
            if self.microsecond==0: outformat='arc'
            else: outformat='arm'
        
        # Treat the output format as a string, force it lowercase,
        # and only look at the first three characters.
        if outformat=='file_name':
            reduced_outformat='fin' # Translate to one of the more obscure three-letter abbreviations.
        elif outformat.startswith('file_short'):
            reduced_outformat='fis' # Translate to one of the more obscure three-letter abbreviations.
        elif outformat.startswith('file_hyphen'):
            reduced_outformat='fih' # Translate to one of the more obscure three-letter abbreviations.
        else:
            reduced_outformat = str(outformat).lower()[0:3]

        # If the format is something with a recognized format string,
        # use that to format the output. Otherwise, check if it could
        # be a Julian Date or related format.
        if reduced_outformat in self._date_formats:
            try:
                date_out = self.strftime(self._date_formats[reduced_outformat])
            except ValueError:
                date_out = "--invalid SptDatetime value--"
        elif reduced_outformat=='uni':
            date_out = datetime_to_mjd(self) + self._mjd_conversions[reduced_outformat]
            date_out *= sec_per_day
            date_out = round(date_out, 5) # Because of digitization/rounding errors, we're only good to 10^-5 s. 
        elif reduced_outformat in _jul_conversions:
            # Base the conversion on the modified Julian day instead of
            # the Julian day to reduce digitization errors.
            date_out = datetime_to_mjd(self) + self._mjd_conversions[reduced_outformat]
        else:
            raise ValueError(str(outformat)+" is an unrecognized output format!")

        # Array output requires some further processing. The array came out
        # as a string. Remove the brackets and commas, then split it into a list.
        if reduced_outformat == 'arr':
            date_out = date_out.translate(None, '[],').split()

        return date_out

    #################################
    @staticmethod
    def now(tz=None):
        """
        OUTPUT
           SptDatetime(datetime.now(tz=tz))
        """
        return SptDatetime(datetime.now(tz=tz))
    
    #################################
    @staticmethod
    def utcnow():
        """
        OUTPUT
           SptDatetime(datetime.utcnow())
        """
        return SptDatetime(datetime.utcnow())

    #################################
    def __getattr__(self, name):
        """ Trying to access an attribute which does not exist, but whose first three
        letters match a valid formatDate format, will give you the output of formatDate. """
        try:
            return self.formatDate(name)
        except ValueError:
            raise AttributeError("Object "+str(type(self))+" has no attribute "+name+"!")

    
    #################################
    def __str__(self):
        """ As a string, the SptDatetime object is represented by the default
        output format in the SptDatetime.formatDate() function. """
        return self.formatDate()

    def __repr__(self):
        """
        The SptDatetime constructor can parse any string output by the SptDatetime.formatDate function.
        """
        return self.__str__()
     
    #################################
    def __fitsformat__(self):
        """
        The sptpol_software.util.fits.addBintabToFits function will call this function to determine
        how to store this object in a FITS file. We'll store it as a modified Julian date. We 
        can easily re-create the SptDatetime object from such a number.
        """
        return self.formatDate('mjd')
    
    #################################
    def __hdf5format__(self, full_object=True):
        """
        The sptpol_software.util.table.writeSptTable function will call this function to determine
        how to store this object in an HDF5 file. We'll store it as a modified Julian date. We 
        can easily re-create the SptDatetime object from such a number.
        """
        if full_object:
            return {'__class__':class_name(self),
                    'date':self.formatDate('mjd')}
        else:
            return self.formatDate('mjd')
    
    #################################
    @classmethod
    def __unpack_hdf5__(cls, data):
        """
        Expect that the data are only a date in SptDatetime-readable format. Use that to create a new SptDatetime object.
        """
        return cls(data['date'])
     
    #################################
    def __add__(self, other):
        """
        Calls datetime's math function, then tries to convert the 
        result to an SptDatetime object. If that fails, return the
        result of the datetime operation.
        """
        result = super(SptDatetime, self).__add__(other)
        try: return SptDatetime(result)
        except TypeError: return result
    def __sub__(self, other):
        """
        Calls datetime's math function, then tries to convert the 
        result to an SptDatetime object. If that fails, return the
        result of the datetime operation.
        """
        result = super(SptDatetime, self).__sub__(other)

        # If we couldn't get a result via the datetime subtraction function,
        # then try converting the input argument into an SptDatetime object
        # and do it again.
        if result == NotImplemented:
            result = super(SptDatetime, self).__sub__( SptDatetime(other) )
        try: return SptDatetime(result)
        except TypeError: return result
    def __radd__(self, other):
        """
        Calls datetime's math function, then tries to convert the 
        result to an SptDatetime object. If that fails, return the
        result of the datetime operation.
        """
        result = super(SptDatetime, self).__radd__(other)
        try: return SptDatetime(result)
        except TypeError: return result
    def __rsub__(self, other):
        """
        Calls datetime's math function, then tries to convert the 
        result to an SptDatetime object. If that fails, return the
        result of the datetime operation.
        """
        result = super(SptDatetime, self).__rsub__(other)
        
        # If we couldn't get a result via the datetime subtraction function,
        # then try converting the input argument into an SptDatetime object
        # and do it again.
        if result == NotImplemented:
            result = super(SptDatetime, self).__rsub__( SptDatetime(other) )
        try: return SptDatetime(result)
        except TypeError: return result
        
    def __lt__(self, other):
        """
        First tries to call datetime's LT function. If that doesn't work,
        tries to return !(other <= self).
        
        AUTHOR
            Stephen Hoover, 29 October 2013
        """
        try:
            result = super(SptDatetime, self).__lt__(other)
        except TypeError:
            result = ~(other <= self)
        return result
        
    def __le__(self, other):
        """
        First tries to call datetime's LE function. If that doesn't work,
        tries to return !(other < self).
        
        AUTHOR
            Stephen Hoover, 29 October 2013
        """
        try:
            result = super(SptDatetime, self).__le__(other)
        except TypeError:
            result = ~(other < self)
        return result
        
    def __gt__(self, other):
        """
        First tries to call datetime's GT function. If that doesn't work,
        tries to return !(other >= self).
        
        AUTHOR
            Stephen Hoover, 29 October 2013
        """
        try:
            result = super(SptDatetime, self).__gt__(other)
        except TypeError:
            result = ~(other >= self)
        return result
        
    def __ge__(self, other):
        """
        First tries to call datetime's GE function. If that doesn't work,
        tries to return !(other > self).
        
        AUTHOR
            Stephen Hoover, 29 October 2013
        """
        try:
            result = super(SptDatetime, self).__ge__(other)
        except TypeError:
            result = ~(other > self)
        return result


    

