"""
This contains general functions that extend typical python/numpy operations on 
basic python data structures.
"""
import inspect, collections, functools, os, traceback, sys
import numpy as np
import scipy.stats
from copy import copy
from spt3g.core.funcconstruct import usefulfunc
from spt3g import core

@usefulfunc
def merge_dictionary(dst, src):
    '''
    Recursively updates dst with the values in src
    -Flat out stolen from an internet post by Manuel Muradas.  
    If you ever meet him, buy him a beer.
    '''

    stack = [(dst, src)]
    while stack:
        current_dst, current_src = stack.pop()
        for key in current_src:
            if key not in current_dst:
                current_dst[key] = current_src[key]
            else:
                if (isinstance(current_src[key], dict) and 
                        isinstance(current_dst[key], dict)):
                    stack.append((current_dst[key], current_src[key]))
                else:
                    current_dst[key] = current_src[key]
    return dst

@usefulfunc
def uniquify_list(seq, preserve_order=False):  
    '''
    Returns a list with all the elements of seq where all the duplicate 
    elements are removed
    '''
    if preserve_order:
        checked = [] 
        for e in seq: 
            if e not in checked: 
                checked.append(e) 
        return checked
        
    return list(set(seq))

@usefulfunc
def sort_two_lists(l1, l2):
    '''
    sorts two lists, where l1 controls the order
    XXX: python built-in sorted() supports this natively
    '''
    if len(l1)==0:
        return [], []
    l3 = list(zip(l1, l2))
    l3.sort()
    l3 = list(zip(*l3))
    l1 = list(l3[0])
    l2 = list(l3[1])
    return l1, l2

@usefulfunc
def sort_n_lists(*args):
    '''
    sorts n lists, where first list controls the order
    XXX: python built-in sorted() supports this natively
    '''
    if len(args[0])==0:
        return tuple([[]] * len(args))
    l3 = list(zip(*args))
    l3.sort()
    return list(zip(*l3))

@usefulfunc
def is_list_sorted(l):
    '''
    returns True if the list is sorted
    '''
    return all(l[i] <= l[i+1] for i in range(len(l)-1))

@usefulfunc
def has_only_unique_elements(seq):
    '''
    returns True if the list has only unique elements
    '''
    for i in range(len(seq)):
        for j in range(i+1, len(seq)):
            if seq[i] == seq[j]:
                return False
    return True

@usefulfunc
def shift(array, shift):
    """
    A wrapper for np.roll. You may specify a shift for each axis, and this
    function will perform each shift at once. This method of specifying the 
    shift is sometimes more convenient.

    Arguments
    ---------
    array : ndarray
        The array which you wish to shift.
    shift : iterable
        The shift of each axis of the input array. `shift`
        must have a number of elements equal to the number of dimensions 
        of `array`!

    Returns
    -------
    A shifted version of the input ndarray.
    """
    # Check the input array sizes.
    if len(array.shape) != len(shift):
        raise ValueError("You must specify a shift for each axis of the array!")

    output = array
    for axis, this_shift in enumerate(shift):
        if this_shift:
            # Only make this shift if it's non-zero.
            output = np.roll(output, this_shift, axis=axis)

    return output

@usefulfunc
def round_to_multiple(n, m):
    '''
    return n rounded to the nearest multiple of m
    '''
    m = float(m)
    if type(n)==type([]):
      return list((np.array(n)/m).round()*m)
    elif type(n)==type(1) or type(n)==type(1.):
      return round(n/m)*m
    else:
      raise Exception('Unknown type %s for rounding to multiple %f' % 
                      (str(type(n)), m))


@usefulfunc
def split_on_numbers(s):
    '''
    Splits the string into a list where the numbers and the characters between 
    numbers are each element
    '''
    prevDig = False
    outList = []
    for char in s:
        if char.isdigit():
            if prevDig:
                outList[-1] += char
            else:
                prevDig = True
                outList.append(char)
        else:
            if not prevDig and len(outList)>0:
                outList[-1] += char
            else:
                prevDig = False
                outList.append(char)

    return outList

@usefulfunc
def is_float(s):
    '''
    returns true if s can be cast to float
    '''
    try:
        float(s)
        return True
    except:
        return False

@usefulfunc
def is_int(s):
    '''
    returns true if s can be cast to integer
    '''
    try:
        int(s)
        return True
    except:
        return False

@usefulfunc
def str_cmp_with_numbers_sorted(str1, str2):
    '''
    Compares two strings where numbers are sorted according to value, 
    so Sq12 ends up after Sq8,  use in sorted function
    '''

    if str1==str2:
        return 0
    split1 = split_on_numbers(str1)
    split2 = split_on_numbers(str2)

    largestStr = 0
    for l in [split1, split2]:
        for s in l:
            if s[0].isdigit():
                largestStr = len(s) if len(s) > largestStr else largestStr

    for l in [split1, split2]:
        for i in range(len(l)):
            if l[i][0].isdigit():
                l[i] =  '0'*(largestStr-len(l[i])) +l[i]

    p1 = reduce(lambda x,y: x+y, split1)
    p2 = reduce(lambda x,y: x+y, split2)
    return -1 if p1<p2 else 1

@usefulfunc
class stdTracer:
    '''
    Instantiate this class and set sys.stdout to it if you want a traceback 
    to be printed whenever a print statement happens.  If you have undesired 
    output this is a way of tracking it down
    '''
    def write(self, msg):
        for line in traceback.format_stack():
            sys.__stdout__.write(line.strip())
        sys.__stdout__.write(msg)

    def flush(self):
        sys.__stdout__.flush()

@usefulfunc
def profile_function_decorator(fn):
    '''
    A decorator for functions you wish to have run information on for 
    optimizing and debugging purposes

    To use it: add @profile_function_decorator to the line above the 
    function declaration
    
    so:

    @profile_function_decorator
    def testcall(taco, fish, k = 'test'):

    will cause the testcall function to always spit out performance profile 
    information after it runs

    '''
    import cProfile
    def wrapped(*args, **kwargs):
        prof = cProfile.Profile()
        retval = prof.runcall(fn, *args, **kwargs)
        prof.print_stats()
        return retval
    return wrapped

@usefulfunc
def get_args_supplied():
    '''
    returns a dictionary with the arguments mapping to their values 
    of the function that is the level above this one
    python magic ho!
    '''
    frame = inspect.currentframe(1)
    args, _, _, values = inspect.getargvalues(frame)
    return values

@usefulfunc
def get_git_commit_version(root_path):
    '''
    returns the git commit hash tag for the git repository at root_path
    
    root_path must be the top level of the git repo
    '''
    import commands
    command_to_run = 'git --git-dir=%s/.git --work-tree=%s/ rev-parse HEAD'%(root_path, root_path)
    return commands.getoutput(command_to_run)

@usefulfunc
def get_git_commit_version_environ(eviron_var):
    '''
    returns the git commit hash tag for the git repository which has a top 
    level directory stored by environ_var
    
    environ_var is the shell's environmental variable
    '''
    import os
    return get_git_commit_version(os.environ[eviron_var])

@usefulfunc
def get_function_args(func):
    ass = inspect.getargspec(func)
    arg_dic = {}
    if ass.defaults != None:
        nkargs = len(ass.defaults)
        for i in range(-1,-1*(len(ass.defaults)+1),-1):
            arg_dic[ ass.args[i] ] = ass.defaults[i]
    return arg_dic

@usefulfunc
def gcd(a,b):
    '''
    returns the greatest common denominator of a and b
    '''
    if a < b:
        a, b = b, a
    if a == b:
        return a
    else:
        return gcd(a-b, b)

@usefulfunc
def gcd_set(s):
    '''
    returns the gcd of the list s

    if s is too long, this function will crash horribly
    I should rewrite it as not a recursive function, but who doesn't 
    love recursion?

    Answer: most people
    '''
    if len(s) == 0:
        raise RuntimeException('No GCD exists of empty set')
    elif len(s) == 1:
        return s[0]
    elif len(s) == 2:
        return gcd(s[0], s[1])
    return gcd( gcd_set( s[len(s)//2:]), gcd_set( s[:len(s)//2]) )

@usefulfunc
def find_closest_index(val, lst):
    '''
    lst and val contain numbers
    finds the index of the element in the list closest to val
    '''
    
    dist = abs(val-lst[0])
    ind = 0
    for i in range(1, len(lst)):
        dn = abs(lst[i]-val)
        if dn < dist:
            ind = i
            dist = dn
    return ind

class MemoizeMutable:
    '''
    Memoization decorator for functions that accept mutable objects
    '''
    def __init__(self, fn):
        self.fn = fn
        self.memo = {}
    def __call__(self, *args, **kwds):
        import cPickle
        str_hash = cPickle.dumps(args, 1)+cPickle.dumps(kwds, 1)
        if not str_hash in self.memo: 
            self.memo[str] = self.fn(*args, **kwds)
        return self.memo[str]


class MemoizedHashableNoKWArgs(object):
   '''
   Decorator. Caches a function's return value each time it is called.
   If called later with the same arguments, the cached value is returned
   (not reevaluated).
   '''
   def __init__(self, func):
      self.func = func
      self.cache = {}
   def __call__(self, *args):
      if not isinstance(args, collections.Hashable):
         # uncacheable. a list, for instance.
         # better to not cache than blow up.
         return self.func(*args)
      if args in self.cache:
         return self.cache[args]
      else:
         value = self.func(*args)
         self.cache[args] = value
         return value
   def __repr__(self):
      '''Return the function's docstring.'''
      return self.func.__doc__
   def __get__(self, obj, objtype):
      '''Support instance methods.'''
      return functools.partial(self.__call__, obj)


class DataPrinter:
    '''
    used as a callback for pylab plots to print out what you are clicking on
    

    #example usage:
    mouse = DataPrinter()
    pl.connect('button_press_event', mouse.mycall)
    x = np.arange(6)
    y = x**2
    y2 = x**3
    mouse.addDataSet('x sq', x, y)
    mouse.addDataSet('x cube', x, y2)
    pl.plot(x,y)
    pl.plot(x,y, '.')
    pl.plot(x,y2)
    pl.plot(x,y2, '.')

    pl.show()
    '''
    event = None
    def __init__(self, frac_screen = .03):
        self.data_set = {}
        self.frac_screen = frac_screen
    def mycall(self, event):
        import pylab as pl
        self.event = event
        if event.xdata == None or event.ydata == None:
            return

        ax = pl.gca()  # get current axis
        def valid_points(x,y):
            xsize = ax.get_xlim()[1]-ax.get_xlim()[0]
            ysize = ax.get_ylim()[1]-ax.get_ylim()[0]
            dist_perc =(((x-event.xdata)/xsize)**2 + ((y-event.ydata)/ysize)**2)**.5
            valid_inds = np.where( dist_perc < self.frac_screen)
            return len(valid_inds[0])>0
        for k in self.data_set:
            if valid_points(self.data_set[k][0],self.data_set[k][1]):
                print(k)
        print('')
    def addDataSet(self, name, x, y ):
        self.data_set[name] = (np.asarray(x),np.asarray(y))

@usefulfunc
def class_name(_object):
    """
    Given any object, this function will return its fully-specified 
    type/class name. This output is suitable for use in "tools.getClass". 
    """
    objname = str(type(_object))
    return objname[objname.find("'")+1 : objname.rfind("'")]

@usefulfunc
def validate_gz_files(file_name_list):
    '''
    returns a list of the files that are corrupt and extant
    '''
    bad_files = []
    for fn in file_name_list:
        print("Checking %s"%fn)
        if not os.path.exists(fn):
            print("Uhhhh %s doesn't exist."%fn)
        elif (os.system('gunzip -t %s' % fn)):
            bad_files.append(fn)
    return bad_files

@usefulfunc
def validate_gz_file( fn):
    '''
    returns the file name if it is corrupt
    '''
    bad_files = []
    print("Checking %s"%fn)
    if not os.path.exists(fn):
        print("Uhhhh %s doesn't exist."%fn)
        return None
    elif (os.system('gunzip -t %s' % fn)):
        return fn
    else:
        return None

@usefulfunc
def price_is_right_index(v, lst, check_sorted = False):
    '''
    Returns the index that is closest without going over.
    Assumes that lst is sorted, or checks if check_sorted is true.

    check_sorted is O(n) not O(lg(n)) so makes things slower, asymptotically

    Remember to spay or neuter your pets.
    '''
    if check_sorted:
        assert(is_list_sorted(lst))
    start_ind = 0
    stop_ind = len(lst)-1
    while (True):
        if stop_ind - start_ind <= 1:
            good_ind =  start_ind if v < lst[stop_ind] else stop_ind
            return good_ind
        guess_ind = start_ind + (stop_ind-start_ind)//2
        if lst[guess_ind] == v:
            return guess_ind
        if lst[guess_ind] < v:
            start_ind = guess_ind
        else:
            stop_ind = guess_ind


@usefulfunc
def downsample_2d_array(arr, ds_fac = 2 ):
    sh = arr.shape
    assert(sh[0] % ds_fac == 0)
    assert(sh[1] % ds_fac == 0)
    new_arr = np.zeros( ( sh[0] / ds_fac, sh[1]/ds_fac))
    for i in range(sh[0]/ds_fac):
        for j in range(sh[1]/ds_fac):
            for i_int in range(ds_fac):
                for j_int in range(ds_fac):
                    new_arr[i,j] += arr[ i * ds_fac + i_int, j * ds_fac + j_int]
    return new_arr

def get_corrupt_g3_files(flst):
    def is_bad(fn):
        print("Checking", fn)
        try:
            reader = core.G3Reader(fn)
            while (reader(None)):
                pass
            return False
        except:
            return True
    return filter(is_bad, flst)

def percent_disk_usage(path):
    '''
    Return disk usage fraction for a given path
    '''
    st = os.statvfs(path)
    total = st.f_blocks * st.f_frsize
    used = (st.f_blocks - st.f_bfree) * st.f_frsize
    return float(used)/float(total)

def binFunction(function_in, ell_in, ell_out=None, delta_ell=None,
                ell_bins=None, do_sum=False):
    """
    Given an array of input values and input coordinates for each of those
    values, average the input values into bins as specified.

    INPUTS
        function_in: (array) The array of values to be binned.

        ell_in: (array) An array of the same length as "function_in", this
            specifies the coordinates of each index of the function_in array.

        ell_out [None]: (array) The central value of each of the output bins.
            If ell_out is specified, you must also give delta_ell.

        delta_ell [None]: (float) The width of an output bin. Used with
            the "ell_out" input.

        ell_bins [None]: (2D array or list) An array or list of 2-element
            arrays or lists. The first element gives the lower edge of
            an output bin (inclusive), and the second element is the upper
            edge (also inclusive).

    OUTPUT
        An array with length equal to either len(ell_out) or len(ell_bins),
        depending on which is supplied. Each value of the output will be
        equal to the average of the input function_in values inside an bin.

    EXCEPTIONS
        ValueError if you input both a delta_ell and an ell_bins input.
        ValueError if function_in and ell_in have different lengths.

    """
    # Raise an error if the user supplies no binning information, or if the
    # user fills all binning information fields.
    if (ell_bins is not None) == (ell_out is not None and
                                  delta_ell is not None):
        raise ValueError("Please supply either an ell binning, or else " +
                         "central ells and an ell range, but not both!")

    # Check that function_in and ell_in are the same length.
    if len(function_in) != len(ell_in):
        raise ValueError("The input function and input coordinate arrays " +
                         "must have the same length!")

    #Make sure ell_in is an array
    ell_in = np.array(ell_in)

    binned_function = []
    if ell_bins is not None:
        for this_bin in ell_bins:
            # Find which elements of the input array fall within the bin.
            ind = (ell_in >= this_bin[0]) & (ell_in <= this_bin[1])

            # Average together all values of the function which fall in the bin
            if not do_sum:
                binned_function.append(np.mean(function_in[ind]))
            else:
                binned_function.append(np.sum(function_in[ind]))
    else:
        for e in ell_out:
            ind = (ell_in > e - delta_ell/2.0) & (ell_in < e + delta_ell/2.0)
            if not do_sum:
                binned_function.append(np.mean(function_in[ind]))
            else:
                binned_function.append(np.sum(function_in[ind]))

    # Convert the ouput from a list to an array and return it.
    binned_function = np.array(binned_function)

    return binned_function

@usefulfunc
def multiway_number_partition(input, partitions=2, method='greedy'):
    '''
    Partition a list of numbers into groups such that all the sums of each
    group are approximately equal. Not guaranteed to find perfect partitions.

    Parameters:
    -----------
    input: list
        The list of numbers to partition

    partitions: int
        The number of partitions desired

    method ['greedy']: str
        The method to use to calculate the partition.
        Currently supported options are 'greedy' for the greedy heuristic
        or 'KK' for the Karmarkar-Karp set differencing method.

    Returns:
    --------
    final_partitions: list
        A list of lists containing the desired partitions.

    References:
    -----------
    Most code copied from a UNCW CSC 380 class paper by
    Lisanti, Fernandez, and Hunt

    R.E. Korf, Artif. Intell., 106 (2) (1998), pp. 181-203
    '''
    input = copy(input)

    if not isinstance(input,list):
        input = list(input)

    if method == 'greedy':
        # Sort the numbers in decreasing order, and always place
        # the next unassigned number in the subset with the smallest
        # sum so far, until all numbers have been assigned.
        final_partitions = [[] for i in range(partitions)]
        input.sort(reverse=True)
        for i in range(partitions):
            final_partitions[i].insert(i,input[0])
            input.remove(input[0])
        ii = 0
        finished = False
        while not finished:
            index_list = []
            while ii < partitions:
                index_list.insert(ii, sum(final_partitions[ii]))
                ii += 1
            minVal = min(index_list)
            insert = index_list.index(minVal)
            final_partitions[insert].insert(0,input[0])
            input.remove(input[0])
            ii = 0
            if len(input) == 0:
                finished = True
        return final_partitions

    elif method == 'KK':
        sub_parts = []
        temp_list = []
        final_partitions = []
        f = 1
        count = 0
        input.sort(reverse=True)
        for num in input:
            if len(temp_list) != partitions:
                temp_list.append(num)
            else:
                sub_parts.append(temp_list)
                temp_list = []
                temp_list.append(num)
        sub_parts.append(temp_list)
        while len(sub_parts[len(sub_parts)-1]) < partitions:
            sub_parts[len(sub_parts)-1].append(0)
        while count != len(sub_parts)-1:
            temp_list = sub_parts[0]
            sub_parts.append(temp_list)
            sub_parts.remove(sub_parts[0])
            count += 1
        for n in sub_parts[0]:
            if n != 0:
                final_partitions.append([n])
            else:
                final_partitions.append([])
        sub_parts.remove(sub_parts[0])
        for sub in sub_parts:
            for num in sub:
                if num != 0:
                    final_partitions[partitions - f].append(num)
                    f += 1
            f = 1
        return final_partitions

    else:
        raise ValueError('%s not a supported option'%method)

