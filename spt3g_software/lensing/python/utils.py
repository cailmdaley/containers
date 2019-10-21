import sys, time, copy, os, io
import numpy as np
import pickle as pk
import hashlib
from spt3g import core
from contextlib import contextmanager

def hash_check(hash1, hash2, ignore=[], keychain=[]):
    """
    Check key values between two pickeled hash files
    """
    keys1 = list(hash1.keys())
    keys2 = list(hash2.keys())

    for key in ignore:
        if key in keys1:
            keys1.remove(key)
        if key in keys2:
            keys2.remove(key)

    for key in set(keys1).union(set(keys2)):
        if key == "idf_def":
            v1 = {
                x: [z.rsplit("/", 1)[-1] for z in y]
                for x, y in hash1[key].items()
            }
            v2 = {
                x: [z.rsplit("/", 1)[-1] for z in y]
                for x, y in hash2[key].items()
            }
        elif key == "idfs":
            v1 = [x.rsplit("/", 1)[-1] for x in hash1[key]]
            v2 = [x.rsplit("/", 1)[-1] for x in hash2[key]]
        else:
            v1 = hash1[key]
            v2 = hash2[key]

        def hashfail(msg=None):
            core.log_warn(
                "ERROR: HASHCHECK FAIL AT KEY = "
                + ":".join(keychain + [str(key)])
            , unit="Utils")
            if msg != None:
                core.log_notice("   ", msg, unit="Utils")
            core.log_notice("   ", "V1 = ", v1, unit="Utils")
            core.log_notice("   ", "V2 = ", v2, unit="Utils")
            assert 0

        if type(v1) != type(v2):
            hashfail("UNEQUAL TYPES")
        elif type(v2) == dict:
            hash_check(v1, v2, ignore=ignore, keychain=keychain + [str(key)])
        elif type(v1) == np.ndarray:
            if not np.allclose(v1, v2):
                hashfail("UNEQUAL ARRAY")
        else:
            if not (v1 == v2):
                hashfail("UNEQUAL VALUES")


def cache_pk(suffix="", trim_lm=True):
    """decorator which can be used to cache the output of
    a function using pickle.
    """
    #suffix and trim_lm arguments
    #here are for backwards compatibility and can be
    #removed at some point. ditto for the format of 'tfname'

    def cache_lm_func(f):
        def cachelm(self, *args, **kwargs):
            fname = f.__name__
            if (
                (trim_lm == True)
                and (len(fname) > 3)
                and (fname[-3:] == "_lm")
            ):
                fname = fname[0:-3]
            tfname = self.lib_dir.format(prefix="temp") + (
                "/cache_lm_%s%s_%s.pk"
                % (
                    suffix,
                    hashlib.sha1(
                        np.ascontiguousarray(args + list(kwargs.items())[0])
                    ).hexdigest(),
                    fname,
                )
            )
            if not os.path.exists(tfname):
                core.log_notice("caching lm: %s" % tfname)
                print(
                    (
                        "lib_dir" + " " +
                        self.lib_dir + " " +
                        self.lib_dir.format(prefix="temp"),
                    )
                )
                lm = f(self, *args, **kwargs)
                pk.dump(lm, open(tfname, "wb"))
            try:
                core.log_notice(tfname)
                return pk.load(open(tfname, "rb"))
            except:
                core.log_notice("caching lm: %s" % tfname)
                print(
                    (
                        "lib_dir" + " " +
                        self.lib_dir + " " +
                        self.lib_dir.format(prefix="temp"),
                    )
                )
                lm = f(self, *args, **kwargs)
                pk.dump(lm, open(tfname, "wb"))
                #core.log_notice("wrap_cachelm:: failed to load ", tfname)
                #assert 0

        return cachelm

    return cache_lm_func


def enumerate_progress(list, label="", clear=False):
    """ implementation of python's enumerate built-in which
    prints a progress bar as it yields elements of a list. """
    t0 = time.time()
    ni = len(list)
    for i, v in enumerate(list):
        yield i, v
        ppct = int(100.0 * (i - 1) / ni)
        cpct = int(100.0 * (i + 0) / ni)
        if cpct > ppct:
            dt = time.time() - t0
            dh = np.floor(dt / 3600.0)
            dm = np.floor(np.mod(dt, 3600.0) / 60.0)
            ds = np.floor(np.mod(dt, 60))
            sys.stdout.write(
                "\r ["
                + ("%02d:%02d:%02d" % (dh, dm, ds))
                + "] "
                + label
                + " "
                + int(10.0 * cpct / 100) * "-"
                + "> "
                + ("%02d" % cpct)
                + r"%"
            )
            sys.stdout.flush()
    if clear == True:
        sys.stdout.write("\r")
        sys.stdout.write("\033[K")
        sys.stdout.flush()
    else:
        sys.stdout.write("\n")
        sys.stdout.flush()


class SumObjects(object):
    """ helper class to contain the sum of a number of objets. """

    def __init__(self, clone=copy.deepcopy):
        self.__dict__["__clone"] = clone
        self.__dict__["__sum"] = None

    def __iadd__(self, obj):
        self.add(obj)
        return self

    def __getattr__(self, attr):
        return getattr(self.__dict__["__sum"], attr)

    def __setattr__(self, attr, val):
        setattr(self.__dict__["__sum"], attr, val)

    def add(self, obj):
        if self.__dict__["__sum"] is None:
            self.__dict__["__sum"] = self.__dict__["__clone"](obj)
        else:
            self.__dict__["__sum"] += obj

    def get(self):
        return self.__dict__["__sum"]


class AverageObjects(object):
    """ helper class to contain the average of a number of objets. """

    def __init__(self, clone=copy.deepcopy):
        self.__dict__["__clone"] = clone
        self.__dict__["__sum"] = None
        self.__dict__["__num"] = 0

    def __iadd__(self, obj):
        self.add(obj)
        return self

    def __getattr__(self, attr):
        return getattr(self.get(), attr)

    def add(self, obj):
        if self.__dict__["__sum"] is None:
            self.__dict__["__sum"] = self.__dict__["__clone"](obj)
        else:
            self.__dict__["__sum"] += obj

        self.__dict__["__num"] += 1

    def get(self):
        return self.__dict__["__sum"] / self.__dict__["__num"]


class TimeDifference(object):
    """ helper class to contain / print a time difference. """

    def __init__(self, _dt):
        self.dt = _dt

    def __str__(self):
        return "%02d:%02d:%02d" % (
            np.floor(self.dt / 60 / 60),
            np.floor(np.mod(self.dt, 60 * 60) / 60),
            np.floor(np.mod(self.dt, 60)),
        )

    def __int__(self):
        return int(self.dt)


class DeltaTime(object):
    ''' helper class to contain / print a time difference. '''

    def __init__(self, _dt):
        self.dt = _dt

    def __str__(self):
        return ('%02d:%02d:%02d'% (np.floor(self.dt / 60 / 60),
                                 np.floor(np.mod(self.dt, 60*60) / 60 ),
                                 np.floor(np.mod(self.dt, 60))))
    def __int__(self):
        return int(self.dt)


class StopWatch(object):
    """ simple stopwatch timer class. """

    def __init__(self):
        self.st = time.time()
        self.lt = self.st

    def lap(self):
        """ return the time since start and the time since last call to lap or elapsed. """

        lt = time.time()
        ret = (DeltaTime(lt - self.st), DeltaTime(lt - self.lt))
        self.lt = lt
        return ret

    def elapsed(self):
        """ return the time since initialization. """

        lt = time.time()
        ret = DeltaTime(lt - self.st)
        self.lt = lt
        return ret

# -- below cribbed from
# http://stackoverflow.com/questions/5081657/how-do-i-prevent-a-c-shared-library-to-print-on-stdout-in-python

@contextmanager
def stdout_redirected(to=os.devnull):

    try:
        __IPYTHON__
    except NameError:
        pass
    else:
        to = sys.stdout  # do not redirect in python

    if to is not sys.stdout:
        fd = sys.stdout.fileno()

        def _redirect_stdout(nto):
            sys.stdout.close()  # + implicit flush()
            os.dup2(nto.fileno(), fd)  # fd writes to 'to' file
            sys.stdout = os.fdopen(fd, "w")  # Python writes to fd

        with os.fdopen(os.dup(fd), "w") as old_stdout:
            with open(to, "w") as file:
                _redirect_stdout(nto=file)
            try:
                yield  # allow code to be run with the redirected stdout
            finally:
                _redirect_stdout(nto=old_stdout)  # restore stdout.
                # buffering and flags such as
                # CLOEXEC may be different
    else:
        sys.stdout.flush()
        try:
            yield
        finally:
            pass
