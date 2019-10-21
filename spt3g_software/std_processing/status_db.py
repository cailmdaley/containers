import pandas as pd
import numpy as np
import os
import time
import fcntl
import re
import warnings
import socket
import pwd
import sys
import traceback
import subprocess as sp
from contextlib import contextmanager

def print_time(msg):
    """
    Print timestamped log messages.
    """
    # timestamp prefix
    ts = '[ {:%c %Z} ] '.format(pd.Timestamp.now('utc'))

    # indent multi-line log message
    msg = '\n    '.join(str(msg).split('\n'))

    # print
    print('{}{}'.format(ts, msg))

def print_warn(msg, e=None):
    if isinstance(e, sp.CalledProcessError) and hasattr(e, 'output'):
        msg = '{}:\n{}'.format(msg, e.output.decode() if e.output else None)
    print_time(msg)
    if e is not None:
        tb = traceback.format_exc()
        print(tb)
        print(e)
        warnings.warn('{}:\n{}\n{}'.format(msg, tb, e))
    else:
        warnings.warn(msg)

def obs2date(obsid):
    return str(pd.Timestamp('2017-01-01', tz='utc') +
               pd.Timedelta(int(obsid), unit='s'))

# regex for matching on field scans
field_regex_groups = 'ra([0-9]+)h([0-9\.]*)dec([0-9\.-]+)'
field_regex = field_regex_groups.replace('(', '').replace(')', '')

class StatusMeta(type):
    """
    Modify docstrings for StatusDatabase subclasses.
    """

    def __new__(cls, name, parents, attrs):

        def argstr(cols, dtypes):
            return '\n        '.join(
                ['{} : {}'.format(c, np.dtype(attrs[dtypes][c]).name)
                 for c in attrs[cols]])

        idxstr = argstr('index_columns', 'index_dtypes')
        statstr = argstr('status_columns', 'status_dtypes')

        # update method
        docstr_update = """
        Update the database.

        Arguments
        ---------
        {}
            Index argument(s) that uniquely determine the database entry.

        Keyword Arguments
        -----------------
        {}
            Status arguments
        commit : bool
            If True, update the database on disk as well as in memory.
        """.format(idxstr, statstr)

        try:
            update = attrs['update']
            update.__doc__ = docstr_update

        except KeyError as e:
            def update(self, *args, **kwargs):
                return parents[0].update(self, *args, **kwargs)
            update.__doc__ = docstr_update
            attrs['update'] = update

        # match method
        docstr_match = """
        Return a bool, index array or DataFrame matching the given index and
        status values.

        Arguments
        ---------
        {}
            Index arguments.  Positional, but not all need be supplied.
            Can also be a list of sets of index arguments to match against
            several entries at once.

        Keyword Arguments
        -----------------
        {}
            Status arguments.  Can be a single value or a list of values.
        days : scalar
            If supplied, match all entries newer than the given number of days.
            Use floating point to represent partial days.
        return_df : bool
            If True, return the pandas.DataFrame object.
        return_index : bool
            If True, return a list of indices into the internal DataFrame
            that match the input criteria.  If both return_index and
            return_df are False, a boolean array is returned.
        squeeze : bool
            If True, and `return_df` is True, return the DataFrame object
            as a single row if only one match is found.

        Returns
        -------
        out : bool array, index array, or pandas.DataFrame object
        """.format(idxstr, statstr)

        try:
            match = attrs['match']
            match.__doc__ = docstr_match

        except KeyError as e:
            def match(self, *args, **kwargs):
                return parents[0].match(self, *args, **kwargs)
            match.__doc__ = docstr_match
            attrs['match'] = match

        # drop method
        docstr_drop = """
        Drop entries from the database that match the input index and status
        values.

        Arguments
        ---------
        {}
            Index arguments.  Positional, but not all need be supplied.

        Keyword Arguments
        -----------------
        {}
            Status arguments.  Can be a single value or a list of values.
        commit : bool
            If True, update the database on disk as well as in memory.
        """.format(idxstr, statstr)

        try:
            drop = attrs['drop']
            drop.__doc__ = docstr_drop

        except KeyError as e:
            def drop(self, *args, **kwargs):
                return parents[0].drop(self, *args, **kwargs)
            drop.__doc__ = docstr_drop
            attrs['drop'] = drop

        # get_entries method
        docstr_get = """
        Get a list of entries matching the requested status column values.

        Arguments
        ---------
        {}
            Index arguments.  Positional, but not all need be supplied.

        Keyword Arguments
        -----------------
        {}
            Status arguments.  Can be a single value or a list of values.
        return_time : bool
            If True, return the modified time of each matched entry

        Returns
        -------
        entries : list of tuples
            List of entries found.  Each entry is a tuple of index values
            corresponding to each of the index columns
        modified : array-like
            The last-modified timestamp for each of the returned entries,
            if `return_time` is True.
        """.format(idxstr, statstr)

        try:
            get_entries = attrs['get_entries']
            get_entries.__doc__ = docstr_get

        except KeyError as e:
            def get_entries(self, *args, **kwargs):
                return parents[0].get_entries(self, *args, **kwargs)
            get_entries.__doc__ = docstr_get
            attrs['get_entries'] = get_entries

        for k in attrs:
            attr = attrs[k]
            if not hasattr(attr, '__call__'):
                continue
            if not getattr(attr, '__doc__', None):
                continue
            attr.__doc__ = attr.__doc__.format(
                index_args=idxstr, status_args=statstr)

        return super(StatusMeta, cls).__new__(
            cls, name, parents, attrs)

from six import with_metaclass
class StatusDatabase(with_metaclass(StatusMeta, object)):
    """
    Generic database for storing status information on running jobs,
    e.g. scanification or file transfer.
    """

    # columns required for indexing into the database
    index_columns = ['item']

    # ordered list of columns to sort the database
    # if None, no sorting will be performed
    sort_columns = None

    # index data types
    index_dtypes = {'item': str}

    # docstrings
    index_docs = {'item': 'Database item'}

    # status columns, all optional
    status_columns = ['status']

    # status column data types
    status_dtypes = {'status': str}

    # docstrings
    status_docs = {'status': 'Status'}

    # default status values
    status_defaults = {}

    # possible status values
    status_values = {}

    def __init__(self, filename='db.txt', sync=False, verbose=False,
                 user='sptdaq', host='anal.spt', read_only=False,
                 match_regex=False):
        """
        Arguments
        ---------
        filename : str
            Database filename
        sync : bool
            If True, data base is closed after writing to release lock,
            and assumes commit=True on any update.  This maintains
            synchronization between the database in memory and on disk
            when edited by multiple processes.
        verbose : bool
            More words.
        user : string
            User that controls the file lock
        hostname : string
            Hostname on which the file is locked
        read_only : bool
            If True, read in the database and do not allow writing to it.
            The host/user check is ignored if this option is used.
            Files that are not writable are opened in read_only mode by
            default.
        match_regex : bool
            If True, allow regular expressions to match entries.
        """

        filename = self.get_filename(filename)
        if os.path.exists(filename) and not os.access(filename, os.W_OK):
            if not read_only:
                warnings.warn('Invalid permissions, opening file '
                              '{} in read-only mode'.format(filename))
            read_only = True

        self._read_only = read_only
        self._match_regex = match_regex

        if not self._read_only:
            # check user/host
            this_host = socket.gethostname()
            this_user = pwd.getpwuid(os.getuid())[0] # this works in cron jobs too
            if this_host != host:
                # raise RuntimeError('Database must be opened on host {}, found {}'.format(
                #     host, this_host))
                pass
            if this_user != user:
                raise RuntimeError('Database must be opened by user {}, found {}'.format(
                    user, this_user))

        # always have a modification time column
        if 'modified' not in self.status_columns:
            self.status_columns.append('modified')
        self.status_dtypes.setdefault('modified', pd.Timestamp)
        self.status_docs.setdefault('modified', 'Last modified')
        if self.sort_columns is None:
            self.sort_columns = ['modified']

        # merged database columns and dtype
        self.dtype = self.index_dtypes.copy()
        self.dtype.update(**self.status_dtypes)
        self.columns = self.index_columns + self.status_columns

        # load data from file ...
        self._filename = filename
        self._sync = sync
        self._verbose = verbose
        self._dirty = False
        self._dirty_sql = False
        self.log('Opening database {}.'.format(filename))
        self.open(filename)
        if self._sync or self._read_only:
            self.close()

    def __repr__(self):
        return repr(self.data)

    @staticmethod
    def get_filename(filename):
        """
        Try to automatically determine the database file root directory
        on commonly used machines.
        """
        if os.path.isabs(filename) or os.path.exists(filename):
            return filename

        db_root = {'scott': '/spt/data/transfer_database',
                   'amundsen': '/spt/data/transfer_database',
                   'spt-buffer': '/buffer/transfer_database',
                   'anal': '/data/autoproc/db'}

        host = socket.gethostname().split('.')[0]
        if host in db_root:
            return os.path.join(db_root[host], filename)

        return filename

    def log(self, msg):
        """
        More words
        """
        if self._verbose:
            print_time(msg)

    def open(self, filename=None, load=True):
        """
        Open and load the database file, using a file lock to control access.

        Arguments
        ---------
        filename : str
            If supplied, this filename is cached
        """

        if filename is None:
            filename = self._filename
        if not hasattr(self, '_file'):
            # self.log('Opening database {}.'.format(filename))

            if self._read_only:
                f = open(filename, 'r')
            else:
                f = open(filename, 'a+')
                try:
                    fcntl.lockf(f, fcntl.LOCK_EX)
                except IOError:
                    # wait and try again
                    time.sleep(5)
                    fcntl.lockf(f, fcntl.LOCK_EX)
            f.seek(0)

            # self.log('Opened database {}.'.format(filename))
            try:
                from pandas.errors import EmptyDataError
            except ImportError:
                from pandas.parser import EmptyDataError

            data = None
            if load:
                data = None
                try:
                    dtype = self.dtype.copy()
                    del dtype['modified']
                    converters = {'modified': self._get_date}
                    for key, dt in dtype.items():
                        if dt == pd.Timestamp:
                            converters[key] = self._get_date
                    for key in converters:
                        if key in dtype:
                            del dtype[key]

                    data = pd.read_csv(f, sep='\t', dtype=dtype, converters=converters)
                except EmptyDataError as e:
                    data = pd.DataFrame()
                    for k in self.columns:
                        data[k] = self._get_series(k, empty=True)
                finally:
                    if data is not None:
                        if len(data) or not hasattr(self, 'data'):
                            self.data = data
                    # self.log('Found {} entries'.format(len(data)))
                for k in self.status_columns:
                    if k not in data.columns:
                        data[k] = self._get_default(k, missing_column=True)
                        data[k] = self._get_series(k, data[k])
                        self._dirty = True
                        self._dirty_sql = True
                self._get_pathstrings(force=True)

            self._file = f
        return self._file

    def close(self):
        """
        Unlock and close the database file, flushing any pending changes to disk.
        """
        if not hasattr(self, '_file'):
            return
        if not self._read_only:
            self.write()
        f = self.open(load=False)
        if not self._read_only:
            fcntl.lockf(f, fcntl.LOCK_UN)
        f.close()
        delattr(self, '_file')
        # self.log('Closed database {}.'.format(self._filename))

    def write(self, force=False, filename=None):
        """
        Write the database to disk.

        Arguments
        ---------
        force : bool
            If True, write to disk even if there are no pending changes.
        filename : str
            If given, write the database to a new file with this name.
            Otherwise, the input filename is used.
        """
        if filename is not None:
            if os.path.abspath(filename) == os.path.abspath(self._filename):
                filename = None
        if filename is None:
            if self._read_only:
                raise IOError('Cannot write to read-only database {}.'.format(
                    self._filename))
            if not self._dirty and not force:
                return
            f = self.open(load=False)
            f.seek(0)
            f.truncate(0)
        else:
            f = open(filename, 'w')
        self.data.to_csv(f, index=False, sep='\t', columns=self.columns, mode='w')
        if filename is None:
            f.flush()
            self.log('Wrote {} entries to {}.'.format(len(self.data), self._filename))
            self._dirty = False
            if self._sync:
                self.close()
        else:
            f.close()
            self.log('Wrote {} entries to {}'.format(len(self.data), filename))

    @contextmanager
    def write_all(self, force=False):
        """
        Context for applying a sequence of database updates at once,
        without unnecessary file I/O between each update.

        Arguments
        ---------
        force : bool
            If True, write to disk even if there are no pending changes.

        Usage
        -----
        The following sequence of commands will load the database from disk,
        apply an update for each entry in the list, then commit the database
        to disk when all entries are complete.

            >>> with db.write_all():
            ...     for entry in entries:
            ...         db.update(*entry, status=None)
        """

        self.open()
        sync = self._sync

        try:
            self._sync = False
            yield
        finally:
            self._sync = sync
            self.write(force=force)

    # for completeness
    commit_all = write_all
    commit = write

    def _to_sql_frame(self):
        data = self.data.copy()
        if 'pathstring' in data.columns:
            data.drop('pathstring', axis=1, inplace=True)
        if 'modified' in data.columns:
            data['modified'] = data['modified'].apply(lambda d: str(d))
        for col in data.columns:
            if self.dtype[col] == pd.Timestamp:
                data[col] = data[col].apply(lambda d: str(d))
        return data

    def write_sql(self, filename=None, force=False):
        """
        Write data to SQL database

        Arguments
        ---------
        filename : string
            Filename where the database should be stored.
            If not supplied, this is stored in a file named identically
            to the input database, but with the `.db` extension.
        force : bool
            If True, force write the SQL database.
        """

        if filename is None:
            filename = os.path.splitext(self._filename)[0] + '.db'
        name = os.path.splitext(os.path.basename(filename))[0]

        # check if database file is out of date
        if not os.path.exists(filename):
            self._dirty_sql = True
        elif os.stat(self._filename).st_mtime > os.stat(filename).st_mtime:
            self._dirty_sql = True

        if not self._dirty_sql and not force:
            return

        import sqlalchemy
        conn = sqlalchemy.create_engine('sqlite:////tmp/{}.db'.format(name))

        # write database to disk
        data = self._to_sql_frame()
        data.to_sql(name, conn, if_exists='replace')

        # cleanup
        conn.dispose()
        sp.check_call('mv /tmp/{}.db {}'.format(name, filename).split())

        self._dirty_sql = False

    def _get_date(self, value):
        """
        Return a timestamp object in UTC
        """
        if isinstance(value, pd.Timestamp):
            return value
        return pd.Timestamp(value, tz='utc')

    def _get_default(self, key, missing_column=False):
        """
        Return default value for the given status key
        """
        if key == 'modified':
            if missing_column:
                return pd.Timestamp('2017-01-01', tz='utc')
            else:
                return pd.Timestamp.now('utc')
        else:
            ret = self.status_defaults.get(key, None)
            if missing_column:
                if self.status_dtypes[key] == str and (ret is None or np.isnan(ret)):
                    return ''
            return ret

    def _get_series(self, key, value=None, empty=False):
        """
        Return a pandas series of the correct type for the given value
        """
        if self.dtype[key] != pd.Timestamp:
            if empty:
                return pd.Series(dtype=self.dtype[key])
            return pd.Series(value, dtype=self.dtype[key])

        if empty:
            return pd.Series(pd.Timestamp('2017-01-01', tz='utc')).drop(0)
        if isinstance(value, pd.Series):
            return value.apply(self._get_date)

        return pd.Series(value)

    def _get_string_col(self, key, df=None):
        """
        Return a data column as a string Series.
        """
        if df is None:
            df = self.data
        if self.dtype[key] != pd.Timestamp:
            return df[key].astype(str)
        return df[key].apply(lambda x: str(x))

    def _get_pathstrings(self, force=False, df=None):
        """
        Update entry pathstring cache
        """
        if df is None:
            df = self.data
        if 'pathstring' in df and not force:
            return df['pathstring']
        ent = self._get_string_col(self.index_columns[0], df)
        for col in self.index_columns[1:]:
            ent += '/' + self._get_string_col(col, df)
        df['pathstring'] = ent
        return ent

    def _match_pathstring(self, *args, **kwargs):
        """
        Match an entry by pathstring.
        """
        df = kwargs.pop('df', None)
        entry = '/'.join(str(x) for x in args)
        match_regex = kwargs.pop('match_regex', self._match_regex)
        if df is None:
            df = self.data
        if match_regex:
            pattern = "^{}$".format(entry)
            return df['pathstring'].str.contains(pattern, regex=True)
        else:
            return df['pathstring'].isin([entry])

    def _check_args(self, *args, **kwargs):
        """
        Check values of index and status keyword arguments before processing.

        Arguments
        ---------
        args : index arguments, all required
        kwargs : status arguments, all optional

        Keyword Arguments
        -----------------
        check_len_args : bool
            If True (default), ensure that args contains the correct number
            of elements

        Returns
        -------
        args, kwargs
        """
        check_len_args = kwargs.pop('check_len_args', True)

        regex = False
        new_args = []
        if args:
            if isinstance(args[0], list):
                args = args[0]
            else:
                args = [args]
            for args1 in args:
                if not isinstance(args1, tuple):
                    args1 = (args1,)
                if len(args1) < len(self.index_columns):
                    if check_len_args:
                        raise ValueError('Malformed index {}'.format(args1))
                    else:
                        regex = True
                        args1 += ('[^/]+',) * (len(self.index_columns) - len(args1))
                elif len(args1) > len(self.index_columns):
                    raise ValueError('Malformed index {}'.format(args1))
                # coerce types
                new_args1 = ()
                for arg, col in zip(args1, self.index_columns):
                    if isinstance(arg, str):
                        if not check_len_args:
                            new_args1 += (arg,)
                        else:
                            # check that argument has no wildcards
                            if arg == '[^/]+':
                                raise ValueError('Malformed index {}'.format(args1))
                            new_args1 += (arg,)
                    else:
                        new_args1 += (self.index_dtypes[col](arg),)
                new_args += [new_args1]
            if not len(args):
                new_args = [()]
            elif len(new_args) == 1:
                new_args = new_args[0]

        for k in kwargs:
            if k not in self.status_columns:
                raise KeyError('Unrecognized key {}'.format(k))

            if k not in self.status_values:
                continue

            # allow None
            if kwargs[k] is None:
                continue

            v = kwargs[k]
            if not isinstance(v, list):
                single = True
                v = [v]
            else:
                single = False

            # make sure that all given values are allowed
            for idx, vv in enumerate(v):
                if vv is None:
                    continue
                try:
                    vlow = vv.lower() # ignore case
                except ValueError:
                    vlow = vv
                if vlow not in self.status_values[k]:
                    raise ValueError(
                        "Unrecognized {} value '{}'.  Choices are {}".format(
                            k, vv, self.status_values[k]))
                v[idx] = vlow

            kwargs[k] = v[0] if single else v

        return regex, new_args, kwargs

    def update(self, *args, **kwargs):

        if self._read_only:
            raise IOError('Cannot modify read-only database {}.'.format(
                self._filename))

        commit = kwargs.pop('commit', True if self._sync else False)
        match_regex = kwargs.pop('match_regex', self._match_regex)
        if commit:
            # reload from disk if sync is True
            f = self.open()
        regex, args, kwargs = self._check_args(*args, **kwargs)
        match_regex |= regex

        # match database index
        idx = self.match(*args, match_regex=match_regex)

        # logging
        entry = '/'.join(str(x) for x in args)
        values = '{' + ', '.join('{}={}'.format(k,v) for k,v in kwargs.items()) + '}'
        logstr = '{} entry {} with values {}'.format(
            'Adding' if not idx.any() else 'Updating', entry, values)
        sort_dirty = False

        if not idx.any():
            # create new frame
            f = pd.DataFrame()
            for k, v in zip(self.index_columns, args):
                f[k] = self._get_series(k, [v])
            for k in self.status_columns:
                f[k] = self._get_series(k, [kwargs.pop(k, self._get_default(k))])
            self._get_pathstrings(force=True, df=f)
            if len(self.data) == 0:
                self.data = f
            else:
                self.data = self.data.append(f, ignore_index=True)[self.columns + ['pathstring']]
                if self.sort_columns is not None and len(self.sort_columns):
                    sort_dirty = True
            self._dirty = True
            self._dirty_sql = True

        # update existing frame
        mtime = kwargs.get('modified', self._get_default('modified'))
        for k in self.status_columns:
            if k in kwargs:
                v = kwargs.pop(k)
                if not self.data[k][idx].isin([v]).all():
                    self.data.loc[idx, k] = v
                    self.data.loc[idx, 'modified'] = mtime
                    if k in self.sort_columns:
                        sort_dirty = True
                    self._dirty = True
                    self._dirty_sql = True

        # sort if necessary
        if sort_dirty:
            self.data.sort_values(by=self.sort_columns, inplace=True)

        # write to disk
        if not self._dirty and self._sync:
            self.close()
        elif self._dirty:
            self.log(logstr)
        if commit:
            self.write()

    def match(self, *args, **kwargs):

        return_df = kwargs.pop('return_df', False)
        return_idx = kwargs.pop('return_index', False)
        days = kwargs.pop('days', None)
        squeeze = kwargs.pop('squeeze', False)
        match_regex = kwargs.pop('match_regex', self._match_regex)

        kwargs.setdefault('check_len_args', False)
        regex, args, kwargs = self._check_args(*args, **kwargs)
        match_regex |= regex

        idx = pd.Series(np.ones(len(self.data), dtype=bool), index=self.data.index)

        if len(self.data) == 0:
            if return_df:
                return self.data
            if return_idx:
                return np.array([], dtype=int)
            return idx

        if args:
            # match any index tuples
            if not isinstance(args, list):
                args = [args]
            idx &= False
            for arg in args:
                if not len(arg):
                    continue

                if any([isinstance(x, int) for x in arg]):
                    df = self.data
                    # match all integer index columns first
                    for col, x in zip(self.index_columns, arg):
                        if not isinstance(x, int):
                            continue
                        df = df[df[col] == x]
                    # check all other columns
                    df = df[self._match_pathstring(*arg, df=df, match_regex=match_regex)]
                    idx[df.index] = True
                else:
                    idx |= self._match_pathstring(*arg, match_regex=match_regex)

        for k, v in kwargs.items():
            # match all keywords
            if not isinstance(v, list):
                if v is None:
                    idx &= self.data[k].isnull()
                else:
                    idx &= self.data[k] == v
            elif None in v:
                v = list(v)
                v.remove(None)
                idx &= self.data[k].isnull() | self.data[k].isin(v)
            else:
                idx &= self.data[k].isin(v)

        if days:
            # match all entries from the last `days` days
            days = float(days)
            ts = self._get_default('modified') - pd.Timedelta(days, unit='D')
            idx &= self.data['modified'] > ts

        if return_df:
            ret = self.data.loc[idx]
            if squeeze and len(ret) == 1:
                return ret.loc[ret.index[0]]
            return ret

        if return_idx:
            return np.where(idx)[0]

        return idx

    def drop(self, *args, **kwargs):

        if self._read_only:
            raise IOError('Cannot modify read-only database {}.'.format(
                self._filename))

        commit = kwargs.pop('commit', True if self._sync else False)
        if commit:
            f = self.open()

        kwargs.update(return_df=True)
        df = self.match(*args, **kwargs)
        if not len(df):
            if self._sync:
                self.close()
            return

        dropped_entries = self.get_entries(*args, **kwargs)

        self.data.drop(df.index, inplace=True)
        self._dirty = True
        self._dirty_sql = True

        self.log('Dropped entries {}'.format(dropped_entries))

        if commit:
            self.write()

    def get_entries(self, *args, **kwargs):

        return_time = kwargs.pop('return_time', False)
        kwargs['return_df'] = True
        df = self.match(*args, **kwargs)

        entries = [np.asarray(df[col]) for col in self.index_columns]
        if len(self.index_columns) == 1:
            ret = (entries[0],)
        else:
            ret = (list(zip(*entries)),)

        if return_time:
            ret += (np.asarray(df['modified']),)

        if len(ret) == 1:
            return ret[0]
        return ret


_status_values = ['error', 'processing', 'downsampling', 'buffering', 'buffered',
                  'stored', 'checking', 'ready', 'uploading', 'uploaded', 'sent',
                  'received', 'verified', 'copied', 'copy_error', 'scanify_error',
                  'storage_error', 'check_error']

class TransferDatabase(StatusDatabase):
    """
    Database of observations for managing downsampling and data transfer
    """

    index_columns = ['source', 'observation']
    index_dtypes = {'source': str,
                    'observation': int}
    index_docs = {'source': 'Observation source',
                  'observation': 'Observation ID'}
    sort_columns = ['observation']

    status_columns = ['status_fullrate', 'status_downsampled',
                      'transfer_fullrate', 'transfer_downsampled']
    status_dtypes = {'status_fullrate': str,
                     'status_downsampled': str,
                     'transfer_fullrate': int,
                     'transfer_downsampled': int}
    status_docs = {'status_fullrate': 'Status of fullrate data processing. Options: {}'.format(_status_values),
                   'status_downsampled': 'Status of downsampled data processing. Options: {}'.format(_status_values),
                   'transfer_fullrate': 'Transfer priority of fullrate data (0=no transfer, 1=primary, 2+=lower priority)',
                   'transfer_downsampled': 'Transfer priority of downsampled data (0=no transfer, 1=primary, 2+=lower priority)'}
    status_defaults = {'status_fullrate': np.NaN,
                       'status_downsampled': np.NaN,
                       'transfer_fullrate': 0,
                       'transfer_downsampled': 0}

    status_values = {'status_fullrate': _status_values,
                     'status_downsampled': _status_values}

    def __init__(self, filename='transfer.txt', **kwargs):

        super(TransferDatabase, self).__init__(filename=filename, **kwargs)

    def _to_sql_frame(self):
        data = super(TransferDatabase, self)._to_sql_frame()
        data['date'] = data['observation'].apply(obs2date)
        return data

_scanify_values = ['processing', 'success', 'error', 'rerun', 'permafail']

class ScanifyDatabase(TransferDatabase):
    """
    Database of observations for managing scanification, downsampling and data transfer
    """

    index_columns = [x for x in TransferDatabase.index_columns]
    index_dtypes = TransferDatabase.index_dtypes.copy()
    index_docs = TransferDatabase.index_docs.copy()
    sort_columns = [x for x in TransferDatabase.sort_columns]

    status_columns = TransferDatabase.status_columns + ['obs_start', 'obs_stop', 'status_scanify']
    status_dtypes = dict(list(TransferDatabase.status_dtypes.items()) +
                         [('obs_start', pd.Timestamp), ('obs_stop', pd.Timestamp),
                          ('status_scanify', str)])
    status_docs = dict(list(TransferDatabase.status_docs.items()) +
                       [('obs_start', 'Time of observation start'),
                        ('obs_stop', 'Time of observation stop'),
                        ('status_scanify', 'Scanifier status. Options: {}'.format(_scanify_values))])
    status_defaults = dict(list(TransferDatabase.status_defaults.items()) +
                           [('obs_start', pd.Timestamp('2017-01-01', tz='utc')),
                            ('obs_stop', pd.Timestamp('2017-01-01', tz='utc')),
                            ('status_scanify', np.NaN)])

    status_values = dict(list(TransferDatabase.status_values.items()) +
                         [('status_scanify', _scanify_values)])

    def __init__(self, filename='scanify.txt', **kwargs):

        super(ScanifyDatabase, self).__init__(filename=filename, **kwargs)

    def get_latest_start(self):
        if len(self.data) == 0:
            raise RuntimeError("Empty database")
        return self.data['obs_start'].max()

    def get_latest_stop(self):
        if len(self.data) == 0:
            raise RuntimeError("Empty database")
        return self.data['obs_stop'].max()

    def get_latest_contiguous_run(self):
        """
        Return the start time of the most recent incomplete observation in the
        database.  If all observations in the database are complete, return the
        stop time of the last observation in the database.

        An observation is considered complete if its fullrate processed data are
        in permanent storage, or if the observation has been marked as
        permanently failed.

        This method is used to determine the timestamp of the newest Timepoint
        frame that can be safely deleted from the temporary data buffer.
        """
        if len(self.data) == 0:
            raise RuntimeError("Empty database")
        idx = self.data['status_scanify'].isin(['success', 'permafail'])
        idx |= self.data['status_fullrate'] \
                   .isin(['ready', 'uploading', 'uploaded', 'sent', 'verified'])
        idx = np.logical_not(idx)
        start_times = self.data['obs_start'][idx]
        if len(start_times) == 0:
            return self.data['obs_stop'].max()
        start_times = start_times.apply(pd.Timestamp)
        last_stop = start_times.min()
        if last_stop is pd.NaT:
            # Somehow, entries without start or stop times
            # got into the database, and caused a couple days of
            # data loss.  This will just cause a hard failure instead.
            i = np.argmin(start_times)
            obsid = self.data['observation'][idx].iat[i]
            print(obsid)
            raise RuntimeError('Job database contains Not a Time (Observation ID: {})'.format(obsid))
        return last_stop

    def check_ready(self, *args, **kwargs):
        """
        Check whether the input observation is ready to be processed.

        Arguments
        ---------
        {index_args}
            Index arguments that uniquely determine the observation.

        Keyword Arguments
        -----------------
        rate : 'fullrate' or 'downsampled'
            The data rate of the observation to check.

        Returns
        -------
        True if the observation has been successfully scanified and verified,
        False otherwise.
        """
        rate = kwargs.pop('rate', 'fullrate')

        df = self.match(*args, return_df=True, squeeze=True, check_len_args=True)
        data_status = df['status_' + rate]

        if socket.gethostname().split('.')[0] == 'anal':
            # at pole
            return data_status in ['ready', 'uploading', 'uploaded', 'sent', 'verified']
        else:
            # north
            return data_status in ['verified']


_file_status_values = ['error', 'uploading', 'uploaded', 'sent', 'received',
                       'verified', 'copied', 'copy_error']

class TransferFileDatabase(StatusDatabase):
    """
    Database of files in each observation for managing data transfer
    """

    index_columns = ['rate', 'source', 'observation', 'file']
    index_dtypes = {'rate': str,
                    'source': str,
                    'observation': int,
                    'file': str}
    index_docs = {'rate': 'Data rate (fullrate or downsampled)',
                  'source': 'Observation source',
                  'observation': 'Observation ID',
                  'file': 'Individual file in the observation'}
    sort_columns = ['observation', 'rate', 'file']

    status_columns = ['size', 'priority', 'status', 'archive']
    status_dtypes = {'size': int,
                     'priority': int,
                     'status': str,
                     'archive': str}
    status_docs = {'size': 'File size',
                   'priority': 'Transfer priority (0=no transfer, 1=primary, 2+=lower priority)',
                   'status': 'Transfer status. Options: {}'.format(_file_status_values),
                   'archive': 'Transfer archive file'}
    status_defaults = {'size': 0,
                       'priority': 1}

    status_values = {'status': _file_status_values}

    def __init__(self, filename='transfer_files.txt', **kwargs):

        super(TransferFileDatabase, self).__init__(filename=filename, **kwargs)

    def _to_sql_frame(self):
        data = super(TransferFileDatabase, self)._to_sql_frame()
        data['date'] = data['observation'].apply(obs2date)
        return data


_aux_type_values = ['pydfmux', 'pydfmuxlog', 'arc', 'rsync', 'tar', 'eht', 'maser', 'map']
_aux_status_values = ['error', 'uploading', 'uploaded', 'sent', 'received',
                      'verified', 'copied', 'copy_error']

class AuxTransferDatabase(StatusDatabase):
    """
    Database of aux files for managing data transfer
    """

    index_columns = ['filename']
    index_dtypes = {'filename': str}
    index_docs = {'filename': 'Absolute path to data file or directory'}
    sort_columns = ['modified', 'type', 'filename']

    status_columns = ['type', 'size', 'status', 'archive']
    status_dtypes = {'type': str,
                     'size': int,
                     'status': str,
                     'archive': str}
    status_defaults = {'size': 0}
    status_docs = {'type': 'File data type.  Options: {}'.format(_aux_type_values),
                   'size': 'File size',
                   'status': 'Transfer status: Options: {}'.format(_aux_status_values),
                   'archive': 'Transfer archive file'}

    status_values = {'status': _aux_status_values,
                     'type': _aux_type_values}

    def __init__(self, filename='aux_transfer.txt', **kwargs):

        super(AuxTransferDatabase, self).__init__(filename=filename, **kwargs)

    def _to_sql_frame(self):
        data = super(AuxTransferDatabase, self)._to_sql_frame()

        def calc_date(entry):
            tp = entry['type']
            fn = entry['filename']
            if tp == 'pydfmux':
                fn_split = os.path.split(fn)
                if fn_split[0]:
                    datestr = '_'.join(fn_split[1].split('_', 2)[:2])
                    fmt = '%Y%m%d_%H%M%S'
                else:
                    datestr = fn_split[1]
                    fmt = '%Y%m%d'
            elif tp == 'pydfmuxlog':
                datestr = fn.split('.')[-1]
                fmt = '%Y-%m-%d'
            elif tp == 'arc':
                datestr = fn.split('.', 1)[0]
                fmt = '%Y%m%d_%H%M%S'
            else:
                return str(pd.Timestamp(entry['modified']))

            try:
                return str(pd.Timestamp.strptime(datestr + ' UTC', fmt + ' %Z'))
            except:
                return str(pd.Timestamp(entry['modified']))

        data['date'] = data.apply(calc_date, axis=1)
        return data

class JobDatabase(StatusDatabase):
    """
    Database of scanification jobs
    """

    index_columns = ['jobid']
    index_dtypes = {'jobid': int}
    index_docs = {'jobid': 'Scanification condor job number'}
    sort_columns = ['jobid']

    status_columns = ['start', 'stop', 'finished', 'success']
    status_dtypes = {'finished': bool,
                     'success': bool,
                    'start': pd.Timestamp,
                    'stop': pd.Timestamp}
    status_docs = {'finished': 'Job completed',
                   'success': 'Job completed successfully',
                   'start': 'Data start time',
                   'stop': 'Data stop time'}
    status_defaults = {'finished': False,
                       'success': False,
                       'start': np.NaN,
                       'stop': np.NaN}

    def __init__(self, filename='job.txt', **kwargs):

        super(JobDatabase, self).__init__(filename=filename, **kwargs)

    def get_latest_contiguous_run(self):
        start_times = self.data['start'][np.bitwise_not(self.data['success'])]
        start_times = [pd.Timestamp(t, tz = 'UTC') for t in start_times]
        if len(start_times) == 0:
            return max([pd.Timestamp(t) for t in self.data['stop']])
        last_stop = min(start_times)
        if last_stop is pd.Timestamp(None):
            # Somehow, entries without start or stop times
            # got into the database, and caused a couple days of
            # data loss.  This will just cause a hard failure instead.
            i = np.argmin(start_times)
            jobid = self.data['jobid'][np.bitwise_not(self.data['success'])].iat[i]
            print(jobid)
            raise RuntimeError('Job database contains Not a Time (job id: {})'.format(jobid))
        return last_stop

    def update(self, *args, **kwargs):

        if not self.match(*args).any():
            if kwargs.get('start', None) is None:
                raise ValueError('start value required')
            if kwargs.get('stop', None) is None:
                raise ValueError('stop value required')

        return super(JobDatabase, self).update(*args, **kwargs)

    def _to_sql_frame(self):
        data = super(JobDatabase, self)._to_sql_frame()
        for k in ['start', 'stop']:
            data[k] = data[k].apply(lambda d: str(d))

class BumpyStorageDatabase(StatusDatabase):
    """
    Database of storage locations for each observation
    """

    index_columns = ['source', 'observation']
    index_dtypes = {'source': str,
                    'observation': int}
    index_docs = {'source': 'Observation source',
                  'observation': 'Observation ID'}
    sort_columns = ['observation']

    status_columns = ['disk']
    status_dtypes = {'disk': str}
    status_docs = {'disk': 'Disk on which the observation is stored'}

    def __init__(self, filename='bumpy_storage.txt', **kwargs):

        super(BumpyStorageDatabase, self).__init__(filename=filename, **kwargs)

_proc_status_values = ['scheduled', 'processing', 'complete', 'failed', 'submitted',
                       'executing', 'stopped', 'permafail', 'permafail-hk']

class AutoprocDatabase(StatusDatabase):
    """
    Database of observation autoprocessing jobs
    """

    index_columns = ['source', 'observation']
    index_dtypes = {'source': str,
                    'observation': int}
    index_docs = {'source': 'Observation source',
                  'observation': 'Observation ID'}
    sort_columns=['observation']

    status_columns = ['status', 'job_id']
    status_dtypes = {'status': str, 'job_id': int}
    status_defaults = {'status': np.NaN, 'job_id': 0}
    status_docs = {'status': 'Processing status. Options: {}'.format(_proc_status_values),
                   'job_id': 'Condor job ID'}

    status_values = {'status': _proc_status_values}

    def __init__(self, filename='autoproc.txt', **kwargs):

        super(AutoprocDatabase, self).__init__(filename=filename, **kwargs)

