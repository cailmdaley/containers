#!/usr/bin/env python

import numpy as np
import pandas as pd
import subprocess as sp
import shlex
import os
import sys
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

localtz = None
_df = None
_dftz = None

def get_local_tz():
    global localtz
    if localtz is None:
        try:
            import tzlocal
            localtz = tzlocal.get_localzone()
        except ImportError:
            print('Cannot determine local timezone, using UTC')
            localtz = 'utc'
    return localtz

def get_tz(tz=None):
    if tz is None:
        tz = get_local_tz()
    if str(tz).lower() == 'pole':
        tz = 'Pacific/Auckland'
    elif str(tz).lower() == 'chicago':
        tz = 'US/Central'
    return tz

def get_now(tz=None):
    tz = get_tz(tz)
    return pd.Timestamp.now(tz)

# parse start time string
def fmt_start(x):
    ts = pd.Timestamp.strptime(x, '%y/%j/%H%M%S').tz_localize('utc')
    return ts

# parse duration
def fmt_duration(x):
    x = str(x)
    hms = ':'.join([xx if xx else '0' for xx in [x[:-4], x[-4:-2], x[-2:]]])
    ts = pd.Timedelta(hms)
    return ts

def load_tdrs_data(tz=None, force=False):
    global _df, _dftz

    tz = get_tz(tz)

    if not force and _df is not None:
        if str(_dftz) == str(tz):
            return _df

    # paths
    data_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'tdrs_data')
    sched_file = 'https://www.usap.gov/USAPgov/technology/documents/cel.txt'
    sched_local = os.path.join(data_dir, 'tdrs_schedule.txt')
    sched_pickle = os.path.join(data_dir, 'tdrs_schedule.py{}.pkl'.format(sys.version_info[0]))
    if not os.path.isdir(data_dir):
        os.makedirs(data_dir)

    now = get_now(tz)

    # load pickle file from disk if present
    update = force
    if not os.path.exists(sched_pickle):
        update = True
        df0 = None
    else:
        df0 = pd.read_pickle(sched_pickle)
        if df0['stop'].max() < now + pd.Timedelta(days=2):
            update = True

    # only try to update file from remote once a day
    if update and os.path.exists(sched_local):
        mtime = os.stat(sched_local).st_mtime
        if pd.Timestamp(mtime, unit='s', tz=tz) > now - pd.Timedelta(days=1):
            update = force

    # only update if a network connection is available
    if update:
        try:
            out = sp.check_output(
                shlex.split('ping -c 1 8.8.8.8'), stderr=sp.PIPE)
        except sp.CalledProcessError as e:
            print('No network connection available, schedule possibly out of date')
            update = force
        else:
            # download the schedule file
            print('Updating schedule file from remote...')
            sp.check_call(
                ['wget', sched_file, '-r', '-q', '-O', sched_local])

    if update and os.path.exists(sched_local):
        print('Parsing schedule file...')

        # strip metadata from text file
        with open(sched_local, 'r') as f:
            data = f.readlines()

            # strip header
            while not data[0].strip().startswith('---'):
                data = data[1:]
            data = data[1:]

            # strip footer
            while 'END OF REPORT' not in data[-1]:
                data = data[:-1]
            data = data[:-1]

            # strip intermediate headers
            empty = [not x.strip() for x in data]
            while any(empty):
                header_first = np.where(empty)[0][0]
                header_last = np.where(
                    [x.strip().startswith('---') for x in data[header_first:]])[0]
                if not len(header_last):
                    header_last = len(data) - 1
                else:
                    header_last = header_last[0] + header_first
                data = data[:header_first] + data[header_last + 1:]
                empty = [not x.strip() for x in data]

            # join
            data = ''.join(data)

        # parse data into frame
        df = pd.read_csv(
            StringIO(data), delim_whitespace=True,
            header=None,  usecols=[3, 5], engine='python',
            names=['c0', 'c1', 'c2', 'start', 'stop', 'duration',
                   'c6', 'c7', 'c8', 'c9', 'c10', 'c11', 'c12'])

        # construct parsed schedule
        start = df['start'].apply(fmt_start)
        duration = df['duration'].apply(fmt_duration)
        stop = start + duration

        df = pd.DataFrame({'start': start, 'stop': stop, 'duration': duration},
                          columns=['start', 'stop', 'duration'])

        # merge latest data from schedule
        if df0 is not None:
            df = pd.concat([df0.loc[df0['stop'] < df['start'].min()], df])
            df['start'] = df['start'].apply(lambda x: x.astimezone('utc'))
            df['stop'] = df['stop'].apply(lambda x: x.astimezone('utc'))
            df.reset_index(inplace=True, drop=True)
        df.to_pickle(sched_pickle)

    else:
        if df0 is not None:
            df = df0
        else:
            raise RuntimeError('No schedule data available')

    df['start'] = df['start'].apply(lambda x: x.astimezone(tz))
    df['stop'] = df['stop'].apply(lambda x: x.astimezone(tz))
    _df = df
    _dftz = tz
    return df

def get_passes(tz=None, force=False):
    if tz is None:
        tz = get_local_tz()
    df = load_tdrs_data(tz=tz, force=force)
    now = get_now(tz)

    nsel = df['stop'] > now
    psel = df['stop'] < now
    day_splits = df['start'].diff() > pd.Timedelta(hours=12)

    if not psel.any():
        idx = 0
        prev_idx = 0
    else:
        idx1 = (psel & day_splits).nonzero()[0]
        if len(idx1):
            prev_idx = idx1[-1]
        else:
            prev_idx = 0

    if not nsel.any():
        idx = len(df)
        next_idx = idx + 1
    else:
        idx = nsel.values.argmax()
        idx1 = (nsel & day_splits).nonzero()[0]
        if len(idx1) > 1:
            next_idx = idx1[1]
        elif len(idx1) == 1:
            next_idx = idx1[0]
        else:
            next_idx = idx + 1

    df_prev = df.loc[prev_idx:idx-1] if psel.any() else None
    df_next = df.loc[idx:next_idx-1] if nsel.any() else None

    return df_prev, df_next

def get_last_start(tz=None, force=False):

    df_prev, _ = get_passes(tz, force)
    if df_prev is None:
        return
    return df_prev.iloc[0]['start']

def get_last_stop(tz=None, force=False):

    df_prev, _ = get_passes(tz, force)
    if df_prev is None:
        return
    return df_prev.iloc[-1]['stop']

def get_next_start(tz=None, force=False):

    _, df_next = get_passes(tz, force)
    if df_next is None:
        return
    return df_next.iloc[0]['start']

def get_next_stop(tz=None, force=False):

    _, df_next = get_passes(tz, force)
    if df_next is None:
        return
    return df_next.iloc[-1]['stop']

if __name__ == "__main__":

    import argparse as ap

    P = ap.ArgumentParser(description="Update and print TDRS schedule",
                          formatter_class=ap.ArgumentDefaultsHelpFormatter)
    P.add_argument('--tz', default=get_local_tz(), help='Schedule timezone')
    P.add_argument('-f', '--force', default=False, action='store_true',
                   help='Force update of database')
    args = P.parse_args()

    df_prev, df_next = get_passes(tz=args.tz, force=args.force)

    print('Current time: {:%c %Z}\n'.format(get_now(args.tz)))
    print('Previous pass:\n{}'.format(df_prev))
    if df_prev is not None:
        print('\nTotal duration: {}\n'.format(df_prev['duration'].sum()))
    else:
        print('')
    print('Next pass:\n{}'.format(df_next))
    if df_next is not None:
        print('\nTotal duration: {}'.format(df_next['duration'].sum()))
