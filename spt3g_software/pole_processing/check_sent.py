#!/usr/bin/env python

import os
import sys
import argparse as ap
import socket
from spt3g.std_processing import status_db as sdb
import pandas as pd
import subprocess as sp
import shlex
from spt3g.std_processing.transfer_tools import SFTPComm
from spt3g.std_processing.transfer_tools import LockedProcess, LockedFileError
from check_tdrs import get_last_start, get_now

hostname = socket.gethostname()
if hostname.startswith('spt-buffer'):
    db_root = '/buffer/transfer_database'
    db_user = 'spt'
    sftp_host = 'sftp.usap.gov'
    sftp_identity='~/.ssh/A-379-S-2'
    hours = None
    north = True
elif hostname.startswith('anal'):
    db_root = '/data/autoproc/db'
    db_user = 'sptdaq'
    sftp_host = 'spfs.southpole.usap.gov'
    sftp_identity = '~/.ssh/id_dsa_spfs'
    hours = 24.
    north = False
elif hostname.startswith('scott.grid') or hostname.startswith('amundsen.grid'):
    db_root = '/spt/data/transfer_database'
    db_user = None
    sftp_host = None
    sftp_identity = None
    hours = None
    north = True

P = ap.ArgumentParser(description="Print transfer summary statistics.",
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('hours', default=hours, type=float, nargs='?',
               help='All files modified up to this many hours before the '
               'current time will be counted')
P.add_argument('-f', '--full', default=False, action='store_true',
               help='Print list of files')
P.add_argument('--database-root', default=db_root, help='Database file root')
P.add_argument('--database-user', default=db_user, help='Database user')
P.add_argument('--sftp-host', default=sftp_host, help='SFTP hostname')
P.add_argument('--sftp-identity', default=sftp_identity, help='SFTP identity')
args = P.parse_args()

# create lock file
try:
    lock = LockedProcess(os.path.abspath(__file__), timeout=5)
except LockedFileError as e:
    raise SystemExit

if args.hours is None:
    pass_str = 'last pass at '
    start = get_last_start('utc')
    args.hours = (get_now('utc') - start).total_seconds() / 3600.
    if args.hours > 30:
        args.hours = 24
        pass_str = ''
        start = get_now('utc') - pd.Timedelta(hours=args.hours)
else:
    pass_str = ''
    start = get_now('utc') - pd.Timedelta(hours=args.hours)
days = float(args.hours) / 24.

sftp = SFTPComm(user='A-379-S', host=args.sftp_host,
                identity=args.sftp_identity, verbose=False)

args.database_root = os.path.realpath(args.database_root)
db_opts = dict(sync=True, verbose=False, user=args.database_user,
               read_only=True)

db = sdb.ScanifyDatabase(
    os.path.join(args.database_root, 'scanify.txt'), **db_opts)
fdb = sdb.TransferFileDatabase(
    os.path.join(args.database_root, 'transfer_files.txt'), **db_opts)
adb = sdb.AuxTransferDatabase(
    os.path.join(args.database_root, 'aux_transfer.txt'), **db_opts)
rfdb = sdb.TransferFileDatabase(
    os.path.join(args.database_root, 'transfer_north_files.txt'), **db_opts)
radb = sdb.AuxTransferDatabase(
    os.path.join(args.database_root, 'aux_transfer_north.txt'), **db_opts)

sdb.print_time('Transfer status since {}{:%c %Z}.'.format(pass_str, start))
out = sp.check_output(shlex.split(
    'ls -l --time-style="+%b %d %H:%M %Z" {}'.format(args.database_root)),
                      env={'TZ': 'UTC'}).decode()
sdb.print_time('Databases last updated:\n{}'.format(out))

def print_summary(stype, sz, num):
    print('  {:>8s}: {:7.3f} GB  ({} files)'.format(
        stype, float(sz) / 1024**3, num))

def print_entry(status, sz, modified, fname):
    print('    {:8s}{:12d}{:>25s}  {}'.format(
        status, sz, modified, fname))

def print_status(statuses, action, north_statuses=None, days=None,
                 ignore=None, ignore_north=None, ignore_days=None):
    print(action)

    ignore_entries = []
    ignore_obs = dict(fullrate=[], downsampled=[])
    if ignore is not None:
        ignore_entries += fdb.get_entries(status=ignore, days=ignore_days)
        for rate in ['fullrate', 'downsampled']:
            ignore_obs[rate] = db.get_entries(
                days=ignore_days, **{'status_{}'.format(rate): ignore})
    if ignore_north is not None:
        ignore_entries += rfdb.get_entries(status=ignore_north, days=ignore_days)

    def print_obs(fdb, statuses):
        obs_entries = fdb.match(status=statuses, days=days, return_df=True)
        num_obs_entries = len(obs_entries)
        sz_obs = obs_entries['size'].sum()
        if num_obs_entries:
            if args.full:
                print('  Obs data files:')
            for idx, df in obs_entries.iterrows():
                entry = tuple(df[x] for x in fdb.index_columns)
                sz = df['size']
                if (df['source'], df['observation']) in ignore_obs[df['rate']]:
                    sz_obs -= sz
                    num_obs_entries -= 1
                    continue
                if entry in ignore_entries:
                    sz_obs -= sz
                    num_obs_entries -= 1
                    continue
                if not args.full:
                    continue
                priority = df['priority']
                status = str(df['status'])
                if status == 'nan':
                    status = 'waiting'
                modified = '{:%Y-%m-%d %H:%M:%S %Z}'.format(df['modified'])
                fname = sftp.local2remote(entry, priority=priority)
                print_entry(status, sz, modified, fname)
        return num_obs_entries, sz_obs

    num_obs_entries, sz_obs = print_obs(fdb, statuses)
    if north_statuses is not None:
        num, sz = print_obs(rfdb, north_statuses)
        num_obs_entries += num
        sz_obs += sz
    print_summary('Obs data', sz_obs, num_obs_entries)

    ignore_entries = []
    if ignore is not None:
        ignore_entries += list(adb.get_entries(status=ignore, days=ignore_days))
    if ignore_north is not None:
        ignore_entries += list(radb.get_entries(status=ignore_north, days=ignore_days))

    def print_aux(adb, statuses):
        aux_entries = adb.match(status=statuses, days=days, return_df=True)
        aux_entries = aux_entries.loc[aux_entries['type'] != 'rsync']
        num_aux_entries = len(aux_entries)
        sz_aux = aux_entries['size'].sum()
        if num_aux_entries:
            if args.full:
                print('  Aux data files:')
            for idx, df in aux_entries.iterrows():
                entry = df['filename']
                sz = df['size']
                if entry in ignore_entries:
                    sz_aux -= sz
                    num_aux_entries -= 1
                    continue
                if not args.full:
                    continue
                ftype = df['type']
                status = str(df['status'])
                if status == 'nan':
                    status = 'waiting'
                modified = '{:%Y-%m-%d %H:%M:%S %Z}'.format(df['modified'])
                fname = sftp.local2remote(entry, aux=True, auxtype=ftype)
                print_entry(status, sz, modified, fname)
        return num_aux_entries, sz_aux

    num_aux_entries, sz_aux = print_aux(adb, statuses)
    if north_statuses is not None:
        num, sz = print_aux(radb, north_statuses)
        num_aux_entries += num
        sz_aux += sz
    print_summary('Aux data', sz_aux, num_aux_entries)
    print_summary('Total', sz_obs + sz_aux, num_obs_entries + num_aux_entries)

def print_remote(status, action):
    print(action)

    f_obs, sz_obs = sftp.list_files('*p[0-9]_*.tar', by_time=True, return_size=True)
    f_aux, sz_aux = sftp.list_files('*aux_*.tar', by_time=True, return_size=True)
    if args.full:
        print('  Obs data files:')
        for f, sz in zip(f_obs, sz_obs):
            print_entry(status, sz, '', f)
    print_summary('Obs data', sum(sz_obs), len(sz_obs))
    if args.full:
        print('  Aux data files:')
        for f, sz in zip(f_aux, sz_aux):
            print_entry(status, sz, '', f)
    print_summary('Aux data', sum(sz_aux), len(sz_aux))
    sz_tot = list(sz_obs) + list(sz_aux)
    print_summary('Total', sum(sz_tot), len(sz_tot))

ignore_north = ['received', 'verified', 'copied']
print_status([None, 'uploading'], 'Waiting at pole:', ignore_north=ignore_north)
print_status('uploaded', 'Queued at pole:', ignore_north=ignore_north)
if not north and args.sftp_host:
    print_remote('queued', 'Queued at {}:'.format(args.sftp_host))
print_status('sent', 'Sent from pole:', ignore_north=ignore_north)
if north and args.sftp_host:
    print_remote('remote', 'Received at {}:'.format(args.sftp_host))
print_status('received', 'Received at north:', 'received')
print_status('copied', 'Verified at north:', 'copied',
             ignore='verified', ignore_days=days)
print_status('verified', 'Verified at pole:', days=days)
