#!/usr/bin/env python

from spt3g.std_processing import status_db as sdb
import argparse as ap
import subprocess as sp
import os
from spt3g.std_processing.transfer_tools import LockedProcess, LockedFileError

P = ap.ArgumentParser(description="Master script for transfer management up north.",
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)

P.add_argument('--file-database',
               default='/buffer/transfer_database/transfer_north_files.txt',
               help='Database file for files transferred from the remote')
P.add_argument('--file-database-south',
               default='/buffer/transfer_database/transfer_files.txt',
               help='Database file for files transfered from remote '
               'rsynced from the south')
P.add_argument('--database-south',
               default='/buffer/transfer_database/scanify.txt',
               help='Transfer database rsynced from the south')
P.add_argument('--aux-database',
               default='/buffer/transfer_database/aux_transfer_north.txt',
               help='Database file for files transferred from the remote')
P.add_argument('--aux-database-south',
               default='/buffer/transfer_database/aux_transfer.txt',
               help='Database file for aux files, rsynced from the south')

args = P.parse_args()

# create transfer lock file
try:
    lock = LockedProcess(os.path.abspath(__file__), timeout=5)
except LockedFileError as e:
    raise SystemExit

print('*' * 80)
sdb.print_time('Starting verification manager')
print('*' * 80)

transfer_root = None

try:

    # load database with options to minimize lock time
    db_opts = dict(sync=True, verbose=True, user='spt', match_regex=True,
                   host='spt-buffer.grid.uchicago.edu')
    db = sdb.ScanifyDatabase(args.database_south, **db_opts)
    fdb = sdb.TransferFileDatabase(args.file_database, **db_opts)
    rfdb = sdb.TransferFileDatabase(args.file_database_south, **db_opts)

    # remove files that have been verified in the south
    sdb.print_time('Clearing verified observations')

    verified_entries = fdb.get_entries(status=['received', 'copied'])
    drop_entries = []
    for rate, src, obs, f in verified_entries:
        entry = (rate, src, obs)
        if entry in drop_entries:
            continue
        status_key = 'status_{}'.format(rate)
        if not db.match(src, obs, **{status_key: 'verified'}).any():
            continue
        files = fdb.match(rate, src, obs, return_df=True)
        if not len(files):
            continue
        if files['status'].isin(['copied']).all():
            drop_entries.append(entry)
        elif files['status'].isin(['received']).all():
            if not rfdb.match(rate, src, obs).any():
                drop_entries.append(entry)
    with fdb.commit_all():
        for entry in drop_entries:
            fdb.drop(*entry)

    # clear re-sent error entries
    sdb.print_time('Clearning re-sent error observations')
    error_entries, error_times = fdb.get_entries(status='error', return_time=True)
    for e, t in zip(error_entries, error_times):
        re, rt = rfdb.get_entries(
            *e, status=[None, 'uploaded', 'uploading', 'sent'], return_time=True)
        # remove if file was uploaded/sent after being marked as error here
        if len(re) and rt > t:
            fdb.drop(*e)

    adb = sdb.AuxTransferDatabase(args.aux_database, **db_opts)
    radb = sdb.AuxTransferDatabase(args.aux_database_south, **db_opts)

    # clear verified entries from northern DB
    sdb.print_time('Clearing verified aux entries')
    verified_entries = radb.get_entries(status='verified')
    copied_entries = adb.get_entries(status='copied')
    rsync_entries = adb.get_entries(type='rsync')
    done_entries = list(set(copied_entries) & set(verified_entries)
                        - set(rsync_entries))
    if done_entries:
        adb.drop(done_entries)

    # clear uploaded entries from northern DB
    # if present in the database, these are likely re-uploads due to error
    sdb.print_time('Clearing re-sent error aux entries')
    error_entries, error_times = adb.get_entries(status='error', return_time=True)
    for e, t in zip(error_entries, error_times):
        if e in rsync_entries:
            continue
        re, rt = radb.get_entries(
            e, status=[None, 'uploaded', 'uploading', 'sent'], return_time=True)
        # remove if file was uploaded/sent after being marked as error here
        if len(re) and rt > t:
            adb.drop(e)

    # store SQL database files
    db.write_sql()
    radb.write_sql()

    # backup
    sp.check_call("rsync -aviP /spt/user/production/autoproc.txt /buffer/rsync/transfer_database/autoproc_north.txt".split())
    sp.check_call("rsync -aviP --exclude='*.db-journal' /buffer/rsync/transfer_database /spt/data/rsync".split())

except Exception as e:

    sdb.print_warn('Verification error', e)
    raise e

finally:

    print('*' * 80)

    # remote transfer lock file
    lock.close()
