#!/usr/bin/env python

from spt3g.std_processing import status_db as sdb
import argparse as ap
import numpy as np
import os
import subprocess as sp
import shutil
import glob
import tempfile
import itertools
import sys
import time

from spt3g.std_processing.transfer_tools import RsyncComm, RsyncError
from spt3g.std_processing.transfer_tools import LockedProcess, LockedFileError

P = ap.ArgumentParser(description="Master script for transfer management.",
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)

P.add_argument('--database', default='/data/autoproc/db/scanify.txt',
               help='Transfer management database file')
P.add_argument('--database-remote', default='/buffer/transfer_database/scanify.txt',
               help='Transfer management database file, northern copy')
P.add_argument('--file-database',
               default='/data/autoproc/db/transfer_files.txt',
               help='Database file for files transfered to the remote')
P.add_argument('--file-database-remote',
               default='/buffer/transfer_database/transfer_files.txt',
               help='Database file for files transfered, rsynced copy form the south')
P.add_argument('--file-database-north',
               default='/data/autoproc/db/transfer_north_files.txt',
               help='Database file for files received north')
P.add_argument('--file-database-north-remote',
               default='/buffer/transfer_database/transfer_north_files.txt',
               help='Remote location for file_database_north')
P.add_argument('--aux-database', default='/data/autoproc/db/aux_transfer.txt',
               help='Database file for aux file transfers')
P.add_argument('--aux-database-remote',
               default='/buffer/transfer_database/aux_transfer.txt',
               help='Database file for aux transfers, rsynced copy from the south')
P.add_argument('--aux-database-north',
               default='/data/autoproc/db/aux_transfer_north.txt',
               help='Database file for aux files received north')
P.add_argument('--aux-database-north-remote',
               default='/buffer/transfer_database/aux_transfer_north.txt',
               help='Database file for aux files received north')
P.add_argument('--rsync-user', default='spt',
               help='Username on north machine')
P.add_argument('--rsync-host', default='spt-buffer.grid.uchicago.edu',
               help='Hostname of north machine')
P.add_argument('--rsync-identity', default='~/.ssh/id_rsa_uc',
               help='SSH secure identity file')
P.add_argument('--no-rsync', dest='rsync', action='store_false', default=True,
               help='Rsync data from/to UC machine')
P.add_argument('--logdir', default='/data/autoproc/logs/transfer',
               help='Logfile directory')
P.add_argument('--logdir-remote', default='/buffer/transfer_logs',
               help='Remote logfile directory')

args = P.parse_args()

# verification lock file
try:
    lock = LockedProcess(os.path.abspath(__file__), timeout=5)
except LockedFileError as e:
    raise SystemExit

print('*' * 80)
sdb.print_time('Starting verification')
print('*' * 80)

try:

    if args.rsync:
        rsync = RsyncComm(args.rsync_user, args.rsync_host,
                          identity=args.rsync_identity)
        if not rsync.isup():
            # check if the satellite is up to synchronize databases
            sdb.print_time('No rsync connection available, exiting.')
            raise SystemExit

        # rsync verified file database from the north
        sdb.print_time('Rsyncing northern file database')
        rsync.get_file(args.file_database_north_remote, args.file_database_north, lock=True)
        rsync.get_file(args.aux_database_north_remote, args.aux_database_north, lock=True)

    # load all the databases
    db_opts = dict(verbose=True)
    db = sdb.ScanifyDatabase(args.database, **db_opts)
    fdb = sdb.TransferFileDatabase(args.file_database, **db_opts)
    rfdb = sdb.TransferFileDatabase(args.file_database_north, **db_opts)
    adb = sdb.AuxTransferDatabase(args.aux_database, **db_opts)
    radb = sdb.AuxTransferDatabase(args.aux_database_north, **db_opts)

    # merge verified files
    sdb.print_time('Updating verified files')
    verified_entries = rfdb.get_entries(status=['verified', 'copied'])
    with fdb.commit_all():
        for entry in verified_entries:
            if fdb.match(*entry, status='sent').any():
                fdb.update(*entry, status='verified', archive=None)

    verified_entries = radb.get_entries(status=['verified', 'copied'])
    with adb.commit_all():
        for entry in verified_entries:
            if adb.match(entry, status='sent').any():
                adb.update(entry, status='verified', archive=None)

    # merge error files
    sdb.print_time('Resetting files with verification errors')
    error_entries = rfdb.get_entries(status='error')
    with fdb.commit_all():
        for entry in error_entries:
            if fdb.match(*entry, status=['sent', 'verified', 'error']).any():
                fdb.update(*entry, status='error')

    error_entries = radb.get_entries(status='error')
    with adb.commit_all():
        for entry in error_entries:
            if adb.match(entry, status=['sent', 'verified', 'error']).any():
                adb.update(entry, status=None)

    # check all sent entries
    for rate in ['downsampled', 'fullrate']:
        status_key = 'status_{}'.format(rate)

        sent_entries = db.get_entries(**{status_key: ['sent', 'verified']})
        for src, obs in sent_entries:
            files = fdb.match(rate, src, obs, return_df=True)

            # already processed
            if not len(files):
                continue

            # all files in the observation have been verified
            if (files['status'] == 'verified').all():
                sdb.print_time('Verified observation {}/{}'.format(src, obs))
                db.update(src, obs, **{status_key: 'verified'})
                fdb.drop(rate, src, obs)

            elif (files['status'] == 'error').any():
                # some files have errors
                sdb.print_time('Observation {}/{} failed verification'.format(src, obs))
                sdb.print_time(files)
                for entry in fdb.get_entries(rate, src, obs, status='error'):
                    fdb.update(*entry, status=None)
                db.update(src, obs, **{status_key: 'ready'})

    db.close()
    fdb.close()
    rfdb.close()
    adb.close()
    radb.close()

    if args.rsync:
        # short wait to make sure file locks are released
        time.sleep(0.1)

        # rsync a copy of the transfer database to the north
        sdb.print_time('Rsyncing local transfer database to the north')

        # rsync database files
        dbdir = os.path.dirname(args.database)
        dbdir_backup = '/poleanalysis/sptdaq/db'
        dbdir_remote = os.path.dirname(args.database_remote)

        sp.check_call(['rsync', '-av', '/spt_data/bolodata/bumpy_storage.txt',
                       os.path.join(dbdir, 'bumpy_storage.txt')])

        # backup, don't complain about known issues here
        proc = sp.Popen(['rsync', '-av', dbdir, dbdir_backup],
                        stdout=sp.PIPE, stderr=sp.STDOUT)
        out = proc.communicate()[0].decode()
        retcode = proc.poll()
        if retcode:
            if 'read errors mapping' not in out:
                sdb.print_time(out)

        files = [args.file_database,
                 args.database,
                 args.aux_database,
                 os.path.join(dbdir, 'bumpy_storage.txt'),
                 os.path.join(dbdir, 'autoproc.txt'),
                 os.path.join(dbdir, 'onlinemaps.txt')]
        remote_files = [args.file_database_remote,
                        args.database_remote,
                        args.aux_database_remote,
                        os.path.join(dbdir_remote, 'bumpy_storage.txt'),
                        os.path.join(dbdir_remote, 'autoproc.txt'),
                        os.path.join(dbdir_remote, 'onlinemaps.txt')]
        for f, rf in zip(files, remote_files):
            try:
                rsync.send_file(f, rf, lock=True)
            except LockedFileError:
                pass

        # rsync log files
        log_files = ['transfer.log',
                     'aux_transfer.log',
                     'rsync.log',
                     'verify.log']
        for f in log_files:
            rsync.send_file(os.path.join(args.logdir, f), args.logdir_remote)

        rsync.get_file(os.path.join(args.logdir_remote, 'transfer_north.log'), args.logdir)
        rsync.get_file(os.path.join(args.logdir_remote, 'aux_transfer_north_copy.log'),
                       os.path.join(args.logdir, 'aux_transfer_north_copy.log'))

except RsyncError as e:

    sdb.print_time('RsyncError: {}'.format(e))

except Exception as e:

    if isinstance(e, IOError) and ("Resource deadlock avoided" in str(e)):
        sdb.print_time(e)
    else:
        sdb.print_warn('Verification error', e)

finally:

    print('*' * 80)

    # remote verification lock file
    lock.close()
