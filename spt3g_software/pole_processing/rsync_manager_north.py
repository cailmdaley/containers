#!/usr/bin/env python

from spt3g.std_processing import status_db as sdb
import argparse as ap
import os
import subprocess as sp
import shutil
import tempfile
import datetime
import glob
from spt3g import gcp
from spt3g.std_processing.transfer_tools import SFTPComm, update_verify_root
from spt3g.std_processing.transfer_tools import LockedProcess, LockedFileError

P = ap.ArgumentParser(description="Master script for rsync management up north.",
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)

P.add_argument('--database',
               default='/buffer/transfer_database/aux_transfer_north.txt',
               help='Database file for files transferred from the remote')
P.add_argument('--database-south',
               default='/buffer/transfer_database/aux_transfer.txt',
               help='Database file for aux files, rsynced from the south')
P.add_argument('--aux-data-root', default='/buffer',
               help='Default directory for aux output files')
P.add_argument('--rsync-data-root', default='/buffer/rsync',
               help='Directory for rsync data')
P.add_argument('--transfer-root', action='store', default='/buffer/3g_temp',
               help='Transfer data root')
P.add_argument('--verify-root', default='/buffer/3g_verified',
               help='Verified transfer data root')
P.add_argument('--copy-root', default='/spt/data',
               help='Copied data root')

args = P.parse_args()

# create transfer lock file
try:
    lock = LockedProcess(os.path.abspath(__file__), timeout=5)
except LockedFileError as e:
    raise SystemExit

print('*' * 80)
sdb.print_time('Starting rsync transfer manager')
print('*' * 80)

sftp = SFTPComm()

transfer_root = None

try:

    # load database with options to minimize lock time
    db_opts = dict(sync=True, verbose=True, user='spt',
                   host='spt-buffer.grid.uchicago.edu')
    db = sdb.AuxTransferDatabase(args.database, **db_opts)
    rdb = sdb.AuxTransferDatabase(args.database_south, **db_opts)

    # create temporary directory for file transfer
    if args.transfer_root:
        transfer_root = args.transfer_root
    else:
        transfer_root = tempfile.mkdtemp()
    sdb.print_time('Storing temporary transfer files in {}'.format(transfer_root))

    # create date-stamped verified files directory
    verify_root = update_verify_root(args.verify_root)

    # check that rsync'ed aux files have been updated recently
    now = float(datetime.datetime.now().strftime('%s'))

    rsync_entries = sorted(set(rdb.get_entries(type='rsync')) |
                           set(db.get_entries(type='rsync')))
    for entry in rsync_entries:
        if entry[-1] == '/':
            base_entry = os.path.basename(entry[:-1])
        else:
            base_entry = os.path.basename(entry)
        local_entry = os.path.join(args.rsync_data_root, base_entry)
        rsync_touch = os.path.join(
            args.rsync_data_root, os.path.splitext(
                sftp.local2remote(entry, aux=True, auxtype='rsync'))[0]
            + '.last_rsync')

        db.update(entry, type='rsync')

        if base_entry in ['transfer_database', 'transfer_logs']:
            sp.check_call(['touch', rsync_touch])

        try:
            mtime = os.stat(rsync_touch).st_mtime
        except OSError as e:
            sdb.print_warn(
                'Rsync timestamp file {} is missing'.format(rsync_touch))
            db.update(entry, status='error')
            continue
        else:
            mdate = datetime.datetime.fromtimestamp(mtime).strftime('%c')
            if (now - mtime) > 86400:
                sdb.print_warn(
                    'Rsync of {} is out of date, last synced {}'.format(
                        entry, mdate))
                db.update(entry, status='error')
                continue
            else:
                sdb.print_time(
                    'Rsync {} verified, last synced {}'.format(entry, mdate))
                db.update(entry, status='verified')

        # permanent storage paths
        copy_data_root = args.rsync_data_root.replace(
            args.aux_data_root, args.copy_root)
        copy_rsync_touch = rsync_touch.replace(
            args.rsync_data_root, copy_data_root)
        copy_entry = os.path.join(copy_data_root, base_entry)
        try:
            copy_mtime = os.stat(copy_rsync_touch).st_mtime
        except OSError as e:
            copy_mtime = 0

        # daily snapshot tarball
        tarfile_base = '{}.{}.tar.gz'.format(
            base_entry, datetime.datetime.now().strftime('%Y-%m-%d'))
        tarfile = os.path.join(transfer_root, tarfile_base)
        verified_tarfile = os.path.join(verify_root, tarfile_base)

        # check if entry has been updated recently
        if ((copy_mtime != mtime) or
            (base_entry in ['transfer_database', 'transfer_logs'])):
            # skip if up to date

            # rsync to permanent storage
            sdb.print_time(
                'Rsyncing {} to {}'.format(local_entry, copy_data_root))
            try:
                sp.check_call(
                    ['rsync', '-aviP'] +
                    ["--exclude='*.db-journal'"] * (base_entry == 'transfer_database') +
                    [local_entry, copy_data_root])
                sp.check_call(
                    ['rsync', '-aviP', rsync_touch, copy_data_root])
            except Exception as e:
                sdb.print_warn(
                    'Error rsyncing {} to {}'.format(
                        local_entry, copy_data_root), e)
                db.update(entry, status='copy_error')
                continue
            else:
                sdb.print_time(
                    'Rsynced {} to {}'.format(local_entry, copy_data_root))

        if ((copy_mtime == mtime) or
            (base_entry in ['transfer_database', 'transfer_logs'])) \
            and not os.path.exists(verified_tarfile):

            # compress and store snapshot, rsync to NFS, once daily
            sdb.print_time(
                'Archiving {} to {}'.format(local_entry, verified_tarfile))
            try:
                cwd = os.getcwd()
                os.chdir(args.rsync_data_root)
                sp.check_call(['tar', '-czf', tarfile, base_entry])
                shutil.move(tarfile, verified_tarfile)
                os.chdir(cwd)
            except Exception as e:
                sdb.print_warn(
                    'Error archiving {} to {}'.format(
                        local_entry, verified_tarfile), e)
                db.update(entry, status='copy_error')
            else:
                sdb.print_time('Archived {} to {}'.format(
                    local_entry, verified_tarfile))
                db.update(entry, status='copied')

except Exception as e:

    sdb.print_warn('Rsync transfer error', e)
    raise e

finally:

    # cleanup tmp
    if transfer_root and not args.transfer_root:
        shutil.rmtree(transfer_root)

    print('*' * 80)

    # remote transfer lock file
    lock.close()
