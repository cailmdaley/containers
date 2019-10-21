#!/usr/bin/env python

from spt3g.std_processing import status_db as sdb
import argparse as ap
import os
import shutil
import tempfile
from spt3g.std_processing.transfer_tools import SFTPComm, SFTPError
from spt3g.std_processing.transfer_tools import LockedProcess, LockedFileError

P = ap.ArgumentParser(description="Master script for transfer management up north.",
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)

P.add_argument('--file-database',
               default='/buffer/transfer_database/transfer_north_files.txt',
               help='Database file for files transferred from the remote')
P.add_argument('--aux-database',
               default='/buffer/transfer_database/aux_transfer_north.txt',
               help='Database file for aux files transferred from the remote')
P.add_argument('--transfer-root', action='store', default='/buffer/3g_temp',
               help='Transfer data root')
P.add_argument('--sftp-user', action='store', default='A-379-S',
               help='SFTP server username')
P.add_argument('--sftp-host', action='store', default='sftp.usap.gov',
               help='SFTP server hostname')
P.add_argument('--sftp-identity', action='store', default='~/.ssh/A-379-S-2',
               help='SFTP secure identity file')
P.add_argument('--maxdata', default=20, type=int,
               help='Default amount of data in GB to download at once from the server.')

args = P.parse_args()

# create transfer lock file
try:
    lock = LockedProcess(os.path.abspath(__file__), timeout=5)
except LockedFileError as e:
    raise SystemExit

print('*' * 80)
sdb.print_time('Starting transfer manager')
print('*' * 80)

# create SFTP communication object
sftp = SFTPComm(user=args.sftp_user, host=args.sftp_host,
                identity=args.sftp_identity)

transfer_root = None

try:

    # create temporary directory for file transfer
    if args.transfer_root:
        transfer_root = args.transfer_root
    else:
        transfer_root = tempfile.mkdtemp()
    sdb.print_time('Storing temporary transfer files in {}'.format(transfer_root))

    # look for files on the remote filesystem
    remote_files, remote_sizes = sftp.list_files('*.tar', by_time=True, return_size=True)
    sdb.print_time("Found {} files on remote".format(len(remote_files)))

    remote_list = [(f, sz) for f, sz in zip(remote_files, remote_sizes)
                   if not f.startswith('corrupt')]

    if not len(remote_files):
        raise SystemExit

    # load database with options to minimize lock time
    db_opts = dict(verbose=True, user='spt', sync=True,
                   host='spt-buffer.grid.uchicago.edu')
    fdb = sdb.TransferFileDatabase(args.file_database, **db_opts)
    adb = sdb.AuxTransferDatabase(args.aux_database, **db_opts)

    total = 0

    # process one file at a time
    for rf, rsz in remote_list:

        if total >= args.maxdata * 1024**3:
            sdb.print_time('Downloaded {} GB, stopping for now'.format(float(total) / 1024**3))
            break

        if not rf.startswith('bundle_'):
            try:
                entry, tag = sftp.remote2local(rf, return_entry=True, return_tag=True)
            except ValueError:
                continue
            aux = (not isinstance(entry, tuple)) or len(entry) == 1

            if aux and tag == 'rsync':
                continue

            if aux:
                edata = adb.match(entry, return_df=True, squeeze=True)
            else:
                edata = fdb.match(*entry, return_df=True, squeeze=True)

            # Skip files already marked as errored
            if len(edata) and edata['status'] == 'error':
                sdb.print_time("Skipping errored file {}".format(rf))
                continue

        sdb.print_time('Processing {}'.format(rf))

        tarfile = os.path.join(transfer_root, rf)
        try:
            files, sizes, bundle = sftp.download_file(rf, tarfile, delete=True)
        except SFTPError as e:
            sdb.print_warn(str(e))
            continue
        except Exception as e:
            sdb.print_warn(str(e))
            if rf.startswith('bundle_'):
                continue
            status = 'error'
            files, sizes, bundle = [rf], [rsz], None
        else:
            status = 'received'

        for f, sz in zip(files, sizes):

            try:
                entry, tag = sftp.remote2local(f, return_entry=True, return_tag=True)
            except ValueError:
                continue
            aux = (not isinstance(entry, tuple)) or len(entry) == 1

            if aux and tag == 'rsync':
                continue

            if status == 'error':
                sdb.print_warn("Error receiving file {}".format(f))
            else:
                sdb.print_time('File {} received successfully'.format(f))
            if aux:
                adb.update(entry, status=status, type=tag, size=sz, archive=bundle)
            else:
                fdb.update(*entry, status=status, priority=tag, size=sz, archive=bundle)
            total += sz
            sdb.print_time('Downloaded {} GB so far'.format(float(total) / 1024**3))

    fdb.close()
    adb.close()

except Exception as e:

    sdb.print_warn('Transfer error', e)
    raise e

finally:

    # cleanup tmp
    if transfer_root and not args.transfer_root:
        shutil.rmtree(transfer_root)

    print('*' * 80)

    # remote transfer lock file
    lock.close()
