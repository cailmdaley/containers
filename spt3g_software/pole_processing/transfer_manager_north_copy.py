#!/usr/bin/env python

from spt3g.std_processing import status_db as sdb
import argparse as ap
import os
import subprocess as sp
import shutil
import tempfile
import datetime
import numpy as np
import glob
from spt3g.std_processing.transfer_tools import SFTPComm, update_verify_root
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
P.add_argument('--data-root', default='/buffer/bolodata',
               help='Directory for scanified observations')
P.add_argument('--transfer-root', action='store', default='/buffer/3g_temp',
               help='Transfer data root')
P.add_argument('--sftp-user', action='store', default='A-379-S',
               help='SFTP server username')
P.add_argument('--sftp-host', action='store', default='sftp.usap.gov',
               help='SFTP server hostname')
P.add_argument('--sftp-identity', action='store', default='~/.ssh/A-379-S-2',
               help='SFTP secure identity file')
P.add_argument('--no-verify', action='store_false', default=True, dest='verify',
               help='Verify data on extraction')
P.add_argument('--delete', action='store_true', default=False,
               help='Delete tarball after unpacking')
P.add_argument('--verify-root', default='/buffer/3g_verified',
               help='Verified transfer data root')
P.add_argument('--copy-root', default='/spt/data/bolodata',
               help='Copied data root')
P.add_argument('--offline-cal-root', default='/spt/user/production/calibration/calframe',
               help='Directory for offline calibration frames')

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

    # load database with options to minimize lock time
    db_opts = dict(sync=True, verbose=True, user='spt',
                   host='spt-buffer.grid.uchicago.edu')
    fdb = sdb.TransferFileDatabase(args.file_database, **db_opts)

    # create temporary directory for file transfer
    if args.transfer_root:
        transfer_root = args.transfer_root
    else:
        transfer_root = tempfile.mkdtemp()
    sdb.print_time('Storing temporary transfer files in {}'.format(transfer_root))

    # create date-stamped verified files directory
    verify_root = update_verify_root(args.verify_root)

    # all tarball paths relative to data root
    if not os.path.exists(args.data_root):
        os.makedirs(args.data_root)
    os.chdir(args.data_root)

    # look for files on the remote filesystem
    remote_files = [os.path.basename(x) for x in
                    glob.glob(os.path.join(transfer_root, 'p*.tar'))]

    # process one file at a time
    for f in remote_files:

        entry, priority = sftp.remote2local(f, return_entry=True, return_tag=True)

        if not fdb.match(*entry, status=['received', 'verified']).any():
            continue

        tarfile = os.path.join(transfer_root, f)

        gfiles = []

        try:
            # unpack into data root
            out = sftp.unarchive(tarfile, unpack=True, delete=False,
                                 output=args.data_root, verify=args.verify)

            # check that only g3 files came out
            if not all([x.endswith('.g3') for x in out]):
                for o in out:
                    if os.path.exists(o):
                        sp.check_call(['rm', o])
                raise RuntimeError('Archive {} is corrupt'.format(tarfile))

            gfiles = out

            # remove tarball
            if args.delete:
                sp.check_call(['rm', tarfile])
        except Exception as e:
            # remove failed files from storage
            if gfiles:
                sp.check_call(['rm'] + gfiles)
            sdb.print_warn('Failed to unarchive file {}'.format(f), e)
            # raise e
            fdb.update(*entry, status='error')
            continue
        else:
            if args.verify:
                sdb.print_time('File {} unarchived and verified'.format(f))
                fdb.update(*entry, status='verified')
            else:
                sdb.print_time('File {} unarchived'.format(f))

        # verified file name and entry
        _, src, obs, _ = entry

        # copy verified files to NFS storage
        copy_output = os.path.join(args.copy_root, os.path.dirname(gfiles[0]))
        offline_cal_file = os.path.join(
            args.offline_cal_root, src, '{}.g3'.format(obs))
        offline_cal_link = os.path.join(copy_output, 'offline_calibration.g3')
        sdb.print_time('Copying file {} to {}'.format(gfiles, copy_output))
        try:
            if not os.path.exists(copy_output):
                os.makedirs(copy_output)
            sp.check_call(['rsync', '-aviP'] + gfiles + [copy_output])
            if not os.path.islink(offline_cal_link):
                sp.check_call(['ln', '-s', offline_cal_file, offline_cal_link])
        except Exception as e:
            sdb.print_warn('Error copying {} to {}'.format(gfiles, copy_output), e)
            fdb.update(*entry, status='copy_error')
            continue
        else:
            sdb.print_time('File {} copied to {}'.format(gfiles, copy_output))
            fdb.update(*entry, status='copied')

        # remove buffer copy
        try:
            sp.check_call(['rm'] + gfiles)
        except Exception as e:
            sdb.print_warn('Error removing {} from buffer'.format(gfiles), e)
        else:
            sdb.print_time('File {} removed from buffer'.format(gfiles))

        # create date-stamped verified files directory
        verify_root = update_verify_root(args.verify_root)

        # archive original download
        sdb.print_time('Archiving verified tarfile {} to {}'.format(
            tarfile, verify_root))
        priority = int(fdb.match(*entry, return_df=True, squeeze=True)['priority'])
        tarfile = sftp.local2remote(gfiles[0], priority=priority, aux=False)
        local_tarfile = os.path.join(transfer_root, tarfile)
        verified_tarfile = os.path.join(verify_root, tarfile)
        shutil.move(local_tarfile, verified_tarfile)
        sdb.print_time('Finished processing entry {}'.format(entry))

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
