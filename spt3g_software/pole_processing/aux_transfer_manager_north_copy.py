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

P = ap.ArgumentParser(description="Master script for transfer management up north.",
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)

P.add_argument('--database',
               default='/buffer/transfer_database/aux_transfer_north.txt',
               help='Database file for files transferred from the remote')
P.add_argument('--database-south',
               default='/buffer/transfer_database/aux_transfer.txt',
               help='Database file for aux files, rsynced from the south')
P.add_argument('--aux-data-root', default='/buffer',
               help='Default directory for aux output files')
P.add_argument('--tar-data-root', default='/buffer/tar',
               help='Directory for tar output files')
P.add_argument('--pydfmux-data-root', default='/buffer/pydfmux_output',
               help='Directory for pydfmux output files')
P.add_argument('--arc-data-root', default='/buffer/arc',
               help='Directory for arc files')
P.add_argument('--eht-data-root', default='/buffer/eht',
               help='Directory for EHT data files')
P.add_argument('--maser-data-root', default='/buffer/rsync/MaserLogger',
               help='Directory for Maser data files')
P.add_argument('--map-data-root', default='/buffer/onlinemaps',
               help='Directory for mapmaking data files')
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
sdb.print_time('Starting aux transfer manager')
print('*' * 80)

# create SFTP communication object
sftp = SFTPComm(user=args.sftp_user, host=args.sftp_host,
                identity=args.sftp_identity)

transfer_root = args.transfer_root

try:

    # load database with options to minimize lock time
    db_opts = dict(sync=True, verbose=True, user='spt',
                   host='spt-buffer.grid.uchicago.edu')
    db = sdb.AuxTransferDatabase(args.database, **db_opts)
    rdb = sdb.AuxTransferDatabase(args.database_south, **db_opts)

    # create temporary directory for file transfer
    sdb.print_time('Storing temporary transfer files in {}'.format(transfer_root))

    # create date-stamped verified files directory
    verify_root = update_verify_root(args.verify_root)

    # all tarball paths relative to data root
    if not os.path.exists(args.pydfmux_data_root):
        os.makedirs(args.pydfmux_data_root)
    if not os.path.exists(args.arc_data_root):
        os.makedirs(args.arc_data_root)
    if not os.path.exists(args.eht_data_root):
        os.makedirs(args.eht_data_root)
    if not os.path.exists(args.maser_data_root):
        os.makedirs(args.maser_data_root)
    if not os.path.exists(args.tar_data_root):
        os.makedirs(args.tar_data_root)
    if not os.path.exists(args.map_data_root):
        os.makedirs(args.map_data_root)

    # look for files on the remote filesystem
    remote_files = [os.path.basename(x) for x in
                    glob.glob(os.path.join(transfer_root, 'aux*.tar'))]

    # process one file at a time
    for f in remote_files:

        floc, ftype = sftp.remote2local(f, return_tag=True)
        if ftype == 'rsync':
            continue

        if not db.match(floc, status=['received', 'verified']).any():
            continue

        db.update(floc, type=ftype)

        sdb.print_time('Downloading {}'.format(f))

        # set output root
        if ftype in ['pydfmux', 'pydfmuxlog']:
            output = args.pydfmux_data_root
        elif ftype == 'arc':
            output = args.arc_data_root
        elif ftype == 'tar':
            output = args.tar_data_root
        elif ftype == 'eht':
            output = args.eht_data_root
        elif ftype == 'maser':
            output = args.maser_data_root
        elif ftype == 'map':
            output = args.map_data_root
        else:
            sdb.print_warn(
                'Error receiving file {}: unrecognized arhive type {}'.format(
                    f, ftype))
            continue

        tarfile = os.path.join(transfer_root, f)

        # extract and verify archive
        try:
            files = sftp.unarchive(
                tarfile, output=output, delete=False, unpack=False,
                return_checksums=True, verify=True)
        except Exception as e:
            sdb.print_warn('Error extracting {}'.format(tarfile), e)
            db.update(floc, status='error')
            continue
        else:
            db.update(floc, status='verified')

        # copy verified files to NFS storage
        copy_output = output.replace(args.aux_data_root, args.copy_root)

        if ftype == 'map':
            copy_output = os.path.join(copy_output, os.path.dirname(floc))
            if not os.path.exists(copy_output):
                os.makedirs(copy_output)

        # divert bad arc files to a subdirectory
        if ftype == 'arc':
            try:
                gcp.ARCFileReader(str(files[0]))
            except RuntimeError as e:
                if '{} truncated'.format(files[0]) in str(e):
                    sdb.print_warn('Truncated GCP file {}'.format(files[0]))
                    copy_output = os.path.join(copy_output, 'bad')
                else:
                    sdb.print_warn(
                        'Error opening GCP file {}'.format(files[0]), e)
                    db.update(floc, status='copy_error')
                    continue

        sdb.print_time('Copying files {} to {}'.format(files, copy_output))
        try:
            # cmd = ['rsync', '-aviP', '--checksum']
            cmd = ['rsync', '-aviP']
            sp.check_call(cmd + files + [copy_output])
        except Exception as e:
            sdb.print_warn('Error copying {} to {}'.format(files, copy_output), e)
            db.update(floc, status='copy_error')
            continue
        else:
            db.update(floc, status='copied')

        # remove buffer copy
        try:
            sp.check_call(['rm', '-r'] + files)
        except Exception as e:
            sdb.print_warn('Error removing {} from buffer'.format(files), e)
        else:
            sdb.print_time('File {} removed from buffer'.format(files))

        # archive original download
        verify_root = update_verify_root(args.verify_root)
        sdb.print_time('Archiving verified tarfile {} to {}'.format(
            tarfile, verify_root))
        shutil.move(tarfile, os.path.join(verify_root, f))
        sdb.print_time('Entry {} archived successfully'.format(floc))

except Exception as e:

    sdb.print_warn('Aux transfer error', e)
    raise e

finally:

    print('*' * 80)

    # remote transfer lock file
    lock.close()
