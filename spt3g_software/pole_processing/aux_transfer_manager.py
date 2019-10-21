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
import warnings
import time
import datetime
import pytz
import pandas as pd
from spt3g.std_processing.transfer_tools import SFTPComm, RemoteFullError, SFTPError
from spt3g.std_processing.transfer_tools import isvalidname, InvalidFilenameError
from spt3g.std_processing.transfer_tools import LockedProcess, LockedFileError

P = ap.ArgumentParser(description="Master script for transfer management.",
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)

P.add_argument('--database', default='/data/autoproc/db/aux_transfer.txt',
               help='Transfer management database file')
P.add_argument('--pydfmux-data-root', default='/poleanalysis/pydfmux_output',
               help='Directory for pydfmux output files')
P.add_argument('--arc-data-root', default='/spt_data/arc',
               help='Directory for raw arcfile output')
P.add_argument('--eht-data-root', default='/buffer/eht',
               help='Directory for EHT data for upload')
P.add_argument('--maser-data-root', default='/data/MaserLogger/backup',
               help='Directory for MaserLogger output')
P.add_argument('--map-data-root', default='/poleanalysis/sptdaq/calresult/onlinemaps',
               help='Directory for field maps')
P.add_argument('--max-transfer-size', default=125, type=int,
               help='Maximum disk usage on sftp server, in GB')
P.add_argument('--max-file-size', default=1024, type=int,
               help='Maximum size of individual files to transfer, MB. '
               'Files larger than this will be split into segments')
P.add_argument('--sftp-user', action='store', default='A-379-S',
               help='SFTP server username')
P.add_argument('--sftp-host', action='store', default='spfs.southpole.usap.gov',
               help='SFTP server hostname')
P.add_argument('--sftp-identity', action='store', default='~/.ssh/id_dsa_spfs',
               help='SFTP secure identity file')
P.add_argument('--no-pydfmux', dest='pydfmux', action='store_false', default=True,
               help='Upload pydfmux outputs')
P.add_argument('--no-arc', dest='arc', action='store_false', default=True,
               help='Upload arc outputs')
P.add_argument('--no-eht', dest='eht', action='store_false', default=True,
               help='Upload EHT data files')
P.add_argument('--no-maser', dest='maser', action='store_false', default=True,
               help='Upload Maser data backup files')
P.add_argument('--no-maps', dest='maps', action='store_false', default=True,
               help='Upload field maps')
P.add_argument('--pydfmux-timezone', default='utc',
               help='Pydfmux data output timezone')

args = P.parse_args()

# aux transfer lock file
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

transfer_root = None

# timezone-aware date stamps
now = time.time()
tz = pytz.timezone(args.pydfmux_timezone)
today_date = datetime.datetime.fromtimestamp(now, tz).date()
one_day = datetime.timedelta(days=1)
today = today_date.strftime('%Y%m%d')
yesterday = (today_date - one_day).strftime('%Y%m%d')
tomorrow = (today_date + one_day).strftime('%Y%m%d')

# return last-modified time of the pydfmux data directory.
# assume it contains a directory of log files
# and return the most recent timestamp
def get_pydfmux_mtime(path):
    if os.path.exists(os.path.join(path, 'logs')):
        logfiles = glob.glob(os.path.join(path, 'logs', '*.txt'))
        allfiles = glob.glob(os.path.join(path, '*'))
        return max([os.stat(x).st_mtime for x in logfiles + allfiles + [path]])
    logfiles = glob.glob(os.path.join(path, '*', 'logs', '*.txt'))
    allfiles = glob.glob(os.path.join(path, '*', '*'))
    return max([os.stat(x).st_mtime for x in logfiles + allfiles + [path]])

db = None

try:

    # load database
    db = sdb.AuxTransferDatabase(args.database, sync=True)

    # look for files on the remote file system
    remote_files, remote_sizes = sftp.list_files('aux*.tar*', return_size=True)
    remote_entries = [sftp.remote2local(fname, return_entry=False)
                      for fname in remote_files]
    sdb.print_time('Found {} aux files on remote'.format(len(remote_files)))

    # look for bundles on the remote file system
    remote_bundles = sftp.list_files('bundle_aux*.tar*', return_size=False)
    remote_bundle_entries = sum([list(db.get_entries(archive=x))
                                 for x in remote_bundles], [])
    sdb.print_time('Found {} aux bundles ({} files) on remote'.format(
        len(remote_bundles), len(remote_bundle_entries)))

    # not uploaded files
    not_uploaded_files = db.get_entries(status=[None, 'uploading'])
    with db.commit_all():
        for entry in not_uploaded_files:
            if entry in remote_entries or entry in remote_bundle_entries:
                db.update(entry, status='uploaded')

    # update file transfer database
    uploaded_files = db.get_entries(status=['uploaded', 'sent'])
    with db.commit_all():
        for entry in uploaded_files:
            if entry not in remote_entries and entry not in remote_bundle_entries:
                db.update(entry, status='sent')
            elif entry in remote_entries:
                idx = remote_entries.index(entry)
                f = remote_files[idx]
                sz = remote_sizes[idx]
                db.update(entry, status='uploaded', size=sz)
            else:
                db.update(entry, status='uploaded')

    if args.pydfmux:

        # pydfmux data
        sdb.print_time('Checking pydfmux directories')

        # find all subdirectories within date directories
        pydfmux_dirs = sorted(glob.glob(os.path.join(
                    args.pydfmux_data_root, '[0-9]'*8, '*')))
        # find all processed entries
        pydfmux_done = db.get_entries(type='pydfmux')

        # restrict to unprocessed subdirectories with valid dates
        def check_pydfmux_dir(f):
            if not os.path.isdir(f):
                return False
            dt = os.path.basename(os.path.dirname(f))
            if dt in pydfmux_done:
                return False
            if dt > tomorrow:
                return False
            if dt < '20170101':
                return False
            relf = os.path.relpath(f, args.pydfmux_data_root)
            if relf in pydfmux_done:
                return False
            if dt < today:
                return os.path.join(args.pydfmux_data_root, dt)
            return False

        pydfmux_dirs = [check_pydfmux_dir(x) for x in pydfmux_dirs]
        pydfmux_dirs = [x for x in pydfmux_dirs if x]
        pydfmux_dirs = list(np.unique(pydfmux_dirs))

        with db.commit_all():
            for f in pydfmux_dirs:
                relf = os.path.relpath(f, args.pydfmux_data_root)
                if not os.path.exists(f):
                    # skip missing directory
                    # (may have been moved while this script was executing)
                    continue
                elif not len(os.listdir(f)):
                    # skip if empty
                    continue
                else:
                    # queue for upload
                    # if the last update was over six hours ago
                    if (now - get_pydfmux_mtime(f)) > 6 * 3600:
                        sz = sftp.check_disk_usage_local(f)
                        db.update(relf, type='pydfmux', size=sz)

        # pydfmux logs
        pydfmux_logs = sorted(glob.glob(os.path.join(
                    args.pydfmux_data_root, '*.log*')))
        for f in pydfmux_logs:
            if f.endswith('.current'):
                continue # don't double-count current log

            relf = os.path.relpath(f, args.pydfmux_data_root)

            # make static copy of active log files every hour
            if f.endswith('.log'):
                fcur = f + '.current'

                mtime = int(os.stat(f).st_mtime)
                now = int(time.time())
                if (now - mtime) > 86400:
                    # skip stale log files
                    tdelta = 0
                    tdelta_current = 0
                else:
                    if not os.path.exists(fcur):
                        # static copy doesn't exist
                        tdelta = np.inf
                        tdelta_current = np.inf
                    else:
                        # compare timestamp of active log and its static copy
                        mtime_current = int(os.stat(fcur).st_mtime)
                        tdelta = mtime - mtime_current
                        tdelta_current = now - mtime_current

                # wait at least an hour to copy
                if (tdelta > 0) and (tdelta_current > 3600):
                    try:
                        sdb.print_time('Creating static copy of {}'.format(f))
                        sp.check_call(['rsync', '-aviP', f, fcur])
                    except Exception as e:
                        sdb.print_warn(
                            'Error making static copy of {}'.format(f), e)
                        continue # skip if there are problems copying
                    else:
                        ready = True
                else:
                    ready = False

                if ready:
                    relf += '.current'

            df = db.match(relf, return_df=True)
            sz = sftp.check_file_size_local(f)
            if not len(df) and (ready or not relf.endswith('.log')):
                # add to database if missing and not a new log file
                db.update(relf, type='pydfmuxlog', size=sz)
            elif f.endswith('.log') and df['status'].isin(['sent', 'verified']).all() and ready:
                # resend recently updated log files
                db.update(relf, status=None, size=sz)

    if args.arc:

        sdb.print_time('Checking arc files')

        # arcfile data
        arc_data = [x.replace('.verified', '') for x in
                    sorted(glob.glob(os.path.join(
                        args.arc_data_root, '*.verified')))]
        arc_entries = db.get_entries(type='arc')
        now = int(time.time())
        with db.commit_all():
            for f in arc_data:
                relf = os.path.relpath(f, args.arc_data_root)
                if relf not in arc_entries:
                    # skip bad files that have already been moved out of the way
                    if not os.path.exists(f):
                        continue
                    mtime = int(os.stat(f).st_mtime)
                    # wait at least an hour to add to queue, in case
                    # the scanifier finds that the file is bad
                    if ((now - mtime) < 3600):
                        continue
                    sz = sftp.check_file_size_local(f)
                    db.update(relf, type='arc', size=sz)

    if args.eht:

        sdb.print_time('Checking EHT files')

        # EHT data
        eht_data = sorted(glob.glob(os.path.join(
            args.eht_data_root, '*.vdif')))
        eht_entries = db.get_entries(type='eht')
        with db.commit_all():
            for f in eht_data:
                relf = os.path.relpath(f, args.eht_data_root)
                if relf not in eht_entries:
                    sz = sftp.check_file_size_local(f)
                    db.update(relf, type='eht', size=sz)

    if args.maser:

        sdb.print_time('Checking Maser files')

        # MaserLogger data
        maser_data = sorted(glob.glob(os.path.join(
            args.maser_data_root, '*.tar.gz')))
        maser_entries = db.get_entries(type='maser')
        with db.commit_all():
            for f in maser_data:
                if os.stat(f).st_mtime > time.time() - 600:
                    continue
                relf = os.path.relpath(f, args.maser_data_root)
                if relf not in maser_entries:
                    sz = sftp.check_file_size_local(f)
                    db.update(relf, type='maser', size=sz)

    if args.maps:

        sdb.print_time('Checking map files')

        with db.commit_all():
            df = db.match(type='map', status=None, return_df=True)
            for idx, m in df.iterrows():
                relf = m['filename']
                sz = m['size']
                if not sz:
                    # entries with no timestamp were just added by the mapmaking manager
                    f = os.path.join(args.map_data_root, relf)
                    if not os.path.exists(f):
                        db.drop(relf)
                        continue
                    if os.stat(f).st_mtime > time.time() - 600:
                        continue
                    sz = sftp.check_file_size_local(f)
                else:
                    # reset timestamps for older data to push them to the bottom of the queue
                    obs = int(os.path.basename(relf).split('_')[0])
                    if obs < 74803353:
                        db.update(relf, size=0)
                    else:
                        # leave newer data as is
                        continue
                db.update(relf, size=sz)

    db.write()

    # find files to upload
    df_up = db.match(status=None, return_df=True)
    sdb.print_time('Found {} aux files to upload'.format(len(df_up)))

    if not len(df_up):
        raise SystemExit

    # check remote file size
    available_size = sftp.remote_available(args.max_transfer_size * 1024**3)
    sdb.print_time('{} GB available on remote'.format(available_size / 1024.**3))

    # make temporary directory for file transfer
    transfer_root = tempfile.mkdtemp()
    sdb.print_time('Storing temporary transfer files in {}'.format(transfer_root))

    for idx, df in df_up.iterrows():

        fname = df['filename']
        ftype = df['type']
        fmod = df['modified']

        status = None
        bundle = None

        try:

            if not isvalidname(fname):
                raise InvalidFilenameError(
                    'Invalid {} filename {}, removing from database'.format(ftype, fname))

            # change root directory
            if ftype in ['pydfmux', 'pydfmuxlog']:
                if not args.pydfmux:
                    continue
                os.chdir(args.pydfmux_data_root)
            elif ftype == 'arc':
                if not args.arc:
                    continue
                os.chdir(args.arc_data_root)
            elif ftype == 'eht':
                if not args.eht:
                    continue
                os.chdir(args.eht_data_root)
            elif ftype == 'maser':
                if not args.maser:
                    continue
                os.chdir(args.maser_data_root)
            elif ftype == 'map':
                if not args.maps:
                    continue
                os.chdir(args.map_data_root)
                if os.stat(fname).st_mtime > time.time() - 600:
                    continue
            elif ftype == 'rsync':
                continue
            else:
                os.chdir(os.path.dirname(fname))

            if not os.path.exists(fname):
                raise OSError('Missing aux file {}'.format(fname))

            # create tarball
            tarfile = sftp.archive(
                fname, output=transfer_root, aux=True, auxtype=ftype,
                relative=True, checksum=ftype not in ['map'],
                compress=ftype not in ['arc', 'eht', 'maser', 'map'],
                verbose=ftype not in ['pydfmux'])

            # upload
            available_size = sftp.remote_available(args.max_transfer_size * 1024**3)
            sz = sftp.check_file_size_local(tarfile)
            if sz > available_size:
                msg = ('Not enough space remaining on remote. '
                       'Found {} MB for file {} of size {} MB')
                raise RemoteFullError(msg.format(
                    available_size / 1024**2, tarfile, sz / 1024**2))
            sz, status, bundle = sftp.upload_file(tarfile, os.path.basename(tarfile),
                                                  aux=True, auxtype=ftype)

            # cleanup
            sp.check_call(['rm', '-f', tarfile])

        except Exception as e:

            if isinstance(e, InvalidFilenameError):
                db.drop(fname)

            if isinstance(e, SFTPError):
                # This could just be a file that takes too long to upload
                sdb.print_time('Error uploading file {}: {}'.format(fname, e))
                continue

            if isinstance(e, RemoteFullError):
                raise e

            sdb.print_warn('Error uploading file {}'.format(fname), e)

        else:

            if len(sftp.list_files(tarfile)) or (status == 'uploading') or \
               (bundle and len(sftp.list_files(bundle))):
                sdb.print_time('Successfully uploaded aux file {}'.format(fname))
                with db.commit_all():
                    db.update(fname, size=sz, status=status, archive=bundle)
                    if bundle is None or bundle == 'current':
                        continue
                    bundle_entries = db.get_entries(status='uploading', archive='current')
                    for bentry in bundle_entries:
                        db.update(bentry, status='uploaded', archive=bundle)
            else:
                sdb.print_time('Aux file {} upload preempted'.format(fname))

        finally:
            sp.check_call('rm -f {}/*'.format(transfer_root), shell=True)

except Exception as e:

    if isinstance(e, (SystemExit, RemoteFullError)):
        sdb.print_time(e)
    elif isinstance(e, SFTPError) and ("Broken pipe" in str(e) or "Timed out" in str(e)):
        sdb.print_time(e)
    else:
        sdb.print_warn('Aux transfer error', e)

finally:

    if db is not None:
        db.write_sql()

    # cleanup tmp
    if transfer_root:
        shutil.rmtree(transfer_root)

    print('*' * 80)

    # remove aux transfer lock file
    lock.close()
