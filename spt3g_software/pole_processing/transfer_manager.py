#!/usr/bin/env python

from spt3g.std_processing import status_db as sdb
import argparse as ap
import numpy as np
import os
import subprocess as sp
import shutil
import glob
import tempfile
import sys
from spt3g.std_processing.transfer_tools import SFTPComm, RemoteFullError, SFTPError
from spt3g.std_processing.transfer_tools import LockedProcess, LockedFileError

P = ap.ArgumentParser(description="Master script for transfer management.",
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)

P.add_argument('--database', default='/data/autoproc/db/scanify.txt',
               help='Transfer management database file')
P.add_argument('--file-database',
               default='/data/autoproc/db/transfer_files.txt',
               help='Database file for files transfered to the remote')
P.add_argument('--data-root', default='/spt_data/bolodata',
               help='Directory for scanified observations')
P.add_argument('--max-transfer-size', default=125, type=int,
               help='Maximum disk usage on sftp server, in GB')
P.add_argument('--sftp-user', action='store', default='A-379-S',
               help='SFTP server username')
P.add_argument('--sftp-host', action='store', default='spfs.southpole.usap.gov',
               help='SFTP server hostname')
P.add_argument('--sftp-identity', action='store', default='~/.ssh/id_dsa_spfs',
               help='SFTP secure identity file')
P.add_argument('--maxdata', default=20, type=int,
               help='Default amount of data in GB to upload at once to the server.')

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
db = None
fdb = None
adb = None

# bundle nominal cal file with first obs file
rootfile = '0000.g3'
calfile = 'nominal_online_cal.g3'
isroot = lambda x: x.endswith(rootfile)
tocal = lambda x: x.replace(rootfile, calfile)

def get_obs_files(dirname):
    files = glob.glob(os.path.join(dirname, '*.g3'))
    files = filter(lambda x: (not x.endswith('offline_calibration.g3') and
                              not x.endswith(calfile)),
                   files)
    return sorted(files)

try:

    # load database to with options to minimize lock time
    db_opts = dict(sync=True, verbose=True)
    db = sdb.ScanifyDatabase(args.database, **db_opts)
    fdb = sdb.TransferFileDatabase(args.file_database, **db_opts)

    # look for files on the remote filesystem
    remote_files, remote_sizes = sftp.list_files('p*.tar', return_size=True)
    remote_entries = [sftp.remote2local(fname, return_entry=True)
                      for fname in remote_files]
    sdb.print_time('Found {} priority files on remote'.format(len(remote_files)))

    # look for bundles on the remote filesystem
    remote_bundles = sftp.list_files('bundle_p*.tar', return_size=False)
    remote_bundle_entries = sum([fdb.get_entries(archive=x)
                                 for x in remote_bundles], [])
    sdb.print_time('Found {} priority bundles ({} files) on remote'.format(
        len(remote_bundles), len(remote_bundle_entries)))

    # update file transfer database
    uploaded_files = fdb.get_entries(status='uploaded', archive=None)
    uploaded_bundles = fdb.match(status='uploaded', return_df=True)['archive'] \
                          .replace('', np.nan).dropna().unique()
    with fdb.commit_all():
        # mark all stuff on the server as uploaded
        for entry, sz in zip(remote_entries, remote_sizes):
            fdb.update(*entry, status='uploaded', size=sz)
        for entry in remote_bundle_entries:
            fdb.update(*entry, status='uploaded')
        # mark all stuff missing from the server as sent
        for entry in uploaded_files:
            if entry not in remote_entries:
                fdb.update(*entry, status='sent')
        for bundle in uploaded_bundles:
            if bundle not in remote_bundles:
                for entry in fdb.get_entries(archive=bundle):
                    fdb.update(*entry, status='sent')

    max_priority = max(db.data['transfer_downsampled'].max(),
                       db.data['transfer_fullrate'].max())
    priorities = list(range(1, max_priority + 1))

    # update obs transfer database
    with db.commit_all():
        for rate in ['downsampled', 'fullrate']:
            status_key = 'status_{}'.format(rate)
            transfer_key = 'transfer_{}'.format(rate)
            uploading_obs = db.get_entries(**{status_key: 'uploading'})
            for src, obs in uploading_obs:
                files = fdb.match(rate, src, obs, return_df=True)
                if not len(files):
                    continue
                # allow errors here, obs will be reprocessed after verification
                if files['status'].isin(['uploaded', 'sent', 'verified', 'error']).all():
                    db.update(src, obs, **{status_key: 'uploaded'})
            uploaded_obs = db.get_entries(**{status_key: 'uploaded'})
            for src, obs in uploaded_obs:
                files = fdb.match(rate, src, obs, return_df=True)
                if not len(files):
                    continue
                # allow errors here, obs will be reprocessed after verification
                if files['status'].isin(['sent', 'verified', 'error']).all():
                    db.update(src, obs, **{status_key: 'sent'})
            sent_obs = db.get_entries(**{status_key: 'sent'})
            for src, obs in sent_obs:
                files = fdb.match(rate, src, obs, return_df=True)
                if not len(files):
                    continue
                if files['status'].isin(['uploaded']).any():
                    db.update(src, obs, **{status_key: 'uploaded'})
                if files['status'].isin(['uploading']).any():
                    db.update(src, obs, **{status_key: 'uploading'})
                if files['status'].isnull().any():
                    db.update(src, obs, **{status_key: 'ready'})

            # find all waiting observations
            ready_obs = db.get_entries(**{status_key: 'ready',
                                          transfer_key: priorities})
            with fdb.commit_all():
                for src, obs in ready_obs:
                    priority = int(db.match(src, obs, return_df=True, squeeze=True)[transfer_key])
                    if fdb.match(rate, src, obs, priority=priority).any():
                        continue
                    obs_files = get_obs_files(os.path.join(
                        args.data_root, rate, src, str(obs)))
                    if not len(obs_files):
                        continue
                    for f in obs_files:
                        sz = sftp.check_file_size_local(f)
                        if isroot(f):
                            sz += sftp.check_file_size_local(tocal(f))
                        fdb.update(
                            rate, src, obs, os.path.basename(f), size=sz, priority=priority)

            # remove observations not marked for transfer
            pull_obs = db.get_entries(**{status_key: 'ready', transfer_key: 0})
            with fdb.commit_all():
                for src, obs in pull_obs:
                    if fdb.match(rate, src, obs).any():
                        fdb.drop(rate, src, obs)

    # check remote file size
    available_size = sftp.remote_available(args.max_transfer_size * 1024**3)
    sdb.print_time('{} GB available on remote'.format(available_size / 1024.**3))

    # make temporary directory for file transfer
    transfer_root = tempfile.mkdtemp()
    sdb.print_time('Storing temporary transfer files in {}'.format(transfer_root))

    # all tarball paths relative to data root
    os.chdir(args.data_root)

    total = 0

    for priority in range(1, max_priority + 1):

        for rate in ['downsampled', 'fullrate']:

            status_key = 'status_{}'.format(rate)
            transfer_key = 'transfer_{}'.format(rate)

            # count observations
            obs_list = db.get_entries(**{status_key: ['ready', 'uploading'],
                                         transfer_key: priority})
            msg = 'Found {} {} observations ready to transfer with priority {}'
            sdb.print_time(msg.format(len(obs_list), rate, priority))

            for src, obs in obs_list:
                if total >= args.maxdata:
                    break

                obs_root = os.path.join(rate, src, str(obs))
                obs_root_abs = os.path.abspath(obs_root)
                obs_files = get_obs_files(obs_root)

                # make sure all entries that are to be processed are recorded
                with fdb.commit_all():
                    for obs_file in obs_files:
                        obs_base = os.path.basename(obs_file)
                        fdb.update(rate, src, obs, obs_base)

                try:
                    obs_entries = []

                    for obs_file in obs_files:
                        if total >= args.maxdata:
                            sdb.print_time('Uploaded {} GB, stopping for now'.format(total))
                            break

                        obs_file_abs = os.path.abspath(obs_file)
                        remote = sftp.local2remote(obs_file, priority=priority)
                        obs_base = os.path.basename(obs_file)
                        tarfile = os.path.join(transfer_root, remote)

                        # skip individual files that have already been sent
                        if fdb.match(rate, src, obs, obs_base,
                                     status=['uploading', 'uploaded', 'sent', 'verified']).any():
                            continue

                        # remove preloaded but not sent nominal cal files
                        if isroot(obs_file):
                            fdb.drop(rate, src, obs, calfile, status=None)

                        available_size = sftp.remote_available(
                            args.max_transfer_size * 1024**3)

                        status = None
                        bundle = None

                        try:
                            # check file size
                            sz = sftp.check_file_size_local(obs_file)
                            if isroot(obs_file):
                                sz += sftp.check_file_size_local(tocal(obs_file))
                            if sz > available_size:
                                msg = ('Not enough space remaining on remote. '
                                       'Found {} MB to transfer file {} of size {} MB')
                                raise RemoteFullError(msg.format(
                                    available_size / 1024**2, obs_file_abs,
                                    sz / 1024**2))

                            if isroot(obs_file):
                                obs_file = [obs_file, tocal(obs_file)]

                            # create tarball
                            os.chdir(args.data_root)
                            tarfile = sftp.archive(
                                obs_file, compress=True, relative=True, output=tarfile)
                            sz = sftp.check_file_size_local(tarfile)

                            # upload
                            sdb.print_time('Sending {}'.format(tarfile))
                            sz, status, bundle = sftp.upload_file(tarfile, remote)

                            total += sz / float(1024**3)

                            sp.check_call(['rm', '-f', tarfile])
                        except Exception as e:
                            # remove from file database on failure
                            if not isinstance(e, RemoteFullError):
                                sdb.print_warn(
                                    'Error uploading {} to remote'.format(tarfile), e)
                            raise e
                        else:
                            # update file database on success
                            msg = 'File {} uploaded successfully'
                            sdb.print_time(msg.format(obs_file_abs))
                            obs_entries.append(((rate, src, obs, obs_base), sz, status, bundle))
                        finally:
                            sp.check_call('rm -f {}/*'.format(transfer_root), shell=True)

                except Exception as e:
                    # don't process any more if no disk space remaining
                    if isinstance(e, (RemoteFullError, SFTPError)):
                        raise e
                    sdb.print_warn(
                        'Failed to transfer observation {}'.format(obs_root), e)
                else:
                    # transfer complete
                    sdb.print_time('Transfer complete for observation {}'.format(obs_root))
                    db.update(src, obs, **{status_key: 'uploading'})
                finally:
                    with fdb.commit_all():
                        for obs_entry, sz, status, bundle in obs_entries:
                            fdb.update(*obs_entry, size=sz, status=status, archive=bundle)
                            if bundle is None or bundle == 'current':
                                continue
                            bundle_entries = fdb.get_entries(status='uploading', archive='current')
                            for bentry in bundle_entries:
                                fdb.update(*bentry, status='uploaded', archive=bundle)

except Exception as e:

    if isinstance(e, RemoteFullError):
        sdb.print_time(e)
    elif isinstance(e, SFTPError) and ("Broken pipe" in str(e) or "Timed out" in str(e)):
        sdb.print_time(e)
    elif isinstance(e, IOError) and ("Resource deadlock avoided" in str(e)):
        sdb.print_time(e)
    else:
        sdb.print_warn('Error transfering data', e)

finally:

    # cleanup transfer database
    if db is not None:
        db.write_sql()

    # cleanup tmp
    if transfer_root:
        shutil.rmtree(transfer_root)

    print('*' * 80)

    # remove transfer lock file
    lock.close()
