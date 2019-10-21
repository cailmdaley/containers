#!/usr/bin/env python

from spt3g.std_processing import status_db as sdb
import argparse as ap
import os
import shutil
import subprocess as sp
import tempfile
from spt3g.std_processing.transfer_tools import RsyncComm, RsyncConnectionError
from spt3g.std_processing.transfer_tools import isvalidname, InvalidFilenameError
from spt3g.std_processing.transfer_tools import LockedProcess, LockedFileError

P = ap.ArgumentParser(description="Master script for rsync management.",
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)

P.add_argument('--database', default='/data/autoproc/db/aux_transfer.txt',
               help='Transfer management database file')
P.add_argument('--rsync-user', default='spt',
               help='Username on north machine')
P.add_argument('--rsync-host', default='spt-buffer.grid.uchicago.edu',
               help='Hostname of north machine')
P.add_argument('--rsync-identity', default='~/.ssh/id_rsa_uc',
               help='SSH secure identity file')
P.add_argument('--rsync-data-root', default='/buffer/rsync',
               help='Remote data root for storing rsynced data files')

args = P.parse_args()

# aux transfer lock file
try:
    lock = LockedProcess(os.path.abspath(__file__), timeout=5)
except LockedFileError as e:
    raise SystemExit

print('*' * 80)
sdb.print_time('Starting rsync manager')
print('*' * 80)

transfer_root = None

try:

    rsync = RsyncComm(user=args.rsync_user, host=args.rsync_host,
                      identity=args.rsync_identity)
    if not rsync.isup():
        raise RsyncConnectionError('No rsync connection available')

    # load database
    db = sdb.AuxTransferDatabase(args.database, sync=True)

    # make temporary directory for file transfer
    transfer_root = tempfile.mkdtemp()
    sdb.print_time('Storing temporary transfer files in {}'.format(transfer_root))

    for fname in db.get_entries(type='rsync'):

        if not isvalidname(fname):
            raise InvalidFilenameError(
                'Invalid {} filename {}, removing from database'.format(ftype, fname))

        db.update(fname, status='uploading')

        rsync_touch = os.path.splitext(
            os.path.join(os.path.dirname(args.database), rsync.local2remote(
                fname, aux=True, auxtype='rsync')))[0] + '.last_rsync'
        sent_bytes = rsync.send_file(fname, args.rsync_data_root, update=True,
                                     copy_dir_links=True)

        # rsync last-modified touch-file if transfer is complete
        if sent_bytes:
            sp.check_call(['touch', rsync_touch])
            touch_bytes = rsync.send_file(rsync_touch, args.rsync_data_root)
            if touch_bytes:
                db.update(fname, status='uploaded', size=sent_bytes)

except Exception as e:

    if isinstance(e, RsyncConnectionError):
        sdb.print_time(e)
    else:
        sdb.print_warn('Rsync transfer error', e)

finally:

    # cleanup tmp
    if transfer_root:
        shutil.rmtree(transfer_root)

    print('*' * 80)

    # remove aux transfer lock file
    lock.close()
