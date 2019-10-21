#!/usr/bin/env python

import argparse as ap
from spt3g.std_processing.transfer_tools import SFTPComm
import os

P = ap.ArgumentParser(description="Master script for transfer management.",
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)

P.add_argument('--north', default=None, action='store_true',
               help='Use north-side SFTP parameters')
P.add_argument('--user', action='store', help='SFTP server username')
P.add_argument('--host', action='store', help='SFTP server hostname')
P.add_argument('--identity', action='store', help='SFTP secure identity file')
P.add_argument('cmd', default="ls -ltr", nargs='?', help='SFTP command')
G = P.add_mutually_exclusive_group()
G.add_argument('-u', '--usage', action='store_true',
               help='Print the amount of data in the queue in GB')
G.add_argument('-l', '--login', action='store_true',
               help='Login to the remote server')

args = P.parse_args()

if args.north is None:
    args.north = os.getenv('HOSTNAME') == 'spt-buffer.grid.uchicago.edu'

if args.north:

    if args.user is None:
        args.user = 'A-379-S'
    if args.host is None:
        args.host = 'sftp.usap.gov'
    if args.identity is None:
        args.identity = '~/.ssh/A-379-S-2'

else:

    if args.user is None:
        args.user = 'A-379-S'
    if args.host is None:
        args.host = 'spfs.southpole.usap.gov'
    if args.identity is None:
        args.identity = '~/.ssh/id_dsa_spfs'

sftp = SFTPComm(user=args.user, host=args.host, identity=args.identity)
if args.usage:
    print('Found {:.3f} GB of data'.format(sftp.remote_used() / 1024. ** 3))
elif args.login:
    sftp.login()
else:
    print(sftp(args.cmd, parse=False).strip())
    print('Total: {:.3f} GB'.format(sftp.remote_used() / 1024. ** 3))
