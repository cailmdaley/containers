#!/usr/bin/env python

from spt3g.std_processing import status_db as sdb
import argparse as ap
import inspect
import numpy as np
import socket

parser_opts = dict(formatter_class=ap.ArgumentDefaultsHelpFormatter)
P = ap.ArgumentParser(description="Database updating utility", **parser_opts)

S = P.add_subparsers(dest='db', metavar='DB', title='databases',
                     help='Database type to modify. For more help, call: '
                     '%(prog)s %(metavar)s -h')

dbs = [k for k in vars(sdb).keys() if k.endswith('Database')
       and k != 'StatusDatabase']

# auto-generate database arguments
for db_name in dbs:

    name = db_name.replace('Database', '').lower()
    db_class = getattr(sdb, db_name)

    try:
        spec = inspect.getfullargspec(db_class.__init__)
    except:
        spec = inspect.getargspec(db_class.__init__)

    db_file = spec.defaults[spec.args.index('filename') - len(spec.args)]
    if socket.gethostname().startswith('spt-buffer'):
        db_file = db_file.replace('/data/autoproc/db/', '/buffer/transfer_database/')
        north = True
        if db_name == 'auxtransfer':
            db_file = db_file.replace('aux_transfer.txt', 'aux_transfer_north.txt')
        elif db_name == 'transferfile':
            db_file = db_file.replace('transfer_files.txt', 'transfer_north_files.txt')
    else:
        north = False

    helptext = db_class.__doc__

    PP = S.add_parser(name, help=helptext, **parser_opts)

    PP.add_argument('-f', '--database-filename', default=db_file,
                    metavar='path/to/db.txt', help='Database filename')
    PP.add_argument('-n', '--north', action='store_true',
                    help='North-side user/host')

    for arg in db_class.index_columns:
        key = arg.replace('_', '-')
        PP.add_argument(key, type=db_class.index_dtypes[arg],
                        help=db_class.index_docs[arg])

    for kwarg in db_class.status_columns:
        key = kwarg.replace('_', '-')
        tp = db_class.status_dtypes[kwarg]
        dtype = np.dtype(tp)
        doc = db_class.status_docs[kwarg]
        if dtype.char == '?':
            G = PP.add_mutually_exclusive_group()
            G.add_argument('--{}'.format(key), action='store_true',
                           help=doc)
            G.add_argument('--no-{}'.format(key), action='store_false',
                           dest=kwarg, help=doc)
        else:
            PP.add_argument('--{}'.format(kwarg.replace('_', '-')), default=None,
                            type=tp, help=doc)

    PP.add_argument('-a', '--add', action='store_true',
                    help='Add the requested entry to the database if missing')
    PP.add_argument('-d', '--drop', action='store_true',
                    help='Drop the requested entry from the database')

    PP.set_defaults(db_class=db_class)

args = P.parse_args()

# be explicit about the database user/host
dbargs = dict(sync=True, user='sptdaq', host='anal.spt', match_regex=True)
if north or args.north:
    dbargs.update(user='spt', host='spt-buffer.grid.uchicago.edu')

# open database
print('Opening {} at {}'.format(args.db_class.__name__, args.database_filename))
db = args.db_class(args.database_filename, **dbargs)

# construct update arguments
index = [getattr(args, k) for k in db.index_columns]
status = dict()
for k in db.status_columns:
    v = getattr(args, k, None)
    if v is None:
        continue

    if str(v).lower() == 'none':
        status[k] = None
    else:
        status[k] = v

if args.drop:

    # drop entry from database
    if db.match(*index).any():
        print('Removing entry {}'.format(index))
        db.drop(*index)
    else:
        raise RuntimeError('Missing entry {}'.format(index))

else:

    # update database
    if db.match(*index).any():
        mode = 'Updating'
    elif not args.add:
        msg = 'Missing entry {}. Check your spelling or use --add'
        raise RuntimeError(msg.format(index))
    else:
        mode = 'Adding'
    if mode:
        print('{} entry {} with status {}'.format(mode, index, status))
        db.update(*index, **status)

# close
db.close()
