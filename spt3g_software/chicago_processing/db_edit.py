#!/usr/bin/env python

from __future__ import print_function
from spt3g.std_processing import status_db
import sys, os, pandas, fnmatch
import argparse as ap
import socket

# Usage set_status.py [-d /path/to/autoprocdb.txt] new_status source id

host = socket.gethostname()
if 'anal' in host:
    database = '/data/autoproc/db/autoproc.txt'
else:
    database = '/spt/user/production/autoproc.txt'

P = ap.ArgumentParser(description='Edit autoprocessing database by hand',
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('status', action='store',
               help='Action to take. The specified observation will be set '
                    'to the status specified. Reasonable choices are: '
                    '"complete", "rerun", "ignore".')
P.add_argument('source', action='store', help='Source to modify.')
P.add_argument('obs', action='store', nargs='?', type=str, default=None,
               help='Observation ID to modify. If unspecified, changes the '
               'status of all observations of the specified source.')
P.add_argument('-d', '--database', action='store', default=database,
               help='Path to database to modify')
P.add_argument('-p', '--previous-status', action='store', default=None,
               help='Only modify observations with a given previous status')
P.add_argument('-y', '--yes', action='store_true',
               default=False, help='Assume answers to all questions are yes')
args = P.parse_args()

if args.status == 'rerun':
    args.status = None
if args.source == 'all':
    args.source = '.*'

host = socket.gethostname()
if 'anal' in host:
    autoproc_db = status_db.AutoprocDatabase(args.database, user='sptdaq',
                                             host=host, match_regex=True)
else:
    autoproc_db = status_db.AutoprocDatabase(args.database, user='nwhitehorn',
                                             host='scott.grid.uchicago.edu',
                                             match_regex=True)

# Use python-3 style 'input' everywhere
try:
    input = raw_input
except NameError:
    pass

extra_match_args={}
if args.previous_status is not None:
    extra_match_args['status'] = args.previous_status

if args.obs is None or args.obs[0] == '>' or args.obs[0] == '<' or (args.obs[0] != '-' and '-' in args.obs):
    entries = autoproc_db.match(args.source, return_df=True, **extra_match_args)
    if args.obs is not None and args.obs[0] == '<':
        obscutoff = int(args.obs[1:])
        print('Setting the status of all %s jobs before %d to %s' %
          (args.source, obscutoff, args.status))
        entries = [x for x in zip(entries['source'], entries['observation'])
          if x[1] < obscutoff]
    elif args.obs is not None and args.obs[0] == '>':
        obscutoff = int(args.obs[1:])
        print('Setting the status of all %s jobs after %d to %s' %
          (args.source, obscutoff, args.status))
        entries = [x for x in zip(entries['source'], entries['observation'])
          if x[1] > obscutoff]
    elif args.obs is not None and args.obs[0] != '-' and '-' in args.obs:
        obsstart = int(args.obs.split('-')[0])
        obsstop = int(args.obs.split('-')[1])
        print('Setting the status of all %s jobs from %d to %d to %s' %
          (args.source, obsstart, obsstop, args.status))
        entries = [x for x in zip(entries['source'], entries['observation'])
          if x[1] >= obsstart and x[1] <= obsstop]
    else:
        print('Setting the status of all %s jobs to %s' %
          (args.source, args.status))
        entries = [x for x in zip(entries['source'], entries['observation'])]
    if len(entries) == 0:
        print('No matching entries, exiting')
        sys.exit(1)
    if args.yes:
        print('%d entries match.' % len(entries))
    else:
        answer = input('%d entries match. Proceed? [y/n]' % len(entries))
        if answer != 'y':
            sys.exit(1)
    for src, obs in entries:
        autoproc_db.update(src, obs, status=args.status, match_regex=False)
else:
    obs = int(args.obs)
    print('Setting the status of %s observation %d to %s' %
      (args.source, obs, args.status))
    entries = autoproc_db.match(args.source, obs, return_df=True, match_regex=False, **extra_match_args)
    if len(entries) == 0 and not args.yes:
        answer = input('No matching entries, create observation %s/%d? [y/n]' %
          (args.source, obs))
        if answer != 'y':
            sys.exit(1)
    autoproc_db.update(args.source, obs, status=args.status, match_regex=False)

autoproc_db.commit()
autoproc_db.close()

