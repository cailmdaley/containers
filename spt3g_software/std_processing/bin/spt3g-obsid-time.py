#!/usr/bin/env python

import pandas as pd
import argparse as ap
from spt3g.std_processing import obsid_to_g3time, time_to_obsid

P = ap.ArgumentParser(
    description='Get timestamps of input observation ID(s), or the reverse',
    formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('obs_or_time', nargs='+', help='Observation ID(s) or timestamp(s)')
P.add_argument('--tz', default='UTC', help='Input or Output timezone')
P.add_argument('-r', '--reverse', default=False, action='store_true',
               help='Treat the input argument(s) as timestamp(s) and return the corresponding '
               'observation IDs.  If this option is not set, the inputs are treated as '
               'observation IDs, and timestamps are returned.')
P.add_argument('-s', '--simple', default=False, action='store_true',
               help='Simplified output for ease of parsing')

args = P.parse_args()

# check timezone
if args.tz.lower() == 'pole':
    args.tz = 'Pacific/Auckland'
elif args.tz.lower() == 'chicago':
    args.tz = 'US/Central'

# print header
if not args.simple:
    print('Observation\tTime ({})'.format(args.tz))

# process each observation
for arg in args.obs_or_time:
    if args.reverse:
        tm = pd.Timestamp(arg, tz=args.tz)
        obs = time_to_obsid(tm.astimezone('utc').isoformat())
    else:
        obs = int(arg)
        tm = obsid_to_g3time(obs).isoformat()
        tm = pd.Timestamp(tm, tz='utc').astimezone(args.tz)
    if args.simple:
        if args.reverse:
            print(obs)
        else:
            print(tm.isoformat())
    else:
        print('{:11d}\t{}'.format(obs, tm.isoformat()))
