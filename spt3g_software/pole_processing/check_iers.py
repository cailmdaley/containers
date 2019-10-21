#!/usr/bin/env python

import warnings
import argparse as ap
import traceback

P = ap.ArgumentParser(description="Check that astropy IERS tables are current")
P.add_argument('-c', '--clear-cache', action='store_true',
               help='Clear download cache and force a database update. Useful if '
               'the downloaded database file has been corrupted.')
P.add_argument('-d', '--auto-max-age', action='store', type=float, default=10.0,
               help='Age in days of the IERS database required to trigger an update')
P.add_argument('-v', '--verbose', action='store_true',
               help='Print status on success or failure. By default, messages '
               'are only printed if the database must be downloaded, or if the '
               'update is unsuccessful.')
args = P.parse_args()

if args.clear_cache:
    if args.verbose:
        print("Clearing astropy download cache")
    from astropy.utils.data import clear_download_cache
    clear_download_cache()

# update the database more often than the default cadence
# this should avoid errors in the scanifier
from astropy.utils import iers
iers.conf.auto_max_age = args.auto_max_age

# computing UT1 will trigger a database update if necessary
from astropy.time import Time
t = Time.now()
try:
    t1 = t.ut1
    if args.verbose:
        print('Current time (UT1):', t.ut1)
        print('IERS database is up to date')
except:
    print(
        "IERS table is out of date.  Run 'check_iers' when the satellite "
        "is up to fix this manually. Traceback information follows.\n\n"
    )
    traceback.print_exc()
