#!/usr/bin/env python

import argparse as ap
import os
from spt3g.std_processing import status_db as sdb
from spt3g.std_processing.transfer_tools import LockedProcess, LockedFileError
from spt3g.std_processing.transfer_tools import update_verify_root

P = ap.ArgumentParser(description="Archive transfered data files for long-term storage.",
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('--verify-root', default='/buffer/3g_verified',
               help='Verified transfer data root')

args = P.parse_args()

# create transfer lock file
try:
    lock = LockedProcess(os.path.abspath(__file__), timeout=5)
except LockedFileError as e:
    raise SystemExit

print('*' * 80)
sdb.print_time('Starting archive manager')
print('*' * 80)

try:

    # create date-stamped verified files directory
    verify_root = update_verify_root(args.verify_root, archive=True)

except Exception as e:

    sdb.print_warn('Archive error', e)
    raise e

finally:

    print('*' * 80)

    # remote transfer lock file
    lock.close()
