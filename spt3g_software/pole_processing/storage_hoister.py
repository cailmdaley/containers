#!/usr/bin/env python

from spt3g.std_processing import status_db
import sys, os, shutil, glob
import argparse as ap
import subprocess as sp
import warnings

# Copy files associated with observations marked 'buffered' in the transfer
# DB to two sets of disks, potentially including multiple volumes and choosing
# the one where space is available, and then mark the observations as 'stored'.
# If insufficient space is available in either pool, will exit with an error.
#
# Usage: storage_hoister.py /path/to/input/files /path/to/first/flat/disk
#           /path/to/first/bumpy/disk /path/to/second/bumpy/disk ...

P = ap.ArgumentParser(description="Copy scanified files to two sets of storage disks")
P.add_argument("bufferpath", help="Path to input buffer directory.  Observations will be copied from bufferpath/source/observation")
P.add_argument("flatoutpath", help="Path to flat disk for storage")
P.add_argument("bumpy_outdiskpaths", nargs='+', help="Path to bump disk(s) for storage")
P.add_argument("--bumpy-database", default="bumpy_storage.txt")
P.add_argument("--scanify-database", default="/data/autoproc/db/scanify.txt")
P.add_argument("--rate", choices=['fullrate', 'downsampled'], help='Observation data rate. If not given, all found observations will be processed')
P.add_argument("--source", help="Observation source.  If not given, all found observations will be processed")
P.add_argument("--observation", type=int, help="Observation ID.  If not given, all found observations will be processed")
P.add_argument("--offline-cal-root", help="Directory for offline calibration frames.  If set, a symlink will be created "
               "in the observation directory to the appropriate calibration frame.")
args = P.parse_args()

if not os.path.isabs(args.bumpy_database):
    args.bumpy_database = os.path.join(args.flatoutpath, args.bumpy_database)

if args.source is None and args.observation is None and args.rate is None:
    print('Opening scanify DB')
    scanify_db = status_db.ScanifyDatabase(args.scanify_database, host='storage01.spt')
    rates = ['fullrate', 'downsampled']
else:
    scanify_db = None
    if args.rate is None:
        raise ValueError("rate argument required")
    rates = [args.rate]
    if args.source is None:
        raise ValueError("source argument required")
    if args.observation is None:
        raise ValueError("observation argument required")

if len(args.bumpy_outdiskpaths) == 1:
    args.bumpy_outdiskpaths = sorted(glob.glob(args.bumpy_outdiskpaths[0]))

print('Opening bumpy disk DB')
bumpy_db = status_db.BumpyStorageDatabase(args.bumpy_database, host='storage01.spt')
print('Done opening DBs')

def copy_files_to_disk(inpath, outpath_leaf, outpath_base_options, headroom=5*(1024**3)):
    '''
    Copy the directory tree specified by inpath to the first of
    outpath_base_options with enough free space to hold the data,
    allowing at least headroom extra free space (by default, 5 GB).
    The output will be placed at the path outpath_base/outpath_leaf.
    Returns the selected outpath_base on exit.
    '''
    # Get total size of data
    total = 0
    for i in os.walk(inpath):
        for f in i[2]:
            total += os.stat(os.path.join(i[0], f)).st_size

    # Find a disk with space
    for disk in outpath_base_options:
        st = os.statvfs(disk)
        freebytes = st.f_bsize*st.f_bfree
        if freebytes <= total + headroom:
            continue

        try:
            # Now make that copy!
            outpath = os.path.join(disk, outpath_leaf)
            paroutpath = os.path.join(*os.path.split(outpath)[:-1])
            print('Copying %s to %s' % (inpath, outpath))
            if not os.path.exists(paroutpath):
                os.makedirs(paroutpath)
    
            # Delete old data if we have been asked to overwrite it
            if os.path.exists(outpath):
                shutil.rmtree(outpath)
            shutil.copytree(inpath, outpath)

        except Exception as e:
            # Handle errors gracefully
            if disk == outpath_base_options[-1]:
                raise
            warnings.warn('Error archiving %s to %s: %s.\nTrying another disk.'
                          % (inpath, outpath, str(e)))
            continue
        else:
            return disk

    else:
        raise IOError('No disk with the %.2f GB available space needed to copy %s' % (total * 1. / 1024**3, inpath))

for subdir in rates:

    k = 'status_{}'.format(subdir)

    if scanify_db is not None:
        queued = scanify_db.get_entries(**{k: 'buffered'})
    else:
        queued = [(args.source, args.observation)]

    for src, obs in queued:

        inpath = os.path.join(args.bufferpath, subdir, src, str(obs))

        print('Storing observation %s/%d' % (src, obs))

        if not os.path.exists(inpath):
            if scanify_db is None:
                raise IOError('Data for observation %s/%d is not currently on disk!' % (src, obs))
            else:
                warnings.warn('Data for observation %s/%d marked buffered but not currently on disk!' % (src, obs))
                continue

        subpath = os.path.join(subdir, src, str(obs))

        # First to the flat storage
        copy_files_to_disk(inpath, subpath, [args.flatoutpath])

        if args.offline_cal_root:
            offline_cal_file = os.path.join(args.offline_cal_root, src, str(obs) + '.g3')
            offline_cal_link = os.path.join(args.flatoutpath, subpath, 'offline_calibration.g3')
            if not os.path.islink(offline_cal_link):
                sp.check_call(['ln', '-s', offline_cal_file, offline_cal_link])

        # Then to the bumpy storage -- result path to bumpy DB
        bumpydisk = copy_files_to_disk(inpath, subpath, args.bumpy_outdiskpaths)

        # Files copied -- delete the originals
        shutil.rmtree(inpath)

        # And mark completion
        bumpy_db.update(src, obs, disk=bumpydisk, commit=True)
        if scanify_db is not None:
            scanify_db.update(src, obs, commit=True, **{k: 'stored'})

bumpy_db.close()
if scanify_db is not None:
    scanify_db.close()
