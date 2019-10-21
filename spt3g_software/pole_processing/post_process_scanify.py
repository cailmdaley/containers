#!/usr/bin/env python

import os, glob, shutil
import argparse
import spt3g.std_processing.status_db as sdb
from spt3g import core
import pandas as pd

def get_file_list(directory, time, extension=''):
    '''
    Find files older than `time` in `directory`.  File times are parsed 
    from the filename (assuming the filenames are '%Y%m%d_%H%M%S.`extension`')
    INPUTS
    ------
    directory: str
        The directory in which to look for files.
    time: datetime
        Files older than this time are returned.
    extension: str ('')
        Only look for files with this extension.
    '''
    files = sorted(glob.glob(os.path.join(directory, '*' + extension)))
    file_list = []
    for f in files:
        timestamp, extension = os.path.basename(f).split('.')
        file_time = core.G3Time(timestamp)
        if file_time < time:
            file_list.append(f)

    # the last file may contain the start of an observation that needs 
    # to be processed by a later run of scanify_bolodata.
    if len(file_list) > 0:
        file_list.pop(-1) 

    return file_list

def find_files_to_del(db_name, ignore_after=None, 
                      bolodata_dir='/buffer/bolodata/'):
    '''
    Find bolodata files that are safe to delete.

    A file is considered safe to delete (at least by this function)
    if a continuous scanify jobs have been run including its stop time.
    INPUTS
    ------
    db_name: str
        Name of the file containing the job database.
    ignore_after: datetime (None)
        Ignore any files created after this time.  Useful if you want to
        delete only files older than 24 hours.  If None, no limit is used.
    bolodata_dir: str ('/buffer/bolodata/')
        Directory containing bolodata files.
    '''
    if ignore_after is None:
        ignore_after = core.G3Time.Now() - 24 * core.G3Units.hours

    db = sdb.ScanifyDatabase(db_name, read_only=True)
    last_stop = core.G3Time(db.get_latest_contiguous_run().isoformat())
    last_obs_stop = core.G3Time(db.get_latest_stop().isoformat())
    if last_obs_stop == last_stop:
        # if the contiguous good data include the last observation in the database,
        # then all data after that are acquired while the telescope is idle.
        # make sure not to overfill the buffer if we are idle for over a day, so
        # return the timestamp of the most recent idle period
        arc_stop_file = os.path.join(os.path.dirname(db_name), 'latest_arc_time')
        try:
            last_arc_stop = core.G3Time(open(arc_stop_file, 'r').read().strip())
            if last_stop < last_arc_stop:
                last_stop = last_arc_stop
        except:
            pass
    db.close()
    last_stop = ignore_after if ignore_after < last_stop else last_stop
    bolofiles = get_file_list(bolodata_dir, last_stop, '.g3')
    return last_stop, bolofiles


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--database', default = '/data/autoproc/db/scanify.txt',
                        help = 'The location of the database storing job information')
    parser.add_argument('--bolodir', default = '/buffer/bolodata',
                        help = 'The directory containing pre-scanified bolometer data')
    parser.add_argument('--keep-after', default=24, type=int,
                        help='Do not delete files newer than this many hours ago')
    parser.add_argument('--dry-run', default=False, action='store_true',
                        help='Print what will be done and exit without deleting files.')
    args = parser.parse_args()

    now = core.G3Time.Now()
    ignore = now - args.keep_after * core.G3Units.hours
    ignore, bolofiles = find_files_to_del(args.database, ignore_after=ignore,
                                          bolodata_dir=args.bolodir)
    sdb.print_time('Removing {} timepoint files older than {}'.format(len(bolofiles), ignore))
    if args.dry_run:
        if len(bolofiles):
            print('\n'.join(bolofiles))
        raise SystemExit

    for bf in bolofiles:
        print(bf)
        try:
            # shutil.move(bf, os.path.join('/poleanalysis/sptdaq/bolodata_deleteme_maybe/',
            #                          os.path.basename(bf)))
            os.remove(bf)
        except OSError:
            pass
        except shutil.Error:
            pass
        except IOError:
            pass
