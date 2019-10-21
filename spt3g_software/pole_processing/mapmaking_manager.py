#!/usr/bin/env python

from __future__ import print_function
import subprocess as sp
import argparse as ap
import os, glob, warnings, shlex, shutil, re, pwd

P = ap.ArgumentParser(description="Find, process and transfer field maps.",
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('--input-root', default='/spt_data/bolodata/fullrate',
               help='Directory of scanified observation data')
P.add_argument('--output-root', default='/poleanalysis/sptdaq/calresult/onlinemaps',
               help='Directory for storing output map frames')

P.add_argument('--database', default='/data/autoproc/db/onlinemaps.txt',
               help='Map database file')
P.add_argument('--autoproc-database', default='/data/autoproc/db/autoproc.txt',
               help='Autoprocessing database')
P.add_argument('--transfer-database', default='/data/autoproc/db/aux_transfer.txt',
               help='Transfer database')

E = P.add_mutually_exclusive_group()
E.add_argument('--rerun-error', default=False, action='store_true',
               help='Rerun observations with status=error')
E.add_argument('--permafail-error', default=False, action='store_true',
               help='Mark all observations with status=error as '
               'permanently failed.  Use instead of --rerun-error for observations '
               'with known and unfixable problems.')

P.add_argument('--condor-log-root', default='/data/autoproc/logs/onlinemaps',
               help='Directory for storing Condor logs and submit scripts')
P.add_argument('--maxjobs', default=20, type=int,
               help='Maximum number of jobs to submit at once')

args = P.parse_args()

# TODO: handle coadds

from spt3g.std_processing.transfer_tools import LockedProcess, LockedFileError
try:
    lock = LockedProcess(os.path.abspath(__file__), timeout=5)
except LockedFileError as e:
    raise SystemExit

from spt3g.std_processing import status_db as sdb
from spt3g import cluster

print('*' * 80)
sdb.print_time('Starting onlime mapmaking manager')

db = sdb.AutoprocDatabase(args.database, verbose=True, sync=False)
adb = sdb.AutoprocDatabase(args.autoproc_database, verbose=True, read_only=True)

# sync observations
sdb.print_time('Syncing with autoprocessing database')
obs_entries = adb.get_entries('{}/calframe'.format(sdb.field_regex), match_regex=True)
for src, obs in obs_entries:
    src = src.split('/')[0]
    db.drop(src, obs)
    for band in ['90GHz', '150GHz', '220GHz']:
        bsrc = '/'.join([src, band])
        if obs < 62679600:
            # don't bother with observations prior to 2019 season
            db.update(bsrc, obs, status='permafail')
        else:
            db.update(bsrc, obs)
db.commit()

# check on status of jobs already in progress
sdb.print_time('Checking status of running jobs')
complete = []
for src, obs in db.get_entries(status=['submitted', 'executing', 'stopped']):
    condor_file_dir = os.path.join(args.condor_log_root, src, str(obs))
    logfile = os.path.join(condor_file_dir, 'condor.log')
    jid = int(db.match(src, obs, return_df=True)['job_id'])
    status, statusline = cluster.parse_job_status(logfile, jid)

    if status is None:
        sdb.print_warn("Unknown job status for %s/%d: %s" % (src, obs, statusline))
    else:
        db.update(src, obs, status=status)

    if status == 'complete':
        complete.append((src, obs))
    elif status == 'failed':
        errlog = logfile.replace('condor.log', 'condor.err')
        txt = '*****\n{}*****\n'.format(cluster.parse_error_log(errlog))
        sdb.print_warn('Job %s/%d failed. Checking error log %s:\n%s'
                       % (src, obs, errlog, txt))
db.commit()

# add completed jobs to transfer queue
if len(complete):
    tdb = sdb.AuxTransferDatabase(args.transfer_database, verbose=True)
    for bsrc, obs in complete:
        # if obs < 74803353:
        if obs < 70000000:
            # don't upload maps prior to sunset 2019 for now
            continue
        src, band = bsrc.split('/')
        fname = '{}/{}_{}_tonly.g3.gz'.format(src, obs, band)
        tdb.update(fname, type='map')
    tdb.commit()
    tdb.close()

# add older completed jobs to transfer queue
if 0:
    tdb = sdb.AuxTransferDatabase(args.transfer_database, verbose=True)
    for bsrc, obs in db.get_entries(status='complete')[::-1]:
        if obs < 63000000 or obs >= 74803353:
            continue
        src, band = bsrc.split('/')
        fname = '{}/{}_{}_tonly.g3.gz'.format(src, obs, band)
        tdb.update(fname, type='map')
    tdb.commit()
    tdb.close()

# handle failed jobs
for src, obs in db.get_entries(status='failed'):
    if args.permafail_error:
        db.update(src, obs, status='permafail')
    elif args.rerun_error:
        pass
    else:
        sdb.print_warn(
            "Job %s/%d has failed.  Rerun or mark as permanently failed." % (src, obs)
        )

if args.maxjobs == 0:
    raise SystemExit

# submit new jobs
job_template = """
executable = onlinemaps_job.sh
arguments = "\'{band}\' \'{output}\' \'{inputs}\'"
log = {logdir}/condor.log
output = {logdir}/condor.out
error = {logdir}/condor.err
notification = never
request_memory = 30GB
getenv = False
run_as_owner = True
transfer_executable = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
queue
"""

os.chdir(os.path.join(os.getenv('SPT3G_SOFTWARE_PATH'), 'pole_processing'))

sdb.print_time('Processing new jobs')
job_count = len(db.get_entries(status=['submitted', 'executing']))
for bsrc, obs in db.get_entries(status='failed' if args.rerun_error else None)[::-1]:
    if args.rerun_error:
        db.update(bsrc, obs, status=None)

    # process maps in reverse chronological order to prioritize new data
    # don't submit too many jobs at once
    if job_count >= args.maxjobs:
        if args.rerun_error:
            continue
        break

    src, band = bsrc.split('/')

    # find input files
    input_root = os.path.join(args.input_root, src, str(obs))
    calframe = os.path.realpath(os.path.join(input_root, 'offline_calibration.g3'))
    if not os.path.exists(calframe):
        continue
    if adb.match('{}/calframe'.format(src), obs, return_df=True, squeeze=True)['status'] != 'complete':
        continue

    # populate job specs
    inputs = sorted(glob.glob(os.path.join(input_root, '[0-9]*.g3')))
    if not len(inputs):
        continue
    inputs = [calframe] + inputs
    output = os.path.join(args.output_root, src, '{}_{}.g3'.format(obs, band))
    if not os.path.exists(os.path.dirname(output)):
        os.makedirs(os.path.dirname(output))
    logdir = os.path.join(args.condor_log_root, bsrc, str(obs))
    if not os.path.exists(logdir):
        os.makedirs(logdir)
    job_specs = job_template.format(
        band=band,
        inputs=' '.join(inputs),
        output=output,
        logdir=logdir,
    )

    # submit job
    jobid, result = cluster.condor_submit_baked(job_specs)

    # bookkeeping
    if jobid is None:
        sdb.print_warn("Error submitting job %s/%d: %s" % (src, obs, result))
    else:
        sdb.print_time("Submitted job %s/%d: %d" % (src, obs, jobid))
        db.update(bsrc, obs, job_id=jobid, status='submitted')
        job_count += 1

# cleanup
db.commit()
db.close()
print('*' * 80)
