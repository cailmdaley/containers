#!/usr/bin/env python

from __future__ import print_function
from spt3g.std_processing import status_db
from spt3g.cluster.condor_tools import parse_job_status, parse_error_log, condor_submit_baked
import sys, os, pandas, numpy, glob, subprocess, re, tempfile, platform, shutil
import argparse
from autoproc_utils import make_scripts_dict, make_obs_type_dict

# Usage scheduler.py /path/to/dbs /path/to/autoprocdb.txt /path/to/calresults /path/to/data (downsampled|fullrate) proxycert [failed] [-v|--verbose]

P = argparse.ArgumentParser(description="Autoprocessing scheduling with dependency tracking")
P.add_argument("db_root", help="Database root")
P.add_argument("autoproc_db", help="Path to autoprocessing database")
P.add_argument("calresults", help="Path to output directory")
P.add_argument("datapath", help="Path to input data directory")
P.add_argument("datatype", choices=["downsampled", "fullrate"],
               help="Type of data to process")
P.add_argument("proxycert", help="Globus certificate (should be \"\" at pole)")
P.add_argument("jobtype", nargs="?", choices=["failed"],
               help="Type of job to handle")
P.add_argument("-v", "--verbose", default=False, action="store_true",
               help="Print job submission status")
P.add_argument('--maxjobs', type=int, default=None,
               help='Maximum number of jobs to run at once')
args = P.parse_args()


db_root = args.db_root
autoproc_db = args.autoproc_db

calresults = args.calresults
datapath = args.datapath
datatype = args.datatype # 'downsampled' or 'fullrate'
proxycert = args.proxycert # Should be '' at pole
failedjobs = (args.jobtype == 'failed')

if args.maxjobs is not None and args.maxjobs < 0:
    args.maxjobs = numpy.inf

# Some switches in terms of file organization if running at pole
if platform.node() == 'anal.spt':
    atpole = True
    toolset = 'py3-v3'
    maxjobs = 20 if args.maxjobs is None else args.maxjobs
else:
    atpole = False
    toolset = 'py3-v2'
    maxjobs = numpy.inf if args.maxjobs is None else args.maxjobs

if atpole:
    from spt3g.std_processing import transfer_tools
    try:
        lock = transfer_tools.LockedProcess(os.path.abspath(__file__), timeout=5)
    except transfer_tools.LockedFileError as e:
        raise SystemExit

# Open DBs without locking if running in the north. A little racy if we catch
# them during transfer, but the rsync that sent them is racy anyway, so that
# can't be helped. We can't hold the locks because this is not on a shared
# filesystem with pole.
scanify_db = status_db.ScanifyDatabase(os.path.join(db_root, 'scanify.txt'), read_only=True)
if atpole:
    autoproc_db = status_db.AutoprocDatabase(autoproc_db, user='sptdaq',
      host='anal.spt')
else:
    autoproc_db = status_db.AutoprocDatabase(autoproc_db, user='nwhitehorn',
      host='scott.grid.uchicago.edu')

# Processing information: maps observation type to list of qualifying source names
# Used, e.g. to build calframes using calibration sources other htan RCW38
obs_types = make_obs_type_dict()

# Processing information: maps source name to processing function and
# (optional) dependencies list
scripts = make_scripts_dict(
    calresults=calresults,
    proxycert=proxycert,
    atpole=atpole,
    toolset=toolset,
)

# Now figure out the list of observations that exist and might need processing
available_observations = scanify_db.get_entries()

# Synchronize two databases
extant_autoproc_obs = set(autoproc_db.get_entries())
for src, obs in available_observations:
    for k in scripts.keys():
        if k.split('/')[0] == src and (k, obs) not in extant_autoproc_obs:
            autoproc_db.update(k, obs)
        if k.startswith('calframe') and (src + '/' + k, obs) not in extant_autoproc_obs:
            autoproc_db.update(src + '/' + k, obs)
autoproc_db.commit()

# Check up on in-progress jobs by spinning through log files
for src, obs in autoproc_db.get_entries(status=['submitted', 'executing', 'stopped']):
    logfile = os.path.join(calresults, 'logs', src, str(obs), 'condor.log')
    jid = int(autoproc_db.match(src, obs, return_df=True)['job_id'])
    status, statusline = parse_job_status(logfile, jid)
    if status is None:
        print('Unknown status for %s/%d: %s' % (src, obs, statusline), file=sys.stderr)
    else:
        autoproc_db.update(src, obs, status=status)
    if status == 'failed':
        print('Job %s/%d failed' % (src, obs), file=sys.stderr)
        errorlog_path = os.path.join(calresults, 'logs', src, str(obs), 'error.log')
        try:
            err = parse_error_log(errorlog_path)
        except FileNotFoundError:
            err = 'Error log file {} not found'.format(errorlog_path)
        if err:
            print('***** {} *****\n'.format(errorlog_path), file=sys.stderr)
            print(err, file=sys.stderr)
            print('*****\n', file=sys.stderr)
autoproc_db.commit()

# Submit any jobs that need submitting
checkstatus = None
if failedjobs:
    checkstatus = 'failed'

if args.verbose:
    print("*" * 50)
    os.system("date -u")

job_count = len(autoproc_db.get_entries(status=['submitted', 'executing']))

for src, obs in autoproc_db.get_entries(status=checkstatus):
    if job_count >= maxjobs:
        print("Submitted {} jobs, stopping for now".format(job_count))
        break

    if src in scripts:
        job = scripts[src]
    elif src.split('/')[-1] in scripts:
        job = scripts.get(src.split('/')[-1])
    else:
        job = scripts['default']

    if job is None:
        continue

    job_script, job_deps = job

    # Check if dependencies are there
    deps = []
    todmissing = False
    if atpole:
        valid_status = ['ready', 'uploading', 'uploaded', 'sent', 'verified']
    else:
        valid_status = ['verified']
    for dep in job_deps:
        if dep == 'tod':
            # Check if data exist
            obs_src = src.split('/')[0]
            df  = scanify_db.match(obs_src, obs, return_df=True)
            # handle multiple or missing entries in scanifier database (weird edge cases)
            if len(df) > 1:
                raise ValueError("Multiple entries found for observation %s/%d:\n%s"
                                 % (obs_src, obs, df))
            if len(df) == 0:
                raise ValueError("Observation %s/%d not found" % (obs_src, obs))
            # otherwise expect a single database entry
            if not df['status_' + datatype].isin(valid_status).all():
                todmissing = True
            if todmissing:
                break
            else:
                continue

        if dep.startswith('calframe'):
            dep = '%s/%s' % (src.split('/')[0], dep)

        if dep in obs_types:
            # dependency is a type of observation (one of a list of source names)
            deplist = obs_types[dep]
        else:
            # dependency is a specific source name
            deplist = [dep]

        deps_avail = \
          autoproc_db.data[numpy.logical_and(autoproc_db.data['source'].isin(deplist),
                                             autoproc_db.data['observation'] <= obs)]
        if len(deps_avail) == 0 or deps_avail['status'].values[-1] != 'complete':
            break
        deps.append(os.path.join(str(deps_avail['source'].values[-1]), str(deps_avail['observation'].values[-1])))
    if len(deps) != len([d for d in job_deps if d != 'tod']) or todmissing:
        # Didn't find some job dependency, so bail
        if checkstatus == 'failed':
            print('Not submitting job %s/%d because of missing dependencies. Needed %s, have %s.' % (src, obs, job_deps, deps), file=sys.stderr)
        continue

    files = []
    for d in deps:
        if 'calframe' in d:
            files.append(os.path.join(calresults, 'calibration', d.split('/')[1], src.split('/')[0], str(obs) + '.g3'))
        else:
            files.append(os.path.join(calresults, 'calibration', d + '.g3'))
    if 'tod' in job_deps:
        todfiles = sorted(glob.glob(os.path.join(datapath, datatype, src.split('/')[0], str(obs), '[0-9]*.g3')))
        if len(todfiles) == 0:
            print('Missing TOD for %s/%d!' % (src.split('/')[0], obs), file=sys.stderr)
            continue
        files += todfiles

    # Start job
    if args.verbose:
        print(src, obs)
        print('\t', files)
    jid = job_script(src, obs, files)
    if isinstance(jid, str):
        if args.verbose:
            print('\tJob locally completed: %s' % jid)
        autoproc_db.update(src, obs, status=jid, job_id=0)
    elif jid is not None:
        if args.verbose:
            print('\tJob %d' % jid)
        autoproc_db.update(src, obs, status='submitted', job_id=jid)
        job_count += 1

autoproc_db.close()
autoproc_db.write_sql()

if args.verbose:
    print("*" * 50)
