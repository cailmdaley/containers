#!/usr/bin/env python

from __future__ import print_function
import subprocess as sp
import argparse as ap
import os, glob, warnings, shlex, shutil, re, pwd

P = ap.ArgumentParser(description='Find, process and transfer observation data.',
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
S = P.add_subparsers(dest='mode', metavar='MODE', title='subcommands',
                     help='Function to perform. For more help call: '
                     '%(prog)s %(metavar)s -h')
parser_opts = dict(formatter_class=ap.ArgumentDefaultsHelpFormatter)

modes = [('update', 'Process any available observations'),
         ('adduser', 'Add a manual (userdata) observation to be processed'),
         ('pre', 'Pre-process a DAG job for a particular observation'),
         ('run', 'Run a DAG job for a particular observation'),
         ('post', 'Post-process a DAG job for a particular observation')]

for mode, helptext in modes:

    PP = S.add_parser(mode, help=helptext, **parser_opts)

    PP.add_argument('--buffer-input', default='/buffer/bolodata',
                    help='Directory containing raw timepoint data')
    PP.add_argument('--arc-input', default='/spt_data/arc',
                    help='Directory containing arcfiles')
    PP.add_argument('--buffer-output', default='/buffer',
                    help='Base directory for scanified output observation files')
    PP.add_argument('--storage-output', default='/spt_data/bolodata')
    PP.add_argument('--database', default='/data/autoproc/db/scanify.txt',
                    help='Observation database file')
    PP.add_argument('--autoproc-root', default='/poleanalysis/sptdaq/calresult',
                    help='Directory for autoprocessing output frames')

    if mode == 'adduser':
        PP.add_argument('start', help='Observation start time in UTC')
        PP.add_argument('stop', help='Observation stop time in UTC')

    if mode == 'update':
        E = PP.add_mutually_exclusive_group()
        E.add_argument('--rerun-error', default=False, action='store_true',
                       help='Rerun observations with status_scanify=error')
        E.add_argument('--permafail-error', default=False, action='store_true',
                       help='Mark all observations with status_scanify=error as '
                       'permanently failed.  Use instead of --rerun-error for observations '
                       'with known and unfixable problems.')
        PP.add_argument('--force', default=False, action='store_true',
                        help='Force rerun observations to overwrite existing rescue files. '
                        'Use this option if the DAG or job files have been modified.  '
                        'Use with --rerun-error or when manually clearing failed observations.')
        PP.add_argument('--rerun-processing', default=False, action='store_true',
                        help='Rerun observations that are stuck with status_scanify=processing. '
                        'This happens when the scanify process is interrupted (e.g. reboots)')
        PP.add_argument('--condor-log-root', default='/data/autoproc/logs/scanify',
                        help='Directory for storing Condor logs and submit scripts')
        PP.add_argument('--transfer-fullrate', default=False, action='store_true',
                        help='Transfer all fullrate data by default')
        PP.add_argument('--transfer-fullrate-cal', default=False, action='store_true',
                        help='Transfer all non-field observations at fullrate')
        PP.add_argument('--transfer-fullrate-fields', default=3, type=int,
                        help='Transfer every Nth field scan at fullrate (0=no transfer)')
        PP.add_argument('--maxjobs', default=10, type=int,
                        help='Maximum number of jobs to submit at once')
        PP.add_argument('--no-find', default=True, action='store_false', dest='find',
                        help='Do not look for new observations to add to the database')
        PP.add_argument('--autoproc-database', default='/data/autoproc/db/autoproc.txt',
                        help='Autoprocessing database file')

    if mode in ['pre', 'run', 'post']:
        PP.add_argument('job_name', help='Name of DAG job to process')
        PP.add_argument('source', help='Observation source')
        PP.add_argument('observation', type=int, help='Observation ID')
        PP.add_argument('--start', help='Observation start time')
        PP.add_argument('--stop', help='Observation stop time')
        PP.add_argument('--downsample-only', default=False, action='store_true',
                        help='Assume fullrate data have already been processed.')
        PP.add_argument('--storage-error', default=False, action='store_true',
                        help='Assume fullrate data have been processed but a problem '
                        'occurred archiving the data to permanent storage.  Use with '
                        '--rerun-error')
        PP.add_argument('--keep-q', default=False, action='store_true',
                        help='Keep Q-phase data when downsampling')
        PP.add_argument('--calframe',
                        help='Calibration frame for rotating the IQ phase when downsampling')
        PP.add_argument('--user-obs', default=False, action='store_true',
                        help="Construct an observation with source name 'userdata' "
                        "and obs ID determined from the start time")

    if mode == 'post':
        PP.add_argument('job_pre_return', type=int, help='Return status of the job pre script')
        PP.add_argument('job_return', type=int, help='Return status of the job')

args = P.parse_args()

# construct common CLI for this script
script_file = os.path.abspath(os.path.join(os.getenv('SPT3G_SOFTWARE_PATH'), 'pole_processing/scanify_manager.py'))
script_path = os.path.dirname(script_file)

script_args = ('--buffer-input {} --arc-input {} --buffer-output {} --storage-output {} '
               '--autoproc-root {} --database {} {{rerun_opts}} {{keep_q}} '
               '{{calframe}} {{user_obs}}'.format(
                   args.buffer_input, args.arc_input, args.buffer_output,
                   args.storage_output, args.autoproc_root, args.database))

pyexec = sp.check_output(['which', 'python']).strip().decode()

# Observations with this source name will be treated as fake
# i.e. created with the given source name and ObservationID set to the given start time,
# even if GCP claims there were no observations during this period.
user_source_name = 'userdata'

# construct command list
scanify_cmd = ('{} {}/scanify_bolodata.py {}/fullrate {} {} --scratch '
               '--start-time {{start}} --stop-time {{stop}} {{user_obs}}'.format(
                   pyexec, script_path, args.buffer_output, args.arc_input, args.buffer_input))

downsample_cmd = ('{} {}/downsample_bolodata.py {}/fullrate/{{source}}/{{observation}} '
                  '{}/downsampled/{{source}}/{{observation}} --scratch {{keep_q}} '
                  '{{calframe}}'.format(
                      pyexec, script_path, args.storage_output, args.buffer_output))

store_cmd = ('ssh storage01.spt \"source ~/.bash_profile; {} {}/storage_hoister.py {} {} {} '
             '--rate {{rate}} --source {{{{source}}}} --observation {{{{observation}}}} '
             '--offline-cal-root {}\"'.format(
                 pyexec, script_path, args.buffer_output, '/spt_disks/primary/bolodata/',
                 '/spt_disks/S*/', os.path.join(args.autoproc_root, 'calibration', 'calframe')))
store_cmd_full = store_cmd.format(rate='fullrate')
store_cmd_down = store_cmd.format(rate='downsampled')

check_cmd = ('{} {}/check_bolodata.py {{rate}} {{{{source}}}} {{{{observation}}}} '
             '--data-root {}'.format(pyexec, script_path, args.storage_output))
check_cmd_full = check_cmd.format(rate='fullrate')
check_cmd_down = check_cmd.format(rate='downsampled')

scripts = {
    'scanify': {'cmd': scanify_cmd,
                'success': {'status_fullrate': 'buffered'},
                'error': {'status_fullrate': 'scanify_error'}},
    'store_full': {'cmd': store_cmd_full,
                   'success': {'status_fullrate': 'stored'},
                   'error': {'status_fullrate': 'storage_error'}},
    'check_full': {'cmd': check_cmd_full,
                   'success': {'status_fullrate': 'ready'},
                   'error': {'status_fullrate': 'check_error'}},
    'downsample': {'cmd': downsample_cmd,
                   'pre': True,
                   'success': {'status_downsampled': 'buffered'},
                   'error': {'status_downsampled': 'scanify_error'}},
    'store_down': {'cmd': store_cmd_down,
                   'success': {'status_downsampled': 'stored'},
                   'error': {'status_downsampled': 'storage_error'}},
    'check_down': {'cmd': check_cmd_down,
                   'success': {'status_downsampled': 'ready'},
                   'error': {'status_downsampled': 'check_error'}},
}


# construct job and DAG file templates to be used for each observation
job_orders = {
    'fullrate': ['scanify', 'store_full', 'check_full'],
    'downsampled': ['downsample', 'store_down', 'check_down'],
}

def write_dag_files(job_order, file_root, **dag_opts):
    job_template = """
Executable = $(job_exec)
arguments = $(job_args)

log_root = scanify_{source}_{observation}_{rate}_$(job_name)_$(Cluster).$(Process)
Log = $(log_root).log
Output = $(log_root).out
Error = $(log_root).err

Notification = never
Request_Memory = {memory}GB
getenv = True

run_as_owner = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
priority = 100

Queue
    """

    dag_template = []
    for idx, job in enumerate(job_order):
        script = scripts[job]

        dag_template.append('JOB {} {{condor_file_root}}.job'.format(job))

        if script.get('pre', False):
            dag_template.append('SCRIPT PRE {} {} {} pre $JOB {{source}} {{observation}} '
                                '--start {{start}} --stop {{stop}} '
                                '{}'.format(job, pyexec, script_file, script_args))
        if script.get('post', False):
            dag_template.append(
                'SCRIPT POST {} {} {} post $JOB {{source}} {{observation}} '
                '--start {{start}} --stop {{stop}} $PRE_SCRIPT_RETURN $RETURN '
                '{}'.format(job, pyexec, script_file, script_args))

        if idx > 0:
            dag_template.append('PARENT {} CHILD {}'.format(job_order[idx-1], job))

    dag_template.append('VARS ALL_NODES job_name=\"$(JOB)\"')
    dag_template.append('VARS ALL_NODES job_exec=\"{}\"'.format(pyexec))
    dag_template.append('VARS ALL_NODES job_args=\"{} run $(JOB) {{source}} {{observation}} '
                        '--start {{start}} --stop {{stop}} {}\"'.format(script_file, script_args))

    dag_template.append('NODE_STATUS_FILE {condor_file_root}.dag.status')
    dag_template = '\n'.join(dag_template)

    with open(file_root + '.job', 'w') as f:
        f.write(job_template.format(**dag_opts))
    with open(file_root + '.dag', 'w') as f:
        f.write(dag_template.format(**dag_opts))


if args.mode == 'adduser':

    print('*'*80)

    from spt3g.std_processing import status_db as sdb
    from spt3g.std_processing import time_to_obsid
    from spt3g.core import G3Time
    from pandas import Timestamp

    start_time = G3Time(args.start).isoformat()
    obsid = time_to_obsid(start_time)
    start_time = Timestamp(start_time, tz='UTC')
    stop_time = Timestamp(G3Time(args.stop).isoformat(), tz='UTC')
    if stop_time < start_time:
        raise ValueError(
            "Invalid observation times! Start: {} Stop: {}".format(start_time, stop_time)
        )

    db = sdb.ScanifyDatabase(args.database, verbose=True)
    db.update(user_source_name, obsid,
              obs_start=start_time, obs_stop=stop_time,
              status_fullrate=None, status_downsampled=None,
              status_scanify=None, commit=True)
    db.close()

    print('*'*80)
    raise SystemExit

if args.mode == 'update':

    print('*'*80)

    from spt3g.std_processing import status_db as sdb
    from spt3g.std_processing.transfer_tools import LockedProcess, LockedFileError
    from spt3g.cluster.condor_tools import parse_dag_status, parse_error_log
    from spt3g.core import G3Time, G3Units
    import pandas as pd
    import numpy as np

    sdb.print_time('Starting scanify manager')

    try:
        lock = LockedProcess(os.path.abspath(__file__), lock=False)
    except LockedFileError as e:
        print('Scanify manager is locked by another process')
        raise SystemExit

    # Simplify arguments
    if args.permafail_error or args.rerun_error:
        args.find = False

    # update observation database
    try:
        if args.find:
            # read latest processed observation time from file
            # the finder updates this file with the most recent non-observing frame
            # to minimize I/O time on useless data
            first_obs = os.path.join(os.path.dirname(args.database), 'latest_arc_time')
            last_obs = str(G3Time.Now() - 30 * G3Units.minute)

            # if the above file is missing, start by searching after the end of the
            # most recent observation
            if not os.path.exists(first_obs):
                first_obs = 'latest'

            out = sp.check_output(
                '{} {}/find_observations.py {} --first-obs {} --last-obs {} --database {}'
                .format(
                    pyexec, script_path, args.arc_input, first_obs, last_obs, args.database
                ).split(),
                stderr=sp.STDOUT
            ).decode()
            print(out)
            err = parse_error_log(out)
            if err:
                sdb.print_warn('Error in observation finder:\n{}'.format(err))

    except sp.CalledProcessError as e:
        m = re.search('(\S+) truncated; unable to read register map', e.output.decode() if e.output else "")
        if m:
            fname = os.path.basename(m.group(1))
            sdb.print_warn('Found truncated arcfile {}, moving out of the way'.format(fname))
            cmd = 'ssh storage01 mv /spt_disks/primary/arc/{} /spt_disks/primary/arc/bad/.'
            sp.check_call(cmd.format(fname), shell=True)
        else:
            sdb.print_warn('Error in observation finder', e)

    # XXX handle missing G3 data very loudly, but need to make sure enough time is given
    # for the rsync from bolocontrol to complete

    db = sdb.ScanifyDatabase(args.database, verbose=True)
    df_update = db.match(status_scanify=[None, 'error', 'rerun', 'processing'], return_df=True)
    autoproc_db = None

    num_jobs = 0

    for idx, entry in df_update.iterrows():

        src = entry['source']
        obs = entry['observation']

        if entry['status_fullrate'] in ['ready', 'uploading', 'uploaded', 'sent', 'verified']:
            rate = 'downsampled'
        else:
            rate = 'fullrate'

        job_order = job_orders[rate]
        condor_file_dir = os.path.join(args.condor_log_root, src, str(obs))
        condor_file_root = os.path.join(condor_file_dir, 'scanify_{}_{}_{}'.format(src, obs, rate))
        condor_job_file = condor_file_root + '.job'
        condor_dag_file = condor_file_root + '.dag'
        condor_dag_status = condor_dag_file + '.status'

        if entry['status_scanify'] == 'processing':

            # Restart jobs that are stuck in processing (e.g. after computers are rebooted)
            if args.rerun_processing:
                db.update(src, obs, status_scanify='rerun')
                continue

            # Don't process new jobs if rerunning error'ed ones
            if args.permafail_error or args.rerun_error:
                continue

            # status file doesn't exist yet
            if not os.path.exists(condor_dag_status):
                continue

            # read in status file
            status = parse_dag_status(condor_dag_status)

            # check individual node status
            update_opts = {}
            nodes = status['NodeStatus']
            job_error = False

            for job in job_order:
                key = list(scripts[job]['success'].keys())[0]
                if entry[key] in ['uploading', 'uploaded', 'sent', 'verified']:
                    continue
                if nodes[job]['NodeStatus'] == 'done':
                    update_opts.update(**scripts[job]['success'])
                elif nodes[job]['NodeStatus'] == 'error':
                    job_error = True
                    update_opts.update(**scripts[job]['error'])
                    break

            # check overall DAG status
            if status['DagStatus']['DagStatus'] == 'done':
                update_opts.update(status_scanify='success' if rate == 'downsampled' else None)
            elif status['DagStatus']['DagStatus'] == 'error':
                job_error = True
                update_opts.update(status_scanify='error')

            # transfer logic
            if entry['status_fullrate'] == 'ready' or \
               update_opts.get('status_fullrate') == 'ready':
                xfer = False
                if args.transfer_fullrate:
                    xfer = 1
                elif args.transfer_fullrate_fields and \
                     re.match(sdb.field_regex, src):
                    # find the last field observation marked for transfer
                    field_obs = db.match(src, status_scanify='success', return_df=True)
                    idx = np.where(field_obs['transfer_fullrate'] > 0)[0]
                    # transfer every Nth
                    if not len(idx) or idx[-1] <= len(field_obs) - args.transfer_fullrate_fields:
                        # prioritize fullrate calibration data over field scans
                        xfer = 2 if args.transfer_fullrate_cal else 1
                elif args.transfer_fullrate_cal:
                    xfer = 1
                if xfer:
                    update_opts.update(transfer_fullrate=xfer)

            if entry['status_downsampled'] == 'ready' or \
               update_opts.get('status_downsampled') == 'ready':
                xfer = False
                if not args.transfer_fullrate:
                    if args.transfer_fullrate_fields and \
                       re.match(sdb.field_regex, src):
                        # transfer all downsampled field observations
                        xfer = 1
                    elif not args.transfer_fullrate_cal:
                        xfer = 1
                if xfer:
                    update_opts.update(transfer_downsampled=xfer)

            # drop options that are already set
            for k in list(update_opts.keys()):
                if entry[k] == update_opts[k]:
                    update_opts.pop(k)

            # trigger a cron email on errors
            if job_error:
                dag_mtime = os.stat(condor_dag_file).st_mtime
                txt = ''
                for f in glob.glob('{}_*.err'.format(condor_file_root)):
                    bf = os.path.basename(f)
                    st = os.stat(f)
                    if st.st_size == 0:
                        continue
                    if st.st_mtime < dag_mtime:
                        continue
                    # only print the relevant lines from the error file
                    err = parse_error_log(f)
                    if err:
                        txt += '\n***** {} *****\n'.format(bf)
                        txt += err
                        txt += '***** {} *****\n'.format(bf)
                    # automatically rerun scanifier jobs that run into missing data
                    # within one hour of the end of the observation, to give time for the
                    # DAQ rsync to catch up.
                    if 'Rsync from daqcontrol is not up to date' in txt:
                        delta = pd.Timestamp(st.st_mtime, unit='s', tz='utc') - entry['obs_stop']
                        if delta < pd.Timedelta(1, unit='h'):
                            update_opts.update(status_scanify='rerun')
                            job_error = False
            if job_error:
                txt = 'Error {} observation {}/{}, checking logs in {}:\n{}\n'.format(
                    'scanifying' if rate == 'fullrate' else 'downsampling',
                    src, obs, condor_file_dir, txt)
                txt += '*' * 80
                sdb.print_warn(txt)

            if update_opts:
                db.update(src, obs, **update_opts)

        elif num_jobs < args.maxjobs:

            if entry['status_scanify'] == 'error':
                if args.permafail_error:
                    db.update(src, obs, status_scanify='permafail')
                    continue
                elif not args.rerun_error:
                    sdb.print_warn(
                        'Observation {}/{} failed to {}.  Rerun or mark as permanently failed.'
                        .format(src, obs, 'scanify' if rate == 'fullrate' else 'downsample')
                    )
                    continue
                else:
                    entry['status_scanify'] = 'rerun'
            elif args.permafail_error or args.rerun_error:
                # Don't process new jobs if rerunning error'ed ones
                continue

            # check for downsampling dependencies in autoprocessing database
            calframe = None
            if rate == 'downsampled':
                if autoproc_db is None:
                    autoproc_db = sdb.AutoprocDatabase(args.autoproc_database, verbose=False, read_only=True)

                df = autoproc_db.match('elnod', return_df=True)
                df = df[df['observation'] <= obs]
                cal_entry = df.loc[df.index.values[-1]]
                if cal_entry['status'] == 'complete':
                    calframe = os.path.join(args.autoproc_root, 'calibration', cal_entry['source'], str(cal_entry['observation']) + '.g3')
                else:
                    sdb.print_time('Observation {}/{} calframe dependency is not ready.'.format(src, obs))
                    continue

            if not os.path.exists(condor_file_dir):
                os.makedirs(condor_file_dir)

            if src == 'debug-forced-scanify':
                memory = 40
            elif src == 'noise':
                memory = 30
            elif src == 'calibrator':
                memory = 25
            else:
                memory = 10

            rerun_opts = ''
            if entry['status_scanify'] == 'rerun':
                if (entry['status_fullrate'] in
                    ['ready', 'uploading', 'uploaded', 'sent', 'verified']):
                    rerun_opts = '--downsample-only'
                elif entry['status_fullrate'] in ['buffered', 'storage_error']:
                    rerun_opts = '--storage-error'

            if src == 'elnod':
                keep_q = '--keep-q'
            else:
                keep_q = ''

            if src == user_source_name:
                user_obs = '--user-obs'
            else:
                user_obs = ''

            if calframe:
                calframe = '--calframe {}'.format(calframe)
            else:
                calframe = ''

            dag_opts = dict(source=src, observation=obs, rate=rate,
                            start=entry['obs_start'].isoformat(),
                            stop=entry['obs_stop'].isoformat(),
                            condor_file_root=condor_file_root,
                            memory=memory, rerun_opts=rerun_opts,
                            keep_q=keep_q, calframe=calframe,
                            user_obs=user_obs)
            write_dag_files(job_order, condor_file_root, **dag_opts)

            try:
                out = sp.check_output(['condor_submit_dag', '-usedagdir']
                                      + ['-priority', '100']
                                      + ['-force'] * (args.force or entry['status_scanify'] == 'rerun')
                                      + [condor_dag_file],
                                      stderr=sp.STDOUT).decode()
                # sdb.print_time(out)
                m = glob.glob('{}.rescue*'.format(condor_dag_file))
                if len(m):
                    sdb.print_time(
                        'Rerunning {} observation {}/{}, attempt {}'.format(
                            rate, src, obs, len(m)))
                if "Warning: FindLastRescueDagNum() hit maximum" in out:
                    sdb.print_warn(
                        'Observation {}/{} re-running too often, attempt {}'.format(
                            src, obs, len(m)))
            except sp.CalledProcessError as e:
                sdb.print_warn('Error submitting {}'.format(condor_dag_file), e)
            else:
                # update database
                sdb.print_time('Submitted {} jobs for {}/{}'.format(rate, src, obs))
                db.update(src, obs, status_scanify='processing')
                num_jobs += 1

    db.commit()
    db.close()
    print('*'*80)
    raise SystemExit

# Now run each individual job

obs_root = os.path.join(args.source, str(args.observation))
bfull_output = os.path.join(args.buffer_output, 'fullrate', obs_root)
sfull_output = os.path.join(args.storage_output, 'fullrate', obs_root)
bdown_output = os.path.join(args.buffer_output, 'downsampled', obs_root)
sdown_output = os.path.join(args.storage_output, 'downsampled', obs_root)

# construct the command to run this script
if args.mode == 'run':
    script = scripts[args.job_name]
    opts = vars(args)
    opts['keep_q'] = '--keep-q' if args.keep_q else ''
    opts['calframe'] = '--calframe {}'.format(args.calframe) if args.calframe else ''
    opts['downsample_only'] = '--downsample-only' if args.downsample_only else ''
    opts['storage_error'] = '--storage-error' if args.storage_error else ''
    opts['user_obs'] = '--user-obs' if args.user_obs else ''
    cmd = shlex.split(script['cmd'].format(**opts))

def run_cmd(cmd):
    # skip jobs as requested
    if args.downsample_only and args.job_name in ['scanify', 'store_full', 'check_full']:
        print('Fullrate job {} already complete, exiting'.format(args.job_name))
        raise SystemExit

    if args.storage_error and args.job_name in ['scanify']:
        print('Fullrate job {} already complete, exiting'.format(args.job_name))
        raise SystemExit

    try:
        sp.check_call(cmd)
    except sp.CalledProcessError as e:
        raise SystemExit('Error running job {} on observation {}'.format(
            args.job_name, obs_root))

# check return status
if args.mode == 'post':
    success = not args.job_pre_return and not args.job_return

# pre- and post-processing for each job

if args.job_name == 'scanify':
    if args.mode == 'pre':
        pass

    elif args.mode == 'run':
        run_cmd(cmd)

    elif args.mode == 'post' and success:
        pass

elif args.job_name == 'downsample':

    if args.mode == 'pre':
        # create output directory
        if not os.path.exists(bdown_output):
            os.makedirs(bdown_output)

        # copy calibration file
        in_cal = os.path.join(sfull_output, 'nominal_online_cal.g3')
        out_cal = os.path.join(bdown_output, 'nominal_online_cal.g3')
        shutil.copy2(in_cal, out_cal)

        if os.stat(in_cal).st_size != os.stat(out_cal).st_size:
            raise OSError('Failed to copy calibration file {}'.format(in_cal))

    elif args.mode == 'run':
        run_cmd(cmd)

    elif args.mode == 'post' and success:
        pass


elif args.job_name in ['store_full', 'store_down']:

    if args.mode == 'pre':
        pass

    elif args.mode == 'run':
        run_cmd(cmd)

    elif args.mode == 'post' and success:
        pass

elif args.job_name in ['check_full', 'check_down']:

    if args.mode == 'pre':
        pass

    if args.mode == 'run':
        if args.downsample_only and args.job_name == 'check_full':
            print('Fullrate job {} already complete, exiting'.format(args.job_name))
            raise SystemExit

        try:
            out = sp.check_output(cmd, stderr=sp.STDOUT).decode()
        except sp.CalledProcessError as e:

            output = e.output.decode() if e.output else ""
            if 'Scan frame separation' not in output:
                raise SystemExit('Observation {} failed verification:\n{}'.format(
                    obs_root, output))

    elif args.mode == 'post':
        pass


# raise an error if requested
if args.mode == 'post':

    if not success:
        raise SystemExit('Error running job {} on observation {}'.format(args.job_name, obs_root))

# XXX TODO
# modify the transfer_manager to treat the ScanifyDatabase as read-only
# parse TransferFileDatabase here and update ScanifyDatabase as necessary
# move verification logic here
