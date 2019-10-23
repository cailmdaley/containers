import os
import shutil
import subprocess as sp
import glob
import time
import re
from spt3g import core

__all__ = ['condor_submit_baked', 'condor_submit', 'parse_dag_status',
           'parse_job_status', 'parse_error_log']

condor_submit_template = """
executable = {executable}
universe = vanilla

error = {log_root}/{jobname}_$(cluster).$(process).err
output = {log_root}/{jobname}_$(cluster).$(process).out
log = {log_root}/{jobname}_$(cluster).$(process).log
notification = never

Requirements = {requirements}
+WANT_RCC_ciconnect = True
+ProjectName = "spt.all"
+AccountingGroup = "group_opportunistic.spt.all.{user}"
run_as_owner = False

request_cpus = {request_cpus}
request_memory = {request_memory}
request_disk = {request_disk}

use_x509userproxy = {use_proxy}
x509userproxy = {grid_proxy}

should_transfer_files = YES
transfer_input_files = {transfer_input_files}
transfer_output_files = {transfer_output_files}
when_to_transfer_output = {when_to_transfer_output}
transfer_executable = True

{retries}

queue {queue}
"""

condor_script_template = """#!/usr/bin/env bash
# Operating Info
echo "Host: " `hostname`
echo "OS/Arch: " `uname -a`
echo ""

export OMP_NUM_THREADS={omp_threads}
export SHELL=sh
unset PYTHONPATH

env_on_error() {{
	echo 'Environment at error:'
	uname -a
	printenv
}}
trap env_on_error EXIT
set -e

echo "Setting up the environment"
eval `/cvmfs/spt.opensciencegrid.org/{clustertools_version}/setup.sh`

# Download software
{tarball_code}

# Get input files
{globus_input}

{user_code}

echo "Running user command"
ls
{caller} {script}

{user_code_post}

# Send back output
{globus_output}

trap - EXIT
echo "Job complete"
"""

def get_clustertools_version():
    sroot = os.getenv('SROOT')
    arch = os.getenv('OS_ARCH')
    if sroot is None or arch is None:
        raise OSError("Clustertools environment not loaded!")
    m = re.match('/cvmfs/spt.opensciencegrid.org/(.*)/{}'.format(arch), sroot)
    if not m:
        raise OSError("Invalid clustertools environment!")
    return m.group(1)


@core.usefulfunc
def condor_submit_baked(specs):
    """
    Submit an existing job specification and return the job ID.  See
    `condor_submit()` for a more user-friendly interface.

    Arguments
    ---------
    specs : string or list of strings
        A path to an existing condor job specification file, or a list of
        specifications (to pass to condor_submit via STDIN).

    Returns
    -------
    jobid : int
        The cluster ID to which the input job was submitted on success, None otherwise.
    result : string
        Output of the condor_submit call
    """
    if os.path.isfile(str(specs)):
        ret = sp.check_output(['condor_submit', specs])
    else:
        if isinstance(specs, list):
            specs = '\n'.join(specs)
        specs = specs.encode('ascii')
        try:
            ret = sp.check_output('condor_submit', input=specs)
        except TypeError:
            proc = sp.Popen(['condor_submit'], stdout=sp.PIPE, stdin=sp.PIPE)
            ret = proc.communicate(specs)[0]

    if not ret:
        return None, None

    # parse cluster ID
    ret = ret.decode('utf-8')
    m = re.search('submitted to cluster ([0-9]+).', ret)
    if not m:
        return None, ret

    return int(m.group(1)), ret


@core.usefulfunc
def condor_submit(script, args=[], caller='python', jobname=None,
                  log_root=None, output_root=None, input_files=[],
                  input_files_dest=None, output_files=[], grid_proxy=None, 
                  aux_input_files=[], aux_output_files=[], spt3g_env=True,
                  clustertools_version=None, omp_threads=1,
                  user_code='', user_code_post='', retry=True,
                  requirements=None, extra_requirements=None, bigmem=False,
                  request_cpus=1, request_memory=2*core.G3Units.GB,
                  request_disk=1*core.G3Units.GB, queue=1,
                  when_to_transfer_output='ON_EXIT', create_only=False,
                  globus_uri='gsiftp://gridftp.grid.uchicago.edu:2811/cephfs',
                  force_build_tarball=False, verbose=True):
    """
    Construct a condor submit file and run script for a command that requires
    the SPT3G software environment and data.

    Arguments
    ---------
    script : string
        User script to run.
    args : list of strings
        List of arguments to the script.
    caller : string
        Program that runs the script.
    jobname : string
        Name for the submitted job.  If not given, this is constructed
        using the current timestamp.
    log_root : string
        Local directory where the log files are stored.  This should be a
        location on /scratch, since these are typically lots of small files.
    output_root : string
        Local directory where the output files are stored.  This should be a
        location on cephfs (/spt), since output files must be transfered back
        using GridFTP.
    input_files : list of strings
        Files to transfer into the job using GridFTP.  Paths must be absolute
        paths on the cephfs file system (/spt).
    input_files_dest : list of strings
        Destination on grid for reach input file. Must be relative path. If
        None, put all input_files in working directory.
    output_files : list of strings
        Files created by the input `script` that are to be transfered back to the
        `output_root` using GridFTP.  Paths must be relative to the remote
        working directory where `script` is run.
    grid_proxy : string
        Path to a valid grid proxy file.  If None, the $X509_USER_PROXY
        environment variable is checked. Required if any `input_files` or
        `output_files` are supplied. 
    aux_input_files : list of strings
        Small files to transfer with the submit script.  If `script` is a path
        to an existing script, it is automatically added to this list.
        If `grid_proxy` is supplied, it is also automatically added to this
        list. Files must exist and will be transfered in to the remote working
        directory.
    aux_output_files : list of strings
        Small files to transfer to `log_root` on job completion.
        Paths must be relative to the remote working directory.
    spt3g_env : bool
        If True, a tarball of the `spt3g_software` build directory is
        created using `make tarball`, and transfered in to the running job
        using `wget`, and the `env-shell.sh` script is executed when running
        `script`.  The environment variable $SPT3G_BUILD_ROOT must be set
        in the environment in which this `condor_submit` function is called.
    force_build_tarball : bool
        If True and `spt3g_env` is True, rebuild the software tarball that
        is compiled and made available for download by the running job.
        NB: rebuilding the tarball will cause all previously queued condor
        jobs to use the new tarball.  It is up to the user to make sure this
        is intentional behavior.
    clustertools_version : string
        Version of the clustertools environment to use. This only understands
        'py2-v1', 'py3-v1', and 'py3-v2'.  If None (default), the version is
        parsed from the current user environment.
    omp_threads : int
        Number of OpenMP threads to use.  Default to a single thread
        to avoid accidentally grab unallocated resources.
    user_code : string
        If supplied, this text is added as-is to the run script, before
        the python script is called, but after the environment is configured
        and input files transfered.  The code is not sanitized.
    user_code_post : string
        If supplied, this text is added as-is to the run script, after
        the python script is called, but before output files are transfered.
        The code is not sanitized.
    retry: bool
        If True, when job does not complete succesfully, mark it as held and 
        periodically retry it up to 5 times. Not recommended if running with DAG.
    requirements : string
        Computing requirements specification.  If None, this defaults to the
        SPT3G default: 
        ((OSGVO_OS_STRING == "RHEL 7") && (GLIDEIN_ResourceName =!= "NPX")) &&
        (HAS_CVMFS_spt_opensciencegrid_org))
    extra_requirements : string
        Additional computing requirements to add to the defaults.  Use the
        `requirements` argument directly if additional requirements have
        more complicated logic than a simple && with the defaults.
    bigmem : bool
        Add the SPT_BIGMEM requirement to the requirements list.  This submits
        jobs directly to scott/amundsen.  The memory request should be at least
        16GB for a job to match this requirement.   The requirement is added as
        an || type requirement, so that jobs that do not match this can be
        sent out to the regular condor pool.
    request_cpus : int
        Number of CPUs to request for the job
    request_memory : integer
        Amount of memory required.  This should be supplied in units of
        core.G3Units.MB
    request_disk : integer
        Amount of disk space required.  This should be supplied in units of
        core.G3Units.MB
    queue : int or queue spec string
        The job queue specification.  Typically this is an integer
        specifying the number of copies of the job to submit, but can also
        be more complex -- see the Condor documentation for details.
    when_to_transfer_output : string
        When the files listed in `aux_output_files` and logs should be
        transfered back to the local `log_root`.
    create_only : bool
        If True, the submit and shell scripts for running the job are created
        but not submitted.  This mode is useful for debugging.
    globus_uri : string
        URI prefix to use for accessing files on the Ceph pool.

    Returns
    -------
    cluster : int
        If submission is successful, return the cluster ID number.
    """

    # file paths
    user = os.getenv('LOGNAME')
    if jobname is None:
        jobname = 'condor_{}'.format(time.strftime('%Y%m%d_%H%M%S'))
    if log_root is None:
        log_root = os.path.join('/scratch', user)
    log_root = os.path.abspath(log_root)
    if not os.path.exists(log_root):
        os.makedirs(log_root)
    if output_root is None:
        output_root = os.path.join('/spt/user', user)
    output_root = os.path.abspath(output_root)
    if not os.path.exists(output_root):
        os.makedirs(output_root)
    os.chmod(output_root, 0o777)

    executable = os.path.join(log_root, jobname+'.sh')
    if requirements is None:
        requirements = ('((OSGVO_OS_STRING == "RHEL 7") && ' +
                        '(GLIDEIN_ResourceName =!= "NPX") && ' +
                        '(HAS_CVMFS_spt_opensciencegrid_org){})'.format(
                            '&& ({})'.format(extra_requirements)
                            if extra_requirements else ''))
    if bigmem:
        requirements += ' || (SPT_BIGMEM =?= True)'
    if retry:
        retries = """
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)
periodic_release = (NumJobStarts < 6) && ((CurrentTime - EnteredCurrentStatus) >(10*60))
        """
    else:
        retries = ''
    script_absolute = os.path.abspath(script)
    if not os.path.exists(script_absolute):
        raise OSError('Cannot find user script {}'.format(script_absolute))
    aux_input_files = aux_input_files + [script]
    script = ' '.join([os.path.basename(script)] + args)

    # handle GridFTP transfers
    use_proxy = False
    globus_input = ''
    globus_output = ''
    if len(input_files) or len(output_files):
        if grid_proxy is None:
            grid_proxy = os.getenv('X509_USER_PROXY')
        if grid_proxy and not os.path.exists(grid_proxy):
            raise RuntimeError(
                'Missing grid proxy.  Please run voms-proxy-init and set '
                'the proxy location.')
        use_proxy = True
        if input_files_dest is not None:
            assert len(input_files) == len(input_files_dest)
        for i, input_file in enumerate(input_files):
            if not os.path.exists(input_file):
                raise IOError('Missing input file {}'.format(input_file))
            input_file = os.path.abspath(input_file)
            if input_files_dest is None:
                target_file = os.path.basename(input_file)
            else:
                target_file = input_files_dest[i]
            globus_input += 'globus-url-copy ' + \
                            '-rst -cd {}{} file://$PWD/{}\n'.format(
                    globus_uri, input_file, target_file)
                                
        for output_file in output_files:
            output_file_local = os.path.join(output_root,
                                             os.path.basename(output_file))
            globus_output += 'globus-url-copy -rst -stall-timeout 3600 file://$PWD/{} {}{}\n'.format(
                output_file, globus_uri, output_file_local)

    if grid_proxy is None:
        grid_proxy = ''

    transfer_input_files = ', '.join(
        [os.path.abspath(x) for x in aux_input_files])
    if not transfer_input_files:
        transfer_input_files = '" "'
    transfer_output_files = ', '.join(aux_output_files)
    if not transfer_output_files:
        transfer_output_files = '" "'

    request_memory = '{}MB'.format(int(request_memory / core.G3Units.MB))
    request_disk = '{}MB'.format(int(request_disk / core.G3Units.MB))

    if clustertools_version is None:
        clustertools_version = get_clustertools_version()

    tarball_code = ''
    if spt3g_env:
        build_root = os.getenv('SPT3G_BUILD_ROOT')
        if build_root is None:
            build_root = os.getenv('SPT3G_SOFTWARE_BUILD_PATH')
        if build_root is None:
            raise OSError("SPT3G software environment not loaded!")
        tarsource = os.path.join(build_root, 'spt3g_{}.tgz'.format(
            clustertools_version))
        if force_build_tarball or not os.path.exists(tarsource):
            if clustertools_version != get_clustertools_version():
                raise OSError("Clustertools version {} not loaded!".format(
                    clustertools_version))
            sp.check_call('make -C {} tarball'.format(build_root).split())
            sp.check_call('mv {0}/spt3g.tgz {0}/spt3g_{1}.tgz'.format(
                build_root, clustertools_version).split())
        tarfile = os.path.join('/spt/public', user,
                               'spt3g_{}.tgz'.format(clustertools_version))
        tarurl = 'http://stash.ci-connect.net{}'.format(tarfile)
        if force_build_tarball or not os.path.exists(tarfile):
            shutil.copy2(tarsource, tarfile)
        tarball_code = """
echo "Downloading software from {0}"
mkdir code
wget -O code/spt3g.tgz {0} --no-cache
cd code
tar xzf spt3g.tgz
cd ..
""".format(tarurl)
        caller = '$PWD/code/env-shell.sh {}'.format(caller)

    opts = locals()

    condor_submit = condor_submit_template.format(**opts)
    condor_script = condor_script_template.format(**opts)

    condor_submit_file = os.path.join(log_root, jobname+'.submit')
    with open(condor_submit_file, 'w') as f:
        f.write(condor_submit)

    condor_script_file = os.path.join(log_root, jobname+'.sh')
    with open(condor_script_file, 'w') as f:
        f.write(condor_script)
    os.chmod(condor_script_file, 0o755)

    if verbose:
        print('{} job {} for processing {}'.format(
            'Creating' if create_only else 'Submitting',
            condor_submit_file, condor_script_file))

    if create_only:
        return None, condor_submit_file, condor_script_file

    jobid, result = condor_submit_baked(condor_submit_file)
    print(result)
    if jobid:
        return jobid, condor_submit_file, condor_script_file


@core.usefulfunc
def parse_dag_status(filename):
    """
    Parse a DAG status file to determine state of a DAG job and each of its
    nodes.  Returns a dictionary of status records for each node.

    The output dictionary looks like this::

        status = {
            'DagStatus': {
                'DagStatus': 'done',
                ...
            },
            'NodeStatus': {
                'mynode1': {
                    'NodeStatus': 'done',
                    ...
                },
                'mynode2': {
                    'NodeStatus': 'done',
                    ...
                },
            }
        }

    Status values can be: 'not_ready', 'ready', 'prerun', 'submitted',
    'postrun', 'done' or 'error'.

    To generate a status file, add the following line to your DAG script::

        NODE_STATUS_FILE <filename>

    This file will be periodically updated as the job runs, and can be monitored
    using this parser function to check the status of a running job.

    For more details, see section 2.10.13 of the HTCondor DAGMan Applications
    manual:

    http://research.cs.wisc.edu/htcondor/manual/latest/DAGManApplications.html#x22-1100002.10.13
    """

    data = open(filename, 'r').read().strip()

    # remove all comments
    data = re.sub(re.compile("/\*.*?\*/", re.DOTALL), "", data)
    data = re.sub(re.compile("//.*?\n"), "", data)

    # remove new-line characters
    data = re.sub("\n", "", data)

    # replace semicolons commas for item separation
    data = re.sub(";", ", ", data)

    # convert ClassAd blocks into dictionaries
    data = re.sub("\[", "dict(", data)
    data = re.sub("\]", "), ", data)

    # convert ClassAd lists into python lists
    data = re.sub("\{", "[", data)
    data = re.sub("\}", "]", data)

    # evaluate data into a list of dictionaries
    # NB: this will fail for ClassAds with expressions rather than literals
    data = list(eval('[' + data + ']'))

    status_map = {
        0: 'not_ready',
        1: 'ready',
        2: 'prerun',
        3: 'submitted',
        4: 'postrun',
        5: 'done',
        6: 'error',
    }

    status = {}

    for record in data:
        tp = record.pop('Type')
        # convert status number to string description
        for k, v in record.items():
            if k in ['DagStatus', 'NodeStatus'] and v in status_map:
                record[k] = status_map[v]
        # group individual node statuses into a dict keyed by node name
        if tp == 'NodeStatus':
            if tp not in status:
                status[tp] = {}
            node = record.pop('Node')
            status[tp][node] = record
        else:
            status[tp] = record

    return status


@core.usefulfunc
def parse_job_status(logfile, job_id):
    """
    Parse a condor log file to determine the status of a job.  The log file should
    contain lines of the form::

        000 (357866.000.000) 03/17 13:35:38 Job submitted

    Arguments
    ---------
    logfile : string
        Path to the log file

    job_id : int
        ID number of the job to check

    Returns
    -------
    status : string
        Status of the job.  None if the status cannot be determined from the log file.
        Otherwise, one of: ['submitted', 'executing', 'stopped', 'failed', 'complete']
    statusline : string
        Latest entry in the log for the job ID.
    """
    with open(logfile, 'r') as f:
        log = f.readlines()
        statusline = None
        status = None
        for i in range(len(log)):
            # Parse log lines of the form: 000 (357866.000.000) 03/17 13:35:38 Job submitted
            line = log[i].split()
            if len(line) < 6:
                continue
            if line[0].isdigit() and line[4] == 'Job' and int(line[1].split('.')[0][1:]) == job_id:
                status = line[5].strip('.')
                if status == 'submitted':
                    pass
                elif status == 'executing':
                    pass
                elif status == 'disconnected,' or status == 'reconnection':
                    status = 'executing'
                elif status == 'was' and line[6] == 'held.':
                    status = 'stopped'
                elif status == 'was' and line[6] == 'released.':
                    status = 'submitted'
                elif status == 'was' and line[6] == 'evicted.':
                    status = 'stopped'
                elif status == 'was' and line[6] == 'aborted':
                    status = 'failed'
                elif status == 'terminated':
                    if 'Abnormal termination' in log[i+1]:
                        status = 'failed'
                    else:
                        exitcode = re.match('.*\(return value ([0-9]+)\)$', log[i+1])
                        assert exitcode is not None and len(exitcode.groups()) > 0, 'Failure in job %s in %s' % (job_id, logfile)
                        if exitcode.groups()[0] == '0':
                            status = 'complete'
                        else:
                            status = 'failed'
                else:
                    status = None
                statusline = log[i]
        return status, statusline


@core.usefulfunc
def parse_error_log(errlog, tail=20):
    """
    Return the last lines of the error log file that contain the error that
    caused the job to fail.  This includes all lines after a FATAL or ERROR log
    message, or the traceback from a python exception.

    If the log file cannot to be parsed correctly, the last `tail` lines of the
    log are returned instead.

    Arguments
    ---------
    errlog : string
        Path to error log file, or contents thereof.
    tail : int
        Number of log lines to return if the log cannot be parsed correctly.
        If 0, an empty string is returned.

    Returns
    -------
    string
    """

    if '\n' in errlog:
        lines = errlog.splitlines(True)
    else:
        lines = open(errlog, 'r').readlines()

    for idx, line in enumerate(lines):
        if line.startswith('ERROR') or line.startswith('FATAL') or 'Exception' in line \
           or line.startswith('Traceback') or 'failed' in line.lower() \
           or 'segmentation' in line.lower():
            return ''.join(lines[idx:])

    if not tail:
        return ''
    return ''.join(lines[-tail:])
