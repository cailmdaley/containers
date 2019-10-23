import os
import subprocess as sp
import shlex
import re
import glob
import fcntl
import numpy as np
import datetime
import shutil
import signal, errno
import pwd
import socket
from contextlib import contextmanager
from spt3g.std_processing.status_db import print_warn, print_time
import pexpect
import sys

class SFTPError(RuntimeError):
    """Error communicating with SFTP server"""

class RsyncError(RuntimeError):
    """Error Rsyncing files to/from remote"""

class RsyncConnectionError(RuntimeError):
    """Error connecting to rsync remote"""

class RemoteFullError(RuntimeError):
    """Remote disk is at capacity (user quota)"""

class InvalidFilenameError(RuntimeError):
    """Invalid filename"""

class LockedFileError(IOError):
    """File is already locked"""

class LockedProcessError(LockedFileError):
    """Process is already running"""

def isvalidname(fname):
    """
    Only allow transfer of valid filenames
    """
    if re.match("^[a-zA-Z0-9_\.\-\/\*]*$", fname):
        return True
    return False

def update_verify_root(verify_base, archive=False):
    """
    Create date-stamped verified data directory for backing up verified downloads
    """
    today = datetime.datetime.now().strftime('%Y%m%d')
    verify_root = os.path.join(verify_base, today)

    if not os.path.exists(verify_root):

        print_time('Creating new verified directory {}'.format(verify_root))
        os.makedirs(verify_root)

    if not archive:
        return verify_root

    cwd = os.getcwd()
    os.chdir(verify_base)

    for f in glob.glob('[0-9]'*8):

        if f >= today:
            continue

        if not len(os.listdir(f)):
            shutil.rmtree(f)
            continue

        if os.path.exists(f + '.tar'):
            # another process is already making the tarball
            continue

        try:
            print_time('Archiving verified directory {}'.format(f))
            sp.check_call(['tar', 'cf', f + '.tar', f])
        except sp.CalledProcessError as e:
            print_warn('Error archiving verified directory {}'.format(f), e)
            sp.check_call(['rm', '-f', f + '.tar'])
        else:
            print_time('Removing verified directory {}'.format(f))
            sp.check_call(['rm', '-rf', f])

    os.chdir(cwd)

    return verify_root

@contextmanager
def timeout_context(seconds):
    """
    A general-purpose timeout context for interrupting blocking processes.

    Borrowed from
    https://stackoverflow.com/questions/5255220/fcntl-flock-how-to-implement-a-timeout
    """
    if seconds is None:
        yield

    else:
        def timeout_handler(signum, frame):
            pass

        original_handler = signal.signal(signal.SIGALRM, timeout_handler)

        try:
            signal.alarm(seconds)
            yield
        finally:
            signal.alarm(0)
            signal.signal(signal.SIGALRM, original_handler)

def get_pids(pattern):
    try:
        cmd = 'pgrep -f "{}"'.format(pattern)
        pids = sp.check_output(shlex.split(cmd)).decode()
    except sp.CalledProcessError:
        pids = []
    else:
        pids = [int(x) for x in pids.split('\n') if x.strip()]
    return pids

def get_pname(pid=None):
    if pid is None:
        pid = os.getpid()
    try:
        cmd = 'ps -o command {}'.format(pid)
        out = sp.check_output(shlex.split(cmd)).decode().strip().split('\n')[1]
        pname = ' '.join(out.split()[:2])
    except sp.CalledProcessError:
        pname = None
    return pname

def send_cron_email(source, msg, address='cronjobs'):
    """
    Send an cron-like email

    Arguments
    ---------
    source : string
        Source script name for email subject
    msg : string
        Email body
    address : string
        Recipient email address
    """

    subject = 'Cron <{}@{}> {}'.format(
        pwd.getpwuid(os.getuid())[0], socket.gethostname(), source)
    p = sp.Popen(shlex.split('mail -s "{}" {}'.format(subject, address)),
                 stdin=sp.PIPE)
    p.communicate(msg.encode())

class LockedFile(object):

    def __init__(self, filename, mode='a+', buffering=None, lock=True, timeout=None):
        """
        LockedFile(name[, mode[, buffering[, lock[, timeout]]]]) -> file object

        Open and lock a file.  The mode can be 'r', 'w' or 'a' for reading
        (default), writing or appending.  The file will be created if it doesn't
        exist when opened for writing or appending; it will be truncated when
        opened for writing.  Add a 'b' to the mode for binary files.
        Add a '+' to the mode to allow simultaneous reading and writing.
        If the buffering argument is given, 0 means unbuffered, 1 means line
        buffered, and larger numbers specify the buffer size.  The preferred way
        to open a file is with the builtin open() function.
        Add a 'U' to mode to open the file for input with universal newline
        support.  Any line ending in the input file will be seen as a '\n'
        in Python.  Also, a file so opened gains the attribute 'newlines';
        the value for this attribute is one of None (no newline read yet),
        '\r', '\n', '\r\n' or a tuple containing all the newline types seen.

        'U' cannot be combined with 'w' or '+' mode.

        If lock is True (default), the file is locked using
        fcntl.lockf(f, LOCK_EX).  The lock argument can also be any valid value
        for the lockf operation argument.

        If the timeout is an integer (in seconds), then a LockedFileError is
        raised if the file lock cannot be established after the specified
        amount of time.
        """

        if buffering is None:
            self._file = open(filename, mode)
        else:
            self._file = open(filename, mode, buffering)

        self.locked = False
        if lock:
            if lock is True:
                lock = fcntl.LOCK_EX
            with timeout_context(timeout):
                try:
                    fcntl.lockf(self._file, lock)
                except IOError as e:
                    if e.errno != errno.EINTR and e.errno != errno.EDEADLK:
                        raise e
                    raise LockedFileError(
                        "File {} is locked by another process".format(filename))
                else:
                    # lock was successful
                    self.locked = lock

    def close(self):
        """
        close() -> None or (perhaps) an integer.  Unlock and close the file.

        Sets data attribute .locked to False.

        Sets data attribute .closed to True.  A closed file cannot be used for
        further I/O operations.  close() may be called more than once without
        error.  Some kinds of file objects (for example, opened by popen())
        may return an exit status upon closing.
        """
        if self.locked:
            fcntl.lockf(self._file, fcntl.LOCK_UN)
            self.locked = False
        self._file.close()
        self._file = None


class LockedProcess(LockedFile):

    def __init__(self, filename, *args, **kwargs):
        """
        {}

        Additionally, the LockedProcess object raises a LockedProcessError
        if a process using the same filename is already running, as determined
        by `pgrep -f`.
        """

        # this process
        this_pid = os.getpid()
        this_pname = get_pname(this_pid)

        # find all processes run on this file
        file_pids = get_pids(this_pname)
        other_pids = set(file_pids)

        # raise an error if any other such processes are already running
        if other_pids - set([this_pid]):
            raise LockedProcessError(
                "Process {} is already running".format(filename))

        super(LockedProcess, self).__init__(filename, *args, **kwargs)

    __init__.__doc__ = __init__.__doc__.format(LockedFile.__init__.__doc__)

class SFTPComm(object):

    def __init__(self, user=None, host=None, identity=None, verbose=True):

        # defaults
        if host is None and user is None and identity is None:
            user = 'A-379-S'
            if socket.gethostname().startswith('spt-buffer'):
                host = 'sftp.usap.gov'
                identity = '~/.ssh/A-379-S-2'
            else:
                host = 'spfs.southpole.usap.gov'
                identity = '~/.ssh/id_dsa_spfs'

        self.host = host
        if host is None:
            self.route = None
        elif user is None:
            self.route = host
        else:
            self.route = '{}@{}'.format(user, host)
        self.identity = identity
        self.sftp_cmd = 'sftp'
        self.prompt = 'sftp> '
        self.sftp_args = []
        if self.identity:
            self.sftp_args += ['-oIdentityFile={}'.format(self.identity)]
        self.sftp_args += [self.route]
        self.verbose = verbose
        self._proc = None

    def _expect(self):
        if self._proc is not None:
            try:
                self._proc.expect(self.prompt, timeout=1200)
            except pexpect.TIMEOUT as e:
                self._disconnect()
                raise SFTPError('Timed out waiting for prompt from SFTP server {}'.format(self.route))
            except pexpect.EOF as e:
                self._disconnect()
                raise SFTPError('Broken pipe from SFTP server {}'.format(self.route))

    def _connect(self):
        if self._proc is None:
            self.log('{}: Connecting to {} with identity {}'.format(
                self.__class__.__name__, self.route, self.identity))
            try:
                self._proc = pexpect.spawn(self.sftp_cmd, self.sftp_args, timeout=60)
                self._expect()
            except pexpect.TIMEOUT as e:
                self._disconnect()
                raise SFTPError('Timed out connecting to SFTP server {}'.format(self.route))
            except pexpect.EOF as e:
                self._disconnect()
                raise SFTPError('Broken pipe from SFTP server {}'.format(self.route))
        return self._proc

    def _communicate(self, cmd):
        proc = self._connect()
        try:
            proc.sendline(cmd)
            self._expect()
        except pexpect.TIMEOUT as e:
            self._disconnect()
            raise SFTPError('Timed out waiting for response from SFTP server {}'.format(self.route))
        except pexpect.EOF as e:
            self._disconnect()
            raise SFTPError('Broken pipe from SFTP server {}'.format(self.route))
        return proc.before.decode()

    def _disconnect(self):
        if self._proc is not None:
            self._proc.close()
            self._proc = None

    def __del__(self):
        self._disconnect()

    def log(self, msg):
        if self.verbose:
            print_time(msg)

    def __call__(self, cmd, parse=True):
        if not isinstance(cmd, list):
            cmd = shlex.split(cmd)
        cmd = ' '.join(cmd)
        out = self._communicate(cmd)
        # sanitize
        if parse:
            out = [x.strip() for x in ''.join(out.split('\n')).split('\r')][1:]
            out = '\n'.join([o for o in out if o and not o.startswith(self.prompt.strip())])
        return out

    def login(self):
        sp.check_call([self.sftp_cmd] + self.sftp_args)

    def upload_file(self, local, remote, min_upload_size=1024**3,
                    min_bundle_size=2 * 1024**3, priority=1, aux=False,
                    auxtype=None):
        """
        Return file size, upload status, and bundle filename.

        If file size is > min_upload_size, upload directly.
        Otherwise, populate a bundle until it is > min_bundle_size

        Bundle name is None if the file is uploaded directly,
        a tarball name if uploaded once large enough, 'current' otherwise.
        """
        local_sz = self.check_file_size_local(local)

        # if the file is large enough, send as is
        if local_sz >= min_upload_size:
            return self.send_file(local, remote), 'uploaded', None

        # otherwise move the file into the current bundle
        tag = self.local2remote('bundle', priority=priority, aux=aux).split('_')[0]
        bundle_dir = '/tmp/bundle_{}'.format(tag)
        if not os.path.exists(bundle_dir):
            os.makedirs(bundle_dir)
        sp.check_call(['mv', local, bundle_dir])

        # if the bundle is not large enough, return
        bundle_sz = np.sum(self.check_file_size_local(bundle_dir))
        if bundle_sz < min_bundle_size:
            return local_sz, 'uploading', 'current'

        # if the bundle is large enough, tar up and upload
        cwd = os.getcwd()
        os.chdir(bundle_dir)
        stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        tarfile = '/tmp/bundle_{}_{}.tar'.format(tag, stamp)
        files = os.listdir('.')
        sp.check_call(['tar', '-cvf', tarfile] + files)
        bundle_sz = self.send_file(tarfile, os.path.basename(tarfile))
        os.chdir(cwd)
        sp.check_call(['rm', '-rf', bundle_dir, tarfile])
        os.chdir(cwd)
        return local_sz, 'uploaded', os.path.basename(tarfile)

    def send_file(self, local, remote):
        self.log('{}: Sending {} to {}'.format(
            self.__class__.__name__, local, remote))
        if not isvalidname(local):
            raise InvalidFilenameError('Invalid local filename {}'.format(local))
        if not isvalidname(remote):
            raise InvalidFilenameError('Invalid remote filename {}'.format(remote))
        cmd = '"put {} {}"'.format(local, remote)
        try:
            ret = self(cmd)
        except Exception as e:
            if self.list_files(remote):
                self(['rm', remote])
            msg = 'File transfer of {} failed: {}'.format(local, str(e))
            raise SFTPError(msg)
        ret = '\n'.join([x for x in ret.split('\n') if not x.strip().endswith(' ETA')])
        self.log(ret)

        # compare file sizes
        local_sz = self.check_file_size_local(local)
        remote_sz = self.check_file_size_remote(remote)
        if np.sum(remote_sz) == 0 or np.sum(local_sz) == 0 or local_sz != remote_sz:
            msg = 'File transfer of {} failed.  Local size: {}.  Remote size: {}.'
            raise SFTPError(msg.format(local, local_sz, remote_sz))

        return remote_sz

    def download_file(self, remote, local, delete=True):
        """
        Return local filenames, file sizes, and parent bundle
        """
        sz = self.get_file(remote, local, delete=delete)

        # if the file is not a bundle, return
        if not remote.startswith('bundle_'):
            return [os.path.basename(local)], [sz], None

        output = os.path.dirname(local)
        files = self.unarchive(local, output=output, unpack=False,
                               delete=delete, verify=False)
        sizes = [self.check_file_size_local(os.path.join(output, f)) for f in files]
        return files, sizes, remote

    def get_file(self, remote, local, delete=True):
        self.log('{}: Receiving {} from {}'.format(
            self.__class__.__name__, local, remote))
        if not isvalidname(local):
            raise InvalidFilenameError('Invalid local filename {}'.format(local))
        if not isvalidname(remote):
            raise InvalidFilenameError('Invalid remote filename {}'.format(remote))
        cmd = '"get {} {}"'.format(remote, local)
        try:
            ret = self(cmd)
        except Exception as e:
            if 'Bad file descriptor' in str(e):
                raise SFTPError('File transfer of {} failed. {}'.format(remote, e))
            raise e
        ret = '\n'.join([x for x in ret.split('\n') if not x.strip().endswith(' ETA')])
        if 'No such file or directory' in ret:
            msg = 'File transfer of {} failed.  File not found.'
            self('"rename {} {}.corrupt"'.format(remote, remote))
            raise SFTPError(msg.format(remote))
        self.log(ret)

        # compare file sizes
        try:
            local_sz = self.check_file_size_local(local)
        except sp.CalledProcessError:
            local_sz = 0
        remote_sz = self.check_file_size_remote(remote)
        if np.sum(remote_sz) == 0 or np.sum(local_sz) == 0 or local_sz != remote_sz:
            msg = 'File transfer of {} failed.  Remote size: {}.  Local size: {}.'
            raise SFTPError(msg.format(remote, remote_sz, local_sz))

        if delete:
            self(['rm', remote])

        return local_sz

    def list_files(self, files, by_time=False, return_size=False):
        cmd = ['ls']
        if by_time:
            cmd += ['-tr']
        if return_size:
            cmd += ['-l']
        out = self(cmd).strip()
        # filter file list by search criteria
        if not isinstance(files, list):
            files = [files]
        patterns = [glob.fnmatch.translate(f) for f in files]
        out = '\n'.join([o for o in out.split('\n' if return_size else None)
                         if o.strip() and
                         any([re.match(p, o.split()[-1]) for p in patterns])])
        if not out:
            if return_size:
                return [], []
            return []
        if return_size:
            sz = np.atleast_1d(self._parse_size(out))
        out = [o.strip().split()[-1] for o in out.split('\n')]
        if return_size:
            return out, sz
        return out

    @staticmethod
    def _parse_size(llstr):
        llstr = llstr.strip()
        try:
            ret = np.asarray([int(v.split()[4]) for v in llstr.split('\n')
                              if len(v.split()) > 4])
        except ValueError:
            return 0
        if len(ret) == 1:
            return ret[0]
        if len(ret) == 0:
            return 0
        return ret

    @staticmethod
    def _parse_du(dustr):
        dustr = dustr.strip()
        ret = np.asarray([int(v.split()[0]) for v in dustr.split('\n')
                          if len(v.split())])
        if len(ret) == 1:
            return ret[0]
        if len(ret) == 0:
            return 0
        return ret

    def check_disk_usage_local(self, filename):
        return self._parse_du(sp.check_output(
            'du -bs "{}"'.format(filename), shell=True).decode())

    def check_file_size_local(self, filename):
        return self._parse_size(sp.check_output(
            'ls -l "{}"'.format(filename), shell=True).decode())

    def check_file_size_remote(self, filename):
        return self._parse_size(self(['ls', '-l', filename]))

    def remote_used(self):
        remote_size = np.sum(self.check_file_size_remote('*.tar'))
        remote_size += np.sum(self.check_file_size_remote('*.tar.*'))
        return remote_size

    def remote_available(self, max_size):
        remote_size = self.remote_used()
        available_size = max_size - remote_size
        if available_size <= 100 * 1024**2:
            raise RemoteFullError('SFTP directory is at capacity, exiting.')
        return available_size

    @staticmethod
    def local2remote(filename, priority=1, aux=False, auxtype=None):
        if aux:
            tag = 'aux'
            if auxtype:
                tag += '-{}'.format(auxtype)
        else:
            tag = 'p{}'.format(priority)
        if isinstance(filename, tuple):
            filename = os.path.join(*[str(x) for x in filename])
        return '{}_{}.tar'.format(tag, filename.replace(os.path.sep, '__'))

    @staticmethod
    def remote2local(filename, return_entry=False, return_tag=False):
        aux = False
        m = re.match('p([0-9]*)_(.*)', filename)
        if m:
            tag, filename = m.groups()
            tag = int(tag)
            nsep = 3
        else:
            m = re.match('aux-([^_]*)_(.*)', filename)
            if m:
                tag, filename = m.groups()
                aux = True
                if tag == 'pydfmux':
                    nsep = 1
                elif tag == 'pydfmuxlog':
                    nsep = 0
                else:
                    nsep = len(filename)
            else:
                m = re.match('aux_(.*)', filename)
                if m:
                    filename, = m.groups()
                    tag = None
                    aux = True
                else:
                    raise ValueError("Unrecognized filename {}".format(filename))
                nsep = len(filename)
        filename = filename.rsplit('.tar', 1)[0].replace('__', os.path.sep, nsep)
        if aux or not return_entry:
            if return_tag:
                return filename, tag
            return filename

        # parse filename into database entry format
        rate, src, obs, filename = filename.split(os.path.sep)
        obs = int(obs)
        if return_tag:
            return (rate, src, obs, filename), tag
        return (rate, src, obs, filename)

    def archive(self, path, output=None, checksum=False, compress=False, split=None,
                relative=False, priority=1, aux=False, auxtype=None, verbose=True):
        """
        Create file archive using tar

        Arguments
        ---------
        path : str
            filename, directory or list of files to archive
        output : str, optional
            output directory or tar filename
        checksum : bool, optional
            If True, include an md5 checksum of the input data with the archive
        compress: bool, optional
            If True, create an intermediate `.tar.gz` file of the input data
            to archive with the checksum
        split : int, optional
            If nonzero, divide the output tarfile into sections of the given size
            in MB.
        relative : bool, optional
            If True, ensure that input path(s) are relative to the current directory
        priority : int, optional
            Transfer priority
        aux : bool, optional
            If True, input path is an aux file
        auxtype : str, optional
            Type of aux file (for determining file paths)

        Returns
        -------
        tarfile : str or list of str
            Path to resulting tarfile or list of split tarfiles
        """

        if not isinstance(path, list):
            path = [path]

        # construct paths
        if output and not os.path.isdir(output) and \
                os.path.splitext(output)[-1] == '.tar':
            tarfile = output
        else:
            tarfile = self.local2remote(path[0], priority=priority, aux=aux,
                                        auxtype=auxtype)
            if output:
                tarfile = os.path.join(output, tarfile)
        if relative:
            cwd = os.getcwd()
            relpath = [os.path.relpath(x, '.') for x in path]
        else:
            relpath = path

        # archive and/or compress input path(s)
        if compress or (len(path) > 1) or any([os.path.isdir(x) for x in path]):
            tarfile1 = tarfile + '.gz' * compress
            self.log('Compressing {} to {}'.format(path, tarfile1))
            flags = '-c' + 'z' * compress + 'v' * verbose + 'f'
            cmd = ['tar', flags, tarfile1] + relpath
            sp.check_call(cmd)
            if relative:
                os.chdir(os.path.dirname(tarfile1))
                tarfile1 = os.path.basename(tarfile1)
            relpath = [tarfile1]

        # compute checksum
        if checksum:
            rpcs = []
            for rp in relpath:
                # use internal checksum for g3 files
                if rp.endswith('.g3') or rp.endswith('.g3.gz'):
                    continue
                rpc = rp + '.md5sum'
                if not os.path.exists(rpc):
                    self.log('Creating checksum for {}'.format(rp))
                    with open(rpc, 'w') as f:
                        sp.check_call(['md5sum', rp], stdout=f)
                rpcs += [rpc]
            # add to input file list
            relpath += rpcs

        # create output tarball
        self.log('Creating tar file {}'.format(tarfile))
        cmd = ['tar', '-cvf', tarfile] + relpath
        sp.check_call(cmd)

        # remove any intermediate files
        sp.check_call('rm -f {}.*'.format(tarfile), shell=True)

        # split into files of manageable size
        if split:
            self.log('Splitting tar file into chunks of size {} MB'.format(split))
            cmd = 'split -a 3 -d -b {bytes} {tarfile} {tarfile}.'.format(
                bytes=split * 1024**2, tarfile=tarfile).split()
            sp.check_call(cmd)
            # remove whole file
            sp.check_call(['rm', tarfile])
            tarfile = sorted(glob.glob('{}.*'.format(tarfile)))

        if relative:
            os.chdir(cwd)

        # return
        return tarfile

    def unarchive(self, tarfile, output=None, delete=True, unpack=True,
                  verify=True, delete_packed=True, return_checksums=False):
        """
        Extract from tar archive

        Arguments
        ---------
        tarfile : str
            filename or sequential list of filenames of the tar archive
        output : str, optional
            output directory into which the archive will be extracted
        delete : bool, optional
            If True, delete the tar file on successful extraction
        unpack : bool, optional
            If True, unarchive the internal tar(.gz) file on completion
        return_checksums : bool, optional
            If True, return checksum files along with unarchived data filenames

        Returns
        -------
        files : list of strings
            List of extracted files.  If `return_checksums` is True, this
            includes checksum files as well.
        """

        tarfiles = []
        # combine sequential tar files into one
        if isinstance(tarfile, list):
            tarfiles = tarfile
            m = re.search('(.*).tar', tarfiles[0])
            if m:
                tarfile, = m.groups()
            else:
                tarfile = tarfiles[0]
            tarfile += '.tar'

            with open(tarfile, 'wb') as f:
                sp.check_call(['cat'] + tarfiles, stdout=f)

        # tar file paths relative to output directory
        if output:
            os.chdir(output)

        # Extract tar file
        self.log('Extracting tar file {}'.format(tarfile))
        cmd = ['tar', '-xvf', tarfile]
        out = sp.check_output(cmd).decode()
        files = [x for x in out.split('\n') if x]

        # process checksums
        checksums = [f for f in files if f.endswith('.md5sum')] if verify else []
        for checksum in checksums:
            root = os.path.dirname(checksum)
            if root:
                sp.check_call(shlex.split(
                    'sed -i "s#/{}##g" {}'.format(root, checksum)))

            f = os.path.splitext(checksum)[0]
            try:
                sp.check_call(['md5sum', '--status', '-c', checksum])
            except sp.CalledProcessError as e:
                raise RuntimeError('File {} failed checksum:\n{}'.format(
                    f, e))
            else:
                self.log('File {} verified'.format(f))

        if unpack:
            # find all tar and non-tar files extracted from the archive
            subtarfiles = [f for f in files if f.endswith('.tar')
                           or f.endswith('.tar.gz')]
        else:
            subtarfiles = []
        files = [f for f in files if f not in checksums + subtarfiles]

        # decompress sub-tarfiles
        for f in subtarfiles:
            if os.path.dirname(f):
                shutil.move(f, os.path.basename(f))
                if os.path.exists(f + '.md5sum'):
                    shutil.move(f + '.md5sum', os.path.basename(f) + '.md5sum')
                f = os.path.basename(f)
            self.log('Decompressing {}'.format(f))
            compress = os.path.splitext(f)[-1] == '.gz'
            cmd = ['tar', '-x' + 'z' * compress + 'vf', f]
            out = sp.check_output(cmd).decode()
            files += [x for x in out.split('\n') if x]

        # checksum 3g files
        if verify:
            for f in files:
                if not f.endswith('.g3') or not f.endswith('.g3.gz'):
                    continue
                sp.check_call(['spt3g-verify', f])
                self.log('File {} verified'.format(f))

        # remove primary and intermediate tar files
        if delete:
            tarfiles += [tarfile]
            for f in tarfiles:
                if os.path.exists(f):
                    os.remove(f)
        if unpack and delete_packed:
            for f in subtarfiles:
                if os.path.exists(f):
                    os.remove(f)

        if return_checksums:
            return files + checksums
        return files

class RsyncComm(SFTPComm):

    def __init__(self, *args, **kwargs):
        super(RsyncComm, self).__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):

        # parse keywords
        stderr = kwargs.pop('stderr', None)
        return_err = kwargs.pop('return_err', None)
        compress = kwargs.pop('compress', True)
        local_cmd = kwargs.pop('local_cmd', None)
        if not local_cmd:
            local_cmd = 'rsync'
        remote_cmd = kwargs.pop('remote_cmd', None)
        if kwargs:
            raise KeyError('Unrecognized keyword argument(s) {}'.format(
                    kwargs.keys()))

        # exit if connection is broken
        if not self.isup():
            out = 'No route to host {}, connection unexpectedly closed'.format(
                self.host)
            if return_err:
                return out, ''
            return out

        if stderr is None:
            stderr = sp.PIPE
        cmd = shlex.split(local_cmd)
        cmd += ['-aviP' + 'z' * compress, '--checksum', '--timeout=100']
        if self.identity:
            cmd += ['--rsh="ssh -oIdentityFile={}"'.format(self.identity)]
        if remote_cmd:
            cmd += ['--rsync-path="{}"'.format(remote_cmd)]
        cmd += list(args)
        cmd = shlex.split(' '.join(cmd))
        p2 = sp.Popen(cmd, stdout=sp.PIPE, stderr=stderr)
        out, err = p2.communicate()
        # sanitize
        out = '\n'.join([o for o in out.decode().split('\n') if o])
        if return_err:
            return out, err.decode()
        return out

    def isup(self):
        ping_cmd = shlex.split('ping -c 1 -w 10 {}'.format(self.host))

        try:
            out = sp.check_output(ping_cmd, stderr=sp.PIPE).decode()
        except sp.CalledProcessError as e:
            return False

        m = re.search('([0-9]*) received', out)
        if not m:
            return False
        return bool(int(m.group(1)))

    def send_file(self, local, remote, lock=False, update=False, copy_dir_links=False,
                  copy_all_links=False):
        if self.route is not None:
            remote = '{}:{}'.format(self.route, remote)
        self.log('{}: Sending {} to {}'.format(
            self.__class__.__name__, local, remote))
        if lock:
            f = LockedFile(local)
            remote_cmd = 'flock -x {} rsync'.format(remote.split(':')[-1])
        else:
            remote_cmd = None

        local_cmd = 'rsync --copy-unsafe-links'
        if update:
            local_cmd += ' --update'
        if copy_dir_links:
            local_cmd += ' -k'
        if copy_all_links:
            local_cmd += ' -L'

        ret = self(local, remote, stderr=sp.STDOUT, remote_cmd=remote_cmd, local_cmd=local_cmd)
        if lock:
            f.close()

        m = re.search('sent ([0-9]*) bytes', ret)
        if not m:
            if 'connection unexpectedly closed' in ret or 'rsync error: timeout' in ret:
                self.log('{} Error: Timeout'.format(self.__class__.__name__))
                return 0
            raise RsyncError('Error sending {} to {}:\n{}'.format(
                local, remote, ret))

        # check file size
        b = int(m.group(1))
        self.log('{}: Sent {} bytes'.format(self.__class__.__name__, b))
        return b

    def get_file(self, remote, local, delete=False, lock=False):
        if self.route is not None:
            remote = '{}:{}'.format(self.route, remote)
        self.log('{}: Receiving {} from {}'.format(
            self.__class__.__name__, local, remote))
        if lock:
            f = LockedFile(local)
            remote_cmd = 'flock -x {} rsync'.format(remote.split(':')[-1])
        else:
            remote_cmd = None
        ret = self(remote, local, stderr=sp.STDOUT, remote_cmd=remote_cmd)
        if lock:
            f.close()

        m = re.search('received ([0-9]*) bytes', ret)
        if not m:
            if 'connection unexpectedly closed' in ret or 'rsync error: timeout' in ret:
                self.log('{} Error: Timeout'.format(self.__class__.__name__))
                return 0
            raise RsyncError('Error receiving {} from {}:\n{}'.format(
                local, remote, ret))

        # check file size
        b = int(m.group(1))
        self.log('{}: Received {} bytes'.format(self.__class__.__name__, b))
        return b
