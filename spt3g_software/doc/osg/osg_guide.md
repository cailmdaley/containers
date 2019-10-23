# UChicago Computing User Guide for SPT Collaboration

## Open Science Grid

The Open Science Grid (OSG) consists of computing and storage elements at over 100 individual sites spanning the United States. These sites are primarily at universities and national labs and range in size from a few hundred to tens of thousands of CPU cores. The distributed nature of these resource providers allows users submit their jobs at a single entry point and have them execute at whatever resource is available. Sharing is a core principle of the OSG. Over 100 million CPU hours used on the OSG in the past year were utilized opportunistically, i.e. wherever there were resources available that would otherwise have remained idle. This aspect of the OSG is what allows individual researchers gain access to large computing resources.

## Available Machines and Storage

The resources available to SPT are split between compute and storage. SPT has two dedicated machines, the "login" nodes, hosted at the ATLAS Midwest Tier 2 at UChicago (MWT2) that are for general purpose computing and access to MWT2/OSG resources. For storage, there is currently about ~50 TB of storage available on our OSG Connect Stash Ceph pool. Separately, MWT2 has setup a VM which will act as a data gateway from Denver to the storage and to NERSC and/or ANL.

![alt text](https://github.com/SouthPoleTelescope/spt3g_software/blob/master/doc/osg/Processing_flows.png?raw=true "SPT UChicago Infrastructure")

### Login Nodes

`amundsen.grid.uchicago.edu` and `scott.grid.uchicago.edu` are the two SPT login nodes available through MWT2. They are meant as a gateway to submit jobs to MWT2/OSG, for general purpose data analysis for the collaboration, and access point for SPT 3G data stored in the Ceph pool. The access to MWT2/OSG resources is handled through a dedicated HTCondor queue, more details below, that is available on both machines that target MWT2/OSG. For general purpose data analysis tasks, the machines can be used to create plots, etc. The machines have a [JuypterHub](http://jupyter.org/) server running at https://scott.grid.uchicago.edu:8000 and https://amundsen.grid.uchicago.edu:8000, respectively,  that will spawn Python2 or Python3 Jupyter notebooks with the SPT 3G software dependencies loaded by default. Both machines have the SPT storage area in OSG Connect Stash mounted at `/spt`. More details on the data organization are available below. 

### Ceph

The OSG Connect Stash storage pool is based on the [Ceph distributed file system](http://ceph.com/). It has a total of 2.8 PB total raw storage. For SPT, there is an `/spt` area that has ~50 TB of storage available at this point. It is three times replicated across the area for data safety. 

The storage area is organized as follows:

* `/spt/data` - Space for the 3G data coming from the satellite. Only readable by normal user.
* `/spt/user` - User storage area meant for large intermediate analysis outputs, or final analysis outputs before they get moved. Every user has their own directory. 
* `/spt/public/` - Small storage area meant items that need to be http access, i.e. code tarballs. Every user has their own directory. 
* `/spt/analysis` - Meant as a storage area for final level data of analyses. Only readable by normal users. Write access upon request and approval. 
* `/spt/simulation` -Storage area for common simulation outputs. Only readable by normal user.


### Buffer

Separate from the login nodes and the storage pool, we have setup a storage buffer machine at `spt-buffer.grid.uchicago.edu`. This machine is meant as a data ingest machine that will grab the data from USAP server's in Denver and distribute it to the various SPT-specific storage end points. The plan is for it to `SFTP` the data from Denver to the local 4 TB storage array. From there it will copied to the `/spt/data` storage area and then copied to the tape backup endpoints, i.e. NERSC and ANL. 

## SPT Connect

SPT Connect is part of the MWT2 OSG Connect infrastructure. We have set up an SPT top group inside our OSG Connect group hierachy to allow for better separation between SPT and other OSG Connect users.

## Getting an account

You can sign up for an account at [OSG Connect](http://www.osgconnect.net/). When signing up for an account please make sure that your username is your Globus ID (typically `<your_username>@globusid.org`), if your username does not end in `@globusid.org` we cannot create your account. Add yourself to the `spt` and `spt.all` group. Once your account has been approved, your account will be provisioned on `scott.grid.uchicago.edu` and `amundsen.grid.uchicago.edu` and a personal storage area at `/spt/user/<username>` and HTTP-accessible area at `/spt/public/<username>` will be created. The HTTP-accessible area will be created at `http://stash.ci-connect.net/spt/public/<username>/`

### Getting a grid certificate

To get a grid certificiate, please follow the instructions [here](https://opensciencegrid.org/docs/security/user-certs/) and [here](https://pole.uchicago.edu/spt3g/index.php/Grid_Certificate_How-To). Once you have the `user_certificate_and_key.*.p12` file in hand. Please add it to your machines authentication system, i.e. KeyChainAccess on macOS or Firefox Certificate Manager. Once you have added to your machines authentication system go to the [SPT VO website](https://spt-mgt.grid.uchicago.edu:8443/voms/SPT) and sign up for the SPT VO. Alternatively, you can sign up for the OSG VO on [their website](https://voms.opensciencegrid.org:8443/voms/osg).

## High-Throughput Computing and HTCondor

OSG Connect uses [High-Throughput Condor](http://research.cs.wisc.edu/htcondor/manual/v8.5/index.html) (HTCondor or historically Condor) as a job scheduler. Most campus clusters, like UChicago's RCC, or High-Performance Computing centers, like NERSC or ANL, use PBS or one of its derivatives (SLURM and LSF, for example). HTCondor, as the name suggests, is based on the [high-throughput computing (HTC) model](https://en.wikipedia.org/wiki/High-throughput_computing). This model is based on the assumption that the compute resources are independent of each other and distributed (ranging from multiple local servers to servers around the world), tasks requiring a single CPU core (can be scaled to single nodes in newer HTCondor versions) and "small" amounts of memory (preferablly ~2 GB/core, starting at 8 GB/core the number of available resources drops precipitously), and that every job is independent of other concurrently running jobs. By contrast, the [high-performance computing (HPC) model](https://en.wikipedia.org/wiki/High-performance_computing) used by supercomputers (for example NERSC's Cori or Edison, ORNL's Titan) and most campus clusters assumes that a given job requires an aggregate of multiple CPU cores and/or nodes working together in a tightly coupled network, i.e. requiring special networking technology such as Infiniband.

## Jobs on OSG

Getting jobs running with HTCondor and OSG can have a fairly step learning curve. The two main issues that you (the user) will face are the distribution of software and transfer of input and/or output data. 

### Software Distribution - CVMFS and Tarballs

To distribute the SPT 3G software dependencies, we have established an external OASIS [CVMFS](https://cernvm.cern.ch/portal/filesystem) repository for SPT. It is available at OSG sites under `/cvmfs/spt.opensciencegrid.org/`. To make changes to the repository, please contact Nathan. 

Distribution the SPT 3G software itself currently requires the use of software tarballs, i.e. generating `*.tar` (or `*.tar.gz`) file that contains the compiled binaries. This usually requires just running `make tarball` instead of a plain `make`. This tarball can then be put into your `/spt/public/<username>` directory and downloaded from within the job within using `wget`, i.e. `wget http://stash.ci-connect.net/spt/public/<username>/<name_of_tarball>.tgz`.

Please note: Workloads using proprietary software, i.e. MATLAB, IDL, etc., will have to be statically compiled (if you don't know what this means, let us know) to be used on OSG. We cannot distribute versions of proprietary software through CVMFS because of licensing issues, i.e. you can not legally run your MATLAB code using the UChicago license on a UC Berkley machine.

### Data Transfer

To transfer data in and out for the job, we recommend using [GridFTP](http://toolkit.globus.org/toolkit/docs/latest-stable/gridftp/). GridFTP is the underlying file transfer protocol for [Globus](https://www.globus.org/). Using GridFTP requires having a grid certificate, and the Uniform Resource Identifier (URI) of the input and/or output file. Instructions on how to obtain a grid certificate can be found [here](https://twiki.grid.iu.edu/bin/view/Documentation/CertificateGetWeb).

The grid certificate acts similar to a SSH key pair with an attached identity. It authenticates you with the server and provides an identity to map you to. To use the grid certificate you need to generate a user proxy. This is a public key that will expire after some time. To generate the user proxy, you need to run `grid-proxy-init`. By default this command will generate the proxy in `/tmp/x509up_u<unix_user_id>` and the expiration date of the proxy will be fairly short. Since don't want to have to renew the proxy every 12 hours, we can extend the lifetime to up to 168 hours (7 days), using the `-valid` option, i.e. `grid-proxy-init -valid 168:00`. Similarly, we do not want to dig around in `/tmp` to find our user proxy file. We can generate a proxy in a different location using the `-out` option. Putting it all together we get something like `grid-proxy-init -valid 168:00 -out my_proxy`. For more details on `grid-proxy-init`, see [here](http://toolkit.globus.org/toolkit/docs/4.0/security/prewsaa/rn01re05.html)

The URI of a file can be most readily described as the GridFTP/POSIX-like URL of the file. A sample file URI is `gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/spt/data/20170101/observation.g3`. It is split into three parts. The first portion defines the transfer protocol, i.e. `gsiftp://`. The second portion provides the server and port the file should be transferred from, i.e. `gridftp.grid.uchicago.edu:2811`. The final portion, i.e. `/cephfs/spt/data/20170101/observation.g3` is the absolute POSIX path of the file on the GridFTP server. To access a local file or define where the file should be put on the local file system the URI becomes, for example, `file:///path/to/file/on/local/system/`.

Putting all of this together some example GridFTP commands are:

To transfer files to the worker node from remote server, i.e. transfer input files:

`globus-url-copy -vb gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/spt/data/20170101/observation.g3 file:///path/to/local/work/input/directory/observation.g3`

To transfer files from worker to remote server, i.e. transfer output files:

`globus-url-copy -vb file:///path/to/local/work/output/directory/output_file.g3 gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/spt/user/foo/processed.g3`

`-vb` provides the transfer rate, useful for debugging sites.

For more detail documentation on `globus-url-copy` see [here](http://toolkit.globus.org/toolkit/docs/latest-stable/gridftp/user/). 

### Job Submission in HTCondor

Job submission in HTCondor is different than in PBS and its derivatives. In PBS, etc. one submits a bash-like script, which defines the job resource requirements and the commands that should be executed in a single script. HTCondor separates the job execution script or program from the submission file. A job submission file is a text file that defines the job requirements through the HTCondor job description language. In broad terms the submit file defines what the job requirements are, what the executable is, what the arguments to the executable are, which files should be transferred through the HTCondor file transfer mechanism (only useful for small files), how many of this type to queue, etc. 

A sample submit file is given below:

```
#!/bin/bash
# executable is the location of the executable
executable = /path/to/executable

# universe is into which queue the job should be submitted into
universe = vanilla

# Location where log files should be put. The variable 
# `$(cluster)` is the unqiue ID assigned by HTCondor
# to the job. 
Error = /path/to/log/files/std_error_$(cluster)_log
Output  = /path/to/log/files/std_output_$(cluster)_log
Log     = /path/to/log/files/condor_job_$(cluster)_log

# Job requirements 
# Example: `(GLIDEIN_ResourceName =!= "NPX")` means that the 
#          advertised resource name cannot match "NPX"
# Example: `(OSGVO_OS_STRING == "RHEL 6")` 
#          means that the advertised `OpSysAndVer` 
#          (or operation system and version) is RHEL6 or one of this derivaties
# Example: `HAS_CVMFS_spt_opensciencegrid_org`
#          means that we require the SPT cvmfs to be present
Requirements = ( (OSGVO_OS_STRING == "RHEL 6") && (GLIDEIN_ResourceName =!= "NPX") && (HAS_CVMFS_spt_opensciencegrid_org))

# Special condition to get access to the 
# Friends of MWT2 queue available from
# login.ci-connect.uchicago.edu
+WANT_RCC_ciconnect = True
# `request_cpus` tells HTCondor how many CPU cores
# this jobs needs, default value is 1
request_cpus = 1
# `request_memory` tells HTCondor how much 
# memory (in MB) are needed. 
# Here 2 GB are requested, which is the
# default for OSG. 
request_memory = 2048 
# `request_disk` tells HTCondor how much 
# memory (in KB) are needed.
# Here 1 GB of disk are requested. 
request_disk = 1048576
# `transfer_input_files` defines which files should
# be transfered to the remote worker from the
# local machine. 
# NOTE: This has to be a comma separated list of files. 
transfer_input_files = /path/to/user_cert,/path/to/code
# `transfer_output_files` similar to `transfer_input_files`
# this allows one to define which files should be transfered
# back to the local host
transfer_output_files = ""
# `when_to_transfer_output` defines when the output should be transferred
when_to_transfer_output = ON_EXIT
# Whether or not to transfer executable to remote location
transfer_executable = True
# This is for accounting purposes and allow you to get matched to certain
# private SPT slots with large memory
+ProjectName = "spt-all"
# Defines possible arguments to the executable defined above
# Here we just pass the input and output file to the executable 
# NOTE: Format for the arguments depends on your executable
arguments = <input_file> <output_file>

# Defines the end of the jb submission file and tells HTCondor
# to queue a job. The number, in this example, 1
# signifies how many copies of the job should be submitted.
# For example, the line `queue 5` will submit 5 copies
# of the job.
queue 1
```

The `executable` can be any piece of code that can be run in the shell, for example, compiled binary, Python script, or bash script. For OSG jobs, it is recommended to use a bash script as they typically require the transfer of input or output data, and setting up an environment.

To submit a your submit file to the HTCondor file simple execute

`condor_submit /path/to/your/submit/file`

HTCondor is very flexible in what requirements you can make of a job. You can, for example, have a job into a "Held" state instead being removed from the queue, have a job retried after failure, or have a site automatically blacklisted. A more advanced submit file that incorporates all of the three example is given below:

```
#!/bin/bash
# executable is the location of the executable
executable = /path/to/executable

# universe is into which queue the job should be submitted into
universe = vanilla

# Location where log files should be put. The variable 
# `$(cluster)` is the unqiue ID assigned by HTCondor
# to the job. 
Error = /path/to/log/files/std_error_$(cluster)_log
Output  = /path/to/log/files/std_output_$(cluster)_log
Log     = /path/to/log/files/condor_job_$(cluster)_log

# Job requirements
# Example: `HAS_CVMFS_spt_opensciencegrid_org`
#          means that we require the SPT cvmfs to be present
# Example: `(OSGVO_OS_STRING == "RHEL 6")` 
#          means that the advertised `OSGVO_OS_STRING` 
#          (or operation system and version) is RHEL6 or one of this derivaties
# Example: `RCC_Factory == "ciconnect"`
#          means that we overwrite the requirements if the job goes through the 
#          MWT2 friends queue
# Example: `TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName1`
#          means that the job will not go site it just failed on again. It is repeated
#          for every retry attempt. We are going to retry 5 times, so we need this 5 times.
Requirements = ( ( HAS_CVMFS_spt_opensciencegrid_org ) && ( ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName1 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName2 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName3 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName4 ) || ( RCC_Factory == "ciconnect" ) ) ) && ( OSGVO_OS_STRING == "RHEL 6" || ( RCC_Factory == "ciconnect" ) ) )

# Special condition to get access to the 
# Friends of MWT2 queue available from
# login.ci-connect.uchicago.edu
+WANT_RCC_ciconnect = True
# `request_cpus` tells HTCondor how many CPU cores
# this jobs needs, default value is 1
request_cpus = 1
# `request_memory` tells HTCondor how much 
# memory (in MB) are needed. 
# Here 2 GB are requested, which is the
# default for OSG. 
request_memory = 2GB 
# `request_disk` tells HTCondor how much 
# memory (in KB) are needed.
# Here 1 GB of disk are requested. 
request_disk = 35GB
# `transfer_input_files` defines which files should
# be transfered to the remote worker from the
# local machine
transfer_input_files = /path/to/user_cert,/path/to/code
# `transfer_output_files` similar to `transfer_input_files`
# this allows one to define which files should be transfered
# back to the local host
transfer_output_files = ""
# `when_to_transfer_output` defines when the output should be transferred
when_to_transfer_output = ON_EXIT
# Whether or not to transfer executable to remote location
transfer_executable = True
# This is for accounting purposes and allow you to get matched to certain
# private SPT slots with large memory
+ProjectName = "spt-all"
# Defines possible arguments to the executable defined above
# Here we just pass the input and output file to the executable 
# NOTE: Format for the arguments depends on your executable
arguments = <input_file> <output_file>

# This will automatically find your X509 user proxy
# and transfer it to the job, and set the environment accordingly
use_x509userproxy = True

# Retrying the job if it fails
# Send the job to Held state on failure. 
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)  

# Periodically retry the jobs for 5 times with an interval 1 hour.
# `NumJobStarts`: Number of times the job has entered the "Running" state
# `CurrentTime`: Current timestamp
# `EnteredCurrentStatus`: Timestamp when the job was put in the current state, in this case "Held".
periodic_release =  (NumJobStarts < 5) && ((CurrentTime - EnteredCurrentStatus) > (60*60))

# Defines the end of the jb submission file and tells HTCondor
# to queue a job. The number, in this example, 1
# signifies how many copies of the job should be submitted.
# For example, the line `queue 5` will submit 5 copies
# of the job.
queue 1
```

### Sample Executable

In general, a executable for an OSG consists of four main parts:

1. Setting up the code environment
2. Transfering input data and, if necessary, executable with dependencies, to the worker node
3. Running code
4. Transfering output data

A generic sample executable using CVMFS would be:

```
#!/bin/bash

# Takes 2 arguments:
# 1 - Input file URI
# 2 - Output file URI

# Printing hostname and environment
# Good to have for debugging, in case a site causes issues
# and needs to blacklisted or to test whether to 
# whitelist a site again. 
echo $HOSTNAME
env

# pseudo code
source environment

start_dir=$PWD
if [ "${OSG_WN_TMP}" == "" ];
then
    OSG_WN_TMP=$PWD
fi

mkdir $PWD/tmp/
work_dir=$PWD/tmp/
cd ${work_dir}
mkdir ${work_dir}/output

globus-url-copy -vb $1 file://${work_dir}

${start_dir}/script.sh ${work_dir}/input.file ${work_dir}/output/output.file

globus-url-copy -vb file://${work_dir}/output/output.file $2
```

A more SPT specific example:

```
#!/bin/bash

# Takes 2 arguments:
# 1 - Input file URI
# 2 - Output file URI

# Printing hostname and environment
# Good to have for debugging, in case a site causes issues
# and needs to blacklisted or to test whether to 
# whitelist a site again. 
echo $HOSTNAME
env

# Just setting things up so the glidein can 
# clean up after us
start_dir=$PWD
if [ "${OSG_WN_TMP}" == "" ];
then
    OSG_WN_TMP=$PWD
fi

# Get SPT3G code
eval `/cvmfs/spt.opensciencegrid.org/py2-v1/setup.sh`
mkdir $PWD/code/
wget http://stash.osg-connect.net/spt/public/<username>/tarball.tgz
tar xzf tarball.tgz -C $PWD/code

# Setting up 3G Software environment
$PWD/code/env-shell.sh

# Making a new workspace so we don't
# accidentally delete our scripts
mkdir $PWD/tmp/
work_dir=$PWD/tmp/
cd ${work_dir}
mkdir ${work_dir}/output

# Transferring input file
# `$1` is the cmd line argument that is the fully qualified
# URI of the file
# `file://${work_dir}/` is where we are putting the file
globus-url-copy -vb $1 file://${work_dir}/

# Running our processing script
# `input.file` is the basename of our input file, i.e. if the input file is
# gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/spt/data/bolodata/downsampled/calibrator/3009897/0000.g3
# the basename is 0000.g3
${start_dir}/script.sh ${work_dir}/$(basename $1) ${work_dir}/output/output.file

# Transferring output file back to storage
globus-url-copy -vb file://${work_dir}/output/output.file $2

# Being nice and cleaning up after outselves
rm -rf $work_dir
```

### Job Monitoring in HTCondor

There are several ways to monitor a job in HTCondor. The most general is the HTCondor `qstat` replacement: `condor_q`. `condor_q` provides more options than `qstat`. To check only your jobs, either running, queued, or held, one can use:

`condor_q <optional username> <optional job id>`

To see only running jobs:

`condor_q -r <optional username> <optional job id>`

To see only idle jobs:

`condor_q -constraint 'JobStatus =?= 1 && Owner =?= "<username>"'`

To see held jobs:

`condor_q -constraint 'JobStatus =?= 5 && Owner =?= "<username>"'`

For a list of `JobStatus` codes see, [here](http://pages.cs.wisc.edu/~adesmet/status.html). There are lots more options using the `-constraint` option to limit the number of jobs that are displayed, such as "all jobs that have been running for x number of hours". 

NOTE: `condor_q` reports virtual memory not resident memory

#### What to do if my job has been in idle state for a long time?

Jobs will not run is most cases because the requirements are too restrictive, i.e. too much RAM, certain software requirements, etc. To figure out if this is the case you can use this command:

`condor_q <job id> -better-analyze -pool uct2-bosco.uchicago.edu:11012?sock=collector`

`condor_q <job id> -better-analyze -pool condor.grid.uchicago.edu`

`condor_q <job id> -better-analyze -pool osg-flock.grid.iu.edu`

This will check your job requirements against the MWT2 Friends queue (`uct2-bosco.uchicago.edu:11012?sock=collector`), MWT2 CoreOS pool (`condor.grid.uchicago.edu`), and OSG pool (`osg-flock.grid.iu.edu`)

#### What to do if my job is held?

Your job may be held for several reasons. Most common reasons are that the submit file has an error in it or that the logs files cannot created due to permissions issues. To find out what the held reason is you simply need to run:

`condor_q -analyze <jod id>`

This will print out the reason why the job was held.

### Best Practices on OSG with HTCondor

The volatile, heterogeneous, and distributed nature of OSG resouces make it hard to run jobs with "non-standard" resource requirements, either in terms of number of CPUs, RAM, disk space, or run time. In general, we recommend the following:

* Input File size: < 10 GB
* Output File size: < 10 GB
* Run time: 10 min to 12 hours
* CPU: preferably 1 core, up to 8 core
* RAM: preferably 2 GB/core, up to 24 GB/core

### Friends of MWT2

We have added a special "queue" for groups that work with MWT2. This will provide about 15% of the MWT2 pool towards users in designted VOs, at the moment this is Xenon1T, SPT, and Veritas. 

### Large Processing Needs - DAG and DAGMan

In a typical HTCondor works flow one deals with with multiple 1000s, 10000s, and even 10000000000s of jobs that may have inter-dependencies. The management, monitoring, and bookkeeping of such a large number of jobs, especially when there are job inter-dependencies, can be a daunting task for any user. Besides this there is also a job scheduler (as with any piece of software) stability aspect to keep in mind. With the large user base of OSG Connect and CI Connect, and each user submitting jobs, submitting a large batch of jobs into the HTCondor queue can lead to stability issues with scheduling jobs or accepting new jobs. To relief the user from the job management, monitoring, and bookkeeping task as well as allow the user to submit large batches of jobs in more manageable increments, HTCondor includes [Directed Acyclic Graph Manager (DAGMan)](http://research.cs.wisc.edu/htcondor/manual/latest/2_10DAGMan_Applications.html). 

Setting up a Directed Acyclic Graph (DAG) in HTCondor is a straightforward task. The most basic DAG file consists of

```
JOB spt.0 spt.submit
JOB spt.1 spt.submit
JOB spt.2 spt.submit
JOB spt.3 spt.submit
...
```

In the above example, job is defined through the `JOB <unqiue job ID> /path/to/your/submit/file` sequence. The `<unique job ID>`, `spt.0`, `spt.1`, and `spt.2` in the above example, only has to be unqiue to the DAG itself and can be reused.

In a normal processing situation there is different input parameters for every job, for example every job has a different input and output file. For this case we want to change both the submit and DAG file. First we want to replace the previously hardcoded values with variable replacement, for example in our example submit file we had the line:

```
arguments = <input_file> <output_file>
```

This can be replaced with:

```
arguments = $(inputfile) $(outputfile) 
```

where the hardcoded values have been replaced with the bash-like variables `$(inputfile)` and `$(outputfile)`. In a DAG file, these input variables can are defined through the `VARS <unqiue job ID> <sequence of input variables>`. This changes the above example to:

```
JOB spt.0 spt.submit
VARS spt.0 inputfile="/path/to/input/file0.zip" outputfile="/path/to/output/file0.root"
JOB spt.1 spt.submit
VARS spt.1 inputfile="/path/to/input/file1.zip" outputfile="/path/to/output/file1.root"
JOB spt.2 spt.submit
VARS spt.2 inputfile="/path/to/input/file2.zip" outputfile="/path/to/output/file2.zip"
JOB spt.3 spt.submit
VARS spt.3 inputfile="/path/to/input/file3.zip" outputfile="/path/to/output/file3.zip"
...
```

Note: The `""` around the variable values are required, even if assigning numerical values.



The heterogeneous nature of OSG resources may cause a job maybe evicted randomly without the job actually failing or a job may fail at one site and not another. We recommend retrying a job multiple times before considering it failed. To add the retry mechanism to a DAG file, the line `RETRY <unqiue job ID> <number of retries>` is added, such that

```
JOB spt.0 spt.submit
VARS spt.0 inputfile="/path/to/input/file0.zip" outputfile="/path/to/output/file0.root"
RETRY spt.0 5
JOB spt.1 spt.submit
VARS spt.1 inputfile="/path/to/input/file1.zip" outputfile="/path/to/output/file1.root"
RETRY spt.1 5
JOB spt.2 spt.submit
VARS spt.2 inputfile="/path/to/input/file2.zip" outputfile="/path/to/output/file2.zip"
RETRY spt.2 5
JOB spt.3 spt.submit
VARS spt.3 inputfile="/path/to/input/file3.zip" outputfile="/path/to/output/file3.zip"
RETRY spt.3 5
...
```

retries every job at least 5 times before considering it failed.

To define job inter-dependencies in a DAG, one needs to add the line `PARENT <list of unqiue job ID parents> CHILD <list of unique job IDs of children>`. For example, `spt.1` and `spt.2` depend on `spt.0` to finish, this would require to add `PARENT spt.0 CHILD spt.1 spt.2` to the DAG file, such that the DAG file looks like:

```
JOB spt.0 spt.submit
VARS spt.0 inputfile="/path/to/input/file0.zip" outputfile="/path/to/output/file0.root"
RETRY spt.0 5
JOB spt.1 spt.submit
VARS spt.1 inputfile="/path/to/input/file1.zip" outputfile="/path/to/output/file1.root"
RETRY spt.1 5
JOB spt.2 spt.submit
VARS spt.2 inputfile="/path/to/input/file2.zip" outputfile="/path/to/output/file2.zip"
RETRY spt.2 5
JOB spt.3 spt.submit
VARS spt.3 inputfile="/path/to/input/file3.zip" outputfile="/path/to/output/file3.zip"
RETRY spt.3 5
PARENT spt.0 CHILD spt.1 spt.2
...
```

Similarly we can say add the requirement that `spt.3` depends on `spt.1` and `spt.2` finishing successfully. Adding the line `PARENT spt.1 spt.2 CHILD spt.3` will accomplish this:

```
JOB spt.0 spt.submit
VARS spt.0 inputfile="/path/to/input/file0.zip" outputfile="/path/to/output/file0.root"
RETRY spt.0 5
JOB spt.1 spt.submit
VARS spt.1 inputfile="/path/to/input/file1.zip" outputfile="/path/to/output/file1.root"
RETRY spt.1 5
JOB spt.2 spt.submit
VARS spt.2 inputfile="/path/to/input/file2.zip" outputfile="/path/to/output/file2.zip"
RETRY spt.2 5
JOB spt.3 spt.submit
VARS spt.3 inputfile="/path/to/input/file3.zip" outputfile="/path/to/output/file3.zip"
RETRY spt.3 5
PARENT spt.0 CHILD spt.1 spt.2
PARENT spt.1 spt.2 CHILD spt.3
...
```

There are changes that you will need to do some work locally before or after the job can be submitted, for example moving files in the correct location for transfer to the job site or to combine output files of several DAGMan jobs at the end. For this case, DAGMan has the `SCRIPT`, `PRE`, and `POST` keywords for every job. For example, if `spt.0` needs to have a script run before executing and `spt.3` needs a script run after being executing, we add `SCRIPT PRE spt.0 pre_script.sh` and `SCRIPT POST spt.3 post_script.sh` to the DAGMan file. The new DAGMan will look something like:

```
JOB spt.0 spt.submit
VARS spt.0 inputfile="/path/to/input/file0.zip" outputfile="/path/to/output/file0.root"
RETRY spt.0 5
SCRIPT PRE spt.0 pre_script.sh
JOB spt.1 spt.submit
VARS spt.1 inputfile="/path/to/input/file1.zip" outputfile="/path/to/output/file1.root"
RETRY spt.1 5
JOB spt.2 spt.submit
VARS spt.2 inputfile="/path/to/input/file2.zip" outputfile="/path/to/output/file2.zip"
RETRY spt.2 5
JOB spt.3 spt.submit
VARS spt.3 inputfile="/path/to/input/file3.zip" outputfile="/path/to/output/file3.zip"
RETRY spt.3 5
SCRIPT POST spt.3 post_script.sh
PARENT spt.0 CHILD spt.1 spt.2
PARENT spt.1 spt.2 CHILD spt.3
...
```

The combination of `PRE`, `POST`, and `JOB` is called a "node" in DAGMan lingo. The success or failure of a node, i.e. if DAGMan considers a task to have succeeded, depends on the state of the last compontnent to have executed. In most cases, this will be the `POST` script. This means that even though the `PRE` or the `JOB` have failed, the `POST` script could return a success, thereby marking the "node" as successful. Before you start using `PRE`, and `POST` scripts review [the relevant section in the DAGMan documentation](http://research.cs.wisc.edu/htcondor/manual/latest/2_10DAGMan_Applications.html#SECTION003102400000000000000) and contact us so we can help with any issues you may encounter.

The order of the `JOB`, `VAR`, `RETRY`, `PRE`, `POST`, and `PARENT/CHILD` in the file does not matter. The ordering in the example is personal preference. For more complicated example for defining job inter-dependency, [see the documentation for DAGMan](http://research.cs.wisc.edu/htcondor/manual/latest/2_10DAGMan_Applications.html). 

To submit a DAGman to the HTCondor queue 

`condor_submit_dag /path/to/your/dag/file`

One can add additional configuration parameters to a DAGman submission either through the command line, use `condor_submit_dag -help` to see the options, or adding a config file. To add a config file, 

`condor_submit_dag -config /path/to/dagman/config/file /path/to/your/dag/file`

where an example DAGman config file is:

```
DAGMAN_MAX_JOBS_SUBMITTED=30000
DAGMAN_MAX_SUBMITS_PER_INTERVAL=10000
DAGMAN_USER_LOG_SCAN_INTERVAL=1
```

In the above example:

* `DAGMAN_MAX_JOBS_SUBMITTED`: Sets the maximum number of jobs this DAG can have in the queue, i.e. 30000. 
* `DAGMAN_MAX_SUBMITS_PER_INTERVAL`: Sets how many jobs can be submitted per interval, i.e. 10000.
* `DAGMAN_USER_LOG_SCAN_INTERVAL`: How often a job submission interval occurs in seconds, i.e. every second.


DAGMan's job monitoring abilities allow it to produce a "rescue" DAG. The "rescue" DAG contains all the jobs failed after the requested number of retries or the DAGman jobs that had not yet completed when the DAG was removed from the queue. The "rescue" DAG follows the naming pattern of `<dag_file_name>.rescueXXX` where "XXX" is simply a counter starting at 1 or in this case 001. To resubmit the jobs that failed one simply has to submit the DAG file again, i.e. `condor_submit_dag -config /path/to/dagman/config/file /path/to/your/dag/file`.

