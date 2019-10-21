#!/bin/bash
# executable is the location of the executable
executable = /home/dyang18/SampleExecutable.sh

# universe is into which queue the job should be submitted into
universe = vanilla

# Location where log files should be put. The variable 
# `$(cluster)` is the unqiue ID assigned by HTCondor
# to the job. 
Error = /scratch/dyang18/osglogs/std_error_$(cluster)_log
Output  = /scratch/dyang18/osglogs/std_output_$(cluster)_log
Log     = /scratch/dyang18/osglogs/condor_job_$(cluster)_log

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
request_memory = 6144
# `request_disk` tells HTCondor how much 
# memory (in KB) are needed.
# Here 1 GB of disk are requested. 
request_disk = 4194304
# `transfer_input_files` defines which files should
# be transfered to the remote worker from the
# local machine. 
# NOTE: This has to be a comma separated list of files. 
transfer_input_files = /home/dyang18/user.cert,/home/dyang18/spt3g/spt3g_software/calibration/scripts/fluxpointcal/makecoadd.py
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
arguments = 0000.g3 output.g3 offline_calibration.g3
# Defines the end of the jb submission file and tells HTCondor
# to queue a job. The number, in this example, 1
# signifies how many copies of the job should be submitted.
# For example, the line `queue 5` will submit 5 copies
# of the job.
queue 1
