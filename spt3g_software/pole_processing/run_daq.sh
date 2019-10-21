#!/bin/bash
# Start and stop bolometer DAQ.

source ~/.bash_profile

# Destination directory for tuning logs
path=/daq_logs

# Timestampe in UTC for logfile name. Determines the date upon execution of the
# script, i.e. corresponds to the start of DAQ. 
timestamp=`date -u +%Y%m%d-%H%M%S`
logfile=$path/daq-${timestamp}.log

# Copy STDOUT and STDERR of all subsequent command into a log file. 
# Copy original file descriptors to 3 and 4 in case we need them back later.
exec 3>&1 4>&2 >>$logfile 2>&1

# Link to most recent log
ln -s -f $logfile $path/daq.log

# Kill all other instances of this script that may be running. The script is
# renamed at exec, so this won't cause suicide.
oldscripts=`pgrep -f -x '^DAQSCRIPT .*'`
if [ -n "$oldscripts" ]; then
	echo Killing old DAQ scripts $oldscripts
	start=`date +%s`
	kill -INT $oldscripts

	# Wait for old processes to terminate
	while [ -n "`pgrep -f -x '^DAQSCRIPT .*'`" ]; do
		now=`date +%s`
		if [ $(($now - $start)) -gt 30 ]; then
			kill -9 $oldscripts
			start=$now
		fi
		sleep 1; 
	done
fi

if [ "$1" == "--stop" ]; then
	exit 0
fi

exec -a DAQSCRIPT python -u `dirname $0`/../dfmux/onlinescripts/record_bolodata.py $@
