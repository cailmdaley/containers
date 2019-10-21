#!/bin/bash
##########################################################################
# A happy little script that records memory use of the DAQ
# NDH
##########################################################################

INTERVAL=10 # time, in s, between records
PROC_STRING="DAQSCRIPT" # record memory use for a process with this name
LOGFILE=/daq_logs/mem_use_`date -u +%Y%m%d_%H%M%S`.log
rm /daq_logs/mem_use.log 2>/dev/null
ln -s $LOGFILE /daq_logs/mem_use.log
echo '               ' `ps v|head -n1` >> $LOGFILE # just get the column names
while true; do
    echo `date -u +%Y%m%d_%H%M%S` `ps av|grep $PROC_STRING|grep -v grep` >> $LOGFILE
    sleep $INTERVAL
done
