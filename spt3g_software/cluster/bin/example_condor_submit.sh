# Bash script wrapper for the python command line script that submits jobs
# to the grid. This is a more convenient way than configuration files of
# keeping track of settings used for data processing jobs.

# submission script info
SUBMIT_SCRIPT=spt3g-condor-submit.py
WORKFLOW='process_obs'
SUBMIT=true
MEMORY=2147483648
DISK=32212254720

# `process_obs` arguments (use defaults, if empty)
JOB_SCRIPT=/home/adama/SPT/spt_analysis/20180831_mapmaking/field_autoproc.py
OBSTYPE=ra0hdec-44.75
OBSIDS='533*'
JOB_SCRIPT_ARGS="-s $OBSTYPE -r 2.0 -x 75 -y 50 -o output.g3"
LOG_ROOT=
OUTPUT_ROOT=

# setup derived variables and flags
if [ "$SUBMIT" = true ]; then
    SUBMIT_ARG="--submit-jobs"
else
    SUBMIT_ARG=""
fi

if [ -z "$LOG_ROOT" ]; then
    LOG_ROOT_FLAG=""
else
    LOG_ROOT_FLAG="--log-root $LOG_ROOT"
fi

if [ -z "$OUTPUT_ROOT" ]; then
    OUTPUT_ROOT_FLAG=""
else
    OUTPUT_ROOT_FLAG="--output-root $OUTPUT_ROOT"
fi

if [ -z "$MEMORY" ]; then
    MEMORY_FLAG=""
else
    MEMORY_FLAG="--request-memory $MEMORY"
fi

if [ -z "$DISK" ]; then
    DISK_FLAG=""
else
    DISK_FLAG="--request-disk $DISK"
fi


# because passing a list of arguments as an argument to a command line script
# generates parsing confusion unless the list is wrapped in quotes
JOB_SCRIPT_ARGS="'$JOB_SCRIPT_ARGS'"

# execute
cmd="python $SUBMIT_SCRIPT $SUBMIT_ARG $MEMORY_FLAG $DISK_FLAG $WORKFLOW $JOB_SCRIPT $OBSTYPE $OBSIDS --args $JOB_SCRIPT_ARGS $LOG_ROOT_FLAG $OUTPUT_ROOT_FLAG"
echo $cmd
eval $cmd
