# user options
SUBMIT_SCRIPT=/home/adama/SPT/spt3g_software/scratch/adama/condor/submit_condor_jobs.py
JOB_SCRIPT=/home/adama/SPT/spt_analysis/20180831_mapmaking/field_autoproc.py
OBSTYPE=ra0hdec-44.75
JOB_SCRIPT_ARGS='-s $OBSTYPE -r 2.0 -x 75 -y 50 -o output.g3'
OBSIDS='533*'
LOG_ROOT=/scratch/adama
SUBMIT=false

# setup derived variables and flags
if [ "$SUBMIT" = true ]; then
    SUBMIT_ARG='--submit-jobs'
else
    SUBMIT_ARG=''
fi

JOB_SCRIPT_ARGS='"$JOB_SCRIPT_ARGS"'

# execute
echo "python $SUBMIT_SCRIPT $JOB_SCRIPT $OBSTYPE $OBSIDS --args $JOB_SCRIPT_ARGS $SUBMIT_ARG"
python $SUBMIT_SCRIPT $JOB_SCRIPT $OBSTYPE $OBSIDS --args $JOB_SCRIPT_ARGS $SUBMIT_ARG
