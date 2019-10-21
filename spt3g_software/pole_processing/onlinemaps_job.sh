#!/bin/sh

set -e

export HOME=/home/sptdaq
unset PYTHONPATH
eval `/software/clustertools/py3-v3/setup.sh`
export SPT3G_SOFTWARE_BUILD_PATH=/software/spt3g_software_onlinemaps
export PYTHONPATH=$SPT3G_SOFTWARE_BUILD_PATH:$PYTHONPATH
export LD_LIBRARY_PATH=$SPT3G_SOFTWARE_BUILD_PATH/spt3g:$LD_LIBRARY_PATH

BAND=$1
OUTPUT=$2
shift 2
INPUTS=$@

# copy inputs to scratch
cd $_CONDOR_SCRATCH_DIR
mkdir data
LOCAL_INPUTS=""
for f in $INPUTS; do
    rsync $f data/.
    LOCAL_INPUTS="$LOCAL_INPUTS data/`basename $f`"
done

# create output directory (for both polarized and unpolarized maps)
mkdir output
LOCAL_OUTPUT="output/`basename $OUTPUT`"
REMOTE_OUTPUT="`dirname $OUTPUT`"

# run the job
python -u $SPT3G_SOFTWARE_BUILD_PATH/spt3g/std_processing/mapmakers/master_field_mapmaker.py $LOCAL_INPUTS -o $LOCAL_OUTPUT -z --bands-to-use $BAND --config-file $SPT3G_SOFTWARE_BUILD_PATH/spt3g/std_processing/mapmakers/pole_autoproc_config.yaml

# copy output to storage
rsync -aviP output/* $REMOTE_OUTPUT
