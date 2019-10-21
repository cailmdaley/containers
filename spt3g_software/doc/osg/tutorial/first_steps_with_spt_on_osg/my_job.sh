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

start_dir=$PWD
if [ "${OSG_WN_TMP}" == "" ];
then
    OSG_WN_TMP=$PWD
fi

# Get SPT3G code
eval `/cvmfs/spt.opensciencegrid.org/py2-v1/setup.sh`
mkdir $PWD/code
# CHANGE THIS LINE BELOW
wget --no-cache http://stash.ci-connect.net/spt/public/<username>/spt3g.tgz 
tar xzf spt3g.tgz -C $PWD/code
# Invoke SPT 3G Software environment
$PWD/code/env-shell.sh 

# Making space to put data
mkdir $PWD/input/
mkdir $PWD/output/

# Getting input data
globus-url-copy -vb $1 file://$PWD/input/inputdata.g3

# Running code
python calculate_mean.py $PWD/input/inputdata.g3 $PWD/output/outputdata

# Pushing data back
globus-url-copy -vb file://$PWD/output/outputdata.g3 $2


