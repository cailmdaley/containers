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
wget http://stash.osgconnect.net/spt/public/dyang18/spt3g.tgz
tar xzf spt3g.tgz -C $PWD/code

# Setting up 3G Software environment
export SPT3G_SOFTWARE_BUILD_PATH=$PWD/code
export PATH=$PWD/code/bin:$PATH
export LD_LIBRARY_PATH=$PWD/code/spt3g:$LD_LIBRARY_PATH
export PYTHONPATH=$PWD/code:$PYTHONPATH
#$PWD/code/env-shell.sh

#SPT3G_BUILD_ROOT=$(cd `dirname $0`; pwd)
#export PATH=${SPT3G_BUILD_ROOT}/bin:$PATH
#export LD_LIBRARY_PATH=${SPT3G_BUILD_ROOT}/spt3g:$LD_LIBRARY_PATH
#export PYTHONPATH=${SPT3G_BUILD_ROOT}:$PYTHONPATH

echo "After env-shell.sh"
echo "PATH"
echo $PATH
echo "LD_LIBRARY_PATH"
echo $LD_LIBRARY_PATH
echo "PYTHONPATH:"
echo $PYTHONPATH

# Setting environment variable for `globus-url-copy`
export X509_USER_PROXY=${start_dir}/user.cert

# Making a new workspace so we don't
# accidentally delete our scripts
mkdir $PWD/tmp
work_dir=$PWD/tmp
cd ${work_dir}
mkdir ${work_dir}/output

# Transferring input file
# `$1` is the cmd line argument that is the fully qualified
# URI of the file
# `file://${work_dir}/` is where we are putting the file
globus-url-copy -vb gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/spt/data/bolodata/downsampled/RCW38-pixelraster/20950941/$1  file://${work_dir}/$1
globus-url-copy -vb gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/spt/user/production/calibration/calframe/RCW38-pixelraster/20950941.g3  file://${work_dir}/$3
# Running our processing script
# `input.file` is the basename of our input file, i.e. if the input file is
# gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/spt/data/bolodata/downsampled/calibrator/3009897/0000.g3
# the basename is 0000.g3
python ${start_dir}/makecoadd.py ${work_dir}/$3 ${work_dir}/$1  -o ${work_dir}/output/$2 -s RCW38

# Transferring output file back to storage
globus-url-copy -vb file://${work_dir}/output/$2 gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/spt/user/dyang18/output/$2

# Being nice and cleaning up after outselves
rm -rf $work_dir
