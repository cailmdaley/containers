#!/bin/bash
# $1 input
# $2 output
# $3 args to script
set -e

echo "input: $1"
echo "output: $2"
echo "extra args: $3"
echo ""

# Some directories
CODE=$_CONDOR_SCRATCH_DIR/code
mkdir $CODE
OUT=$_CONDOR_SCRATCH_DIR/output
mkdir $OUT
OUTPUT=$OUT/`basename $2`
mkdir $OUTPUT
STARTDIR=$PWD
DATA=$_CONDOR_SCRATCH_DIR/data
mkdir $DATA


echo Host: `hostname`
echo OS/Arch: `uname -a`
# set up environment
if [ ! -x /cvmfs/spt.opensciencegrid.org/py3-v3/setup.sh ]; then
    echo "No py3v3" >&2
    ls /cvmfs/
    exit 10
fi
eval `/cvmfs/spt.opensciencegrid.org/py3-v3/setup.sh`
echo OS_ARCH: $OS_ARCH
# Set up proxy if needed
# export connectHTTPProxy=${connectHTTPProxy:-UNAVAILABLE}
# if [ "$connectHTTPProxy" != UNAVAILABLE ]; then
#   export http_proxy=$connectHTTPProxy
# fi
# export OSG_SQUID_LOCATION=${OSG_SQUID_LOCATION:-UNAVAILABLE}
# if [ "$OSG_SQUID_LOCATION" != UNAVAILABLE ]; then
#   export http_proxy=$OSG_SQUID_LOCATION
# fi
# echo OSG proxy:  $OSG_SQUID_LOCATION
# echo Using proxy:  $http_proxy
# echo ""
# Debug/status info
echo "" 

export X509_USER_PROXY=$STARTDIR/x509up_u39371

# Get me some software
cd $CODE
wget --quiet http://stash.ci-connect.net/spt/public/ndhuang/spt3g.tgz
# And set it up
tar xf spt3g.tgz
export SPT3G_SOFTWARE_BUILD_PATH=$CODE
export LD_LIBRARY_PATH=$SPT3G_SOFTWARE_BUILD_PATH/spt3g:$LD_LIBRARY_PATH
export PYTHONPATH=$SPT3G_SOFTWARE_BUILD_PATH:$PYTHONPATH

# Get me some datas
LOCAL_PATH=""
for file in $1; do
    url="gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/$file"
    url=`echo $url|sed s:/sptpol/:/srm/spt/:g`
    LOCAL_PATH="$LOCAL_PATH $DATA/`basename $file`"
    echo Acquiring $url...
    globus-url-copy -rst -v $url file://$DATA/`basename $file`
done
echo ""

run="python $STARTDIR/find_el_glitches.py $LOCAL_PATH -o $OUTPUT/`basename $2`"
echo "Doing $run"
$run
echo ""

# get datas back
echo "Sending $2 home"
for file in `ls $OUTPUT`; do
    echo "Decquiring $file"
    globus-url-copy -cd file://$OUTPUT/$file gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/$2
done
