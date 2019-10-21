#!/bin/sh
# $1 input
# $2 output directory (everything transferred back to this directory)
# $3 args to script

echo "input: $1"
echo "output directory: $2"
echo "extra args: $3"
echo ""

# Set up proxy if needed
export connectHTTPProxy=${connectHTTPProxy:-UNAVAILABLE}
if [ "$connectHTTPProxy" != UNAVAILABLE ]; then
  export http_proxy=$connectHTTPProxy
fi
export OSG_SQUID_LOCATION=${OSG_SQUID_LOCATION:-UNAVAILABLE}
if [ "$OSG_SQUID_LOCATION" != UNAVAILABLE ]; then
  export http_proxy=$OSG_SQUID_LOCATION
fi
echo OSG proxy:  $OSG_SQUID_LOCATION
echo Using proxy:  $http_proxy
echo "" 
set -e

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

# Debug/status info
echo Host: `hostname`
echo OS/Arch: `uname -a`
echo OS_ARCH: $OS_ARCH
echo "" 

# set up environment
eval `/cvmfs/spt.opensciencegrid.org/py3-v2/setup.sh`
# Get me some software
cd $CODE
wget --quiet http://stash.ci-connect.net/spt/public/ndhuang/spt3g.tgz
# And set it up
tar xf spt3g.tgz
export SPT3G_SOFTWARE_BUILD_PATH=$CODE
export LD_LIBRARY_PATH=$SPT3G_SOFTWARE_BUILD_PATH/spt3g:$LD_LIBRARY_PATH
export PYTHONPATH=$SPT3G_SOFTWARE_BUILD_PATH:$PYTHONPATH
cd $STARTDIR

# Get me some datas
LOCAL_PATH=""
for file in $1; do
    url="gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/$file"
    url=`echo $url|sed s:/sptpol/:/srm/spt/:g`
    LOCAL_PATH="$LOCAL_PATH $DATA/`basename $file`"
    echo Acquiring $url...
    gcopy="globus-url-copy -v $url file://$DATA/`basename $file`"
    $gcopy||(sleep $[( $RANDOM % 30 )]s && $gcopy)||(sleep $[( $RANDOM % 60 )]s && $gcopy)|| exit 1
done
echo ""

# Do things?
run="python stuff here $LOCAL_PATH $OUTPUT"
echo "Doing $run"
$run
echo ""

# get datas back
echo "Sending $2 home"
for file in `ls $OUTPUT`; do
    echo "Decquiring $file"
    globus-url-copy -cd file://$OUTPUT/$file gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/$2/$file
done
