#!/bin/sh

echo Executing on host: `hostname`
echo Start time: $1
echo Stop time: $2
echo Output directory: $3
echo Input Files: $4

toolset=py3-v2

# Set SHELL so setup.sh knows to do the right thing (it gets set to tcsh sometimes)
export SHELL=sh

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

env_on_error() {
	echo 'Environment at error:'
	uname -a
	printenv
}

trap env_on_error EXIT

set -e

eval `/cvmfs/spt.opensciencegrid.org/$toolset/setup.sh`

# Stats
echo 'Operating System: ' `uname`
echo 'CPU Architecture: ' `uname -p`
echo 'Platform: ' $OS_ARCH

# Move to scratch dir, copying files as necessary if that is not already where we are
[ -e $_CONDOR_SCRATCH_DIR/convert-sptpol.py ] || cp convert-sptpol.py $_CONDOR_SCRATCH_DIR/
cd $_CONDOR_SCRATCH_DIR

# Get input files
echo 'Transferring input files to scratch directory' $_CONDOR_SCRATCH_DIR
INPUT_URLS=""
LOCAL_PATH=""
for file in $4; do
	INPUT_URLS="$INPUT_URLS gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/srm/spt/`basename $file` "
	LOCAL_PATH="$LOCAL_PATH data/`basename $file`"
done
mkdir data
for url in $INPUT_URLS; do
	echo Acquiring $url...
	globus-url-copy -rst $url file://$_CONDOR_SCRATCH_DIR/data/
done

# Download software
mkdir software
software_url=http://stash.ci-connect.net/spt/public/nwhitehorn/spt3g_software_autoproc_$OS_ARCH.tgz
echo 'Downloading software distribution' $software_url
wget --quiet -O software/spt3g_software.tgz $software_url
cd software
tar xzf spt3g_software.tgz
cd ..

# Run processing
export PYTHONPATH=`pwd`/software/:$PYTHONPATH
export LD_LIBRARY_PATH=`pwd`/software/spt3g:$LD_LIBRARY_PATH
export SPT3G_SOFTWARE_BUILD_PATH=`pwd`/software
mkdir out
python convert-sptpol.py $1 $2 out/ $LOCAL_PATH

# Send back output
for i in out/*; do
	globus-url-copy -rst -r -cd file://$PWD/$i/ gsiftp://gridftp.grid.uchicago.edu:2811/cephfs$3/`basename $i`/
done

trap - EXIT

echo Finished

