#!/bin/sh

echo Executing on host: `hostname`
echo Script to run: $1
echo Output File: $2
echo Toolset: $3
script=$1
output=$2
toolset=$3
shift 3
echo Input Files: $@
input=$@
echo Arguments: $SPT_SCRIPT_ARGS

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

set -e

env_on_error() {
        echo 'Environment at error:'
        uname -a
        printenv
}

trap env_on_error EXIT

export SHELL=sh
if [ -z "$SPT_NO_CVMFS" ]; then
	unset PYTHONPATH # Avoid inheriting nonsense from the environment -- thanks BU
	eval `/cvmfs/spt.opensciencegrid.org/$toolset/setup.sh`
fi

# Stats
echo 'Operating System: ' `uname`
echo 'CPU Architecture: ' `uname -p`
echo 'Platform: ' $OS_ARCH

# Move to scratch dir, copying files as necessary if that is not already where we are
if [ -n "$X509_USER_PROXY" ]; then
	[ -e $_CONDOR_SCRATCH_DIR/`basename $X509_USER_PROXY` ] || cp $X509_USER_PROXY $_CONDOR_SCRATCH_DIR
fi
[ -e $_CONDOR_SCRATCH_DIR/$script ] || cp $script $_CONDOR_SCRATCH_DIR
cd $_CONDOR_SCRATCH_DIR

# Get input files
echo 'Transferring input files to scratch directory' $_CONDOR_SCRATCH_DIR
INPUT_URLS=""
LOCAL_PATH=""
for file in $input; do
	# Detect if the file is locally accessible; if so, use it directly.
	# If not, set it up to be staged in.
	if [ -e $file ]; then
		LOCAL_PATH="$LOCAL_PATH $file"
	else
		INPUT_URLS="$INPUT_URLS gsiftp://gridftp.grid.uchicago.edu:2811/cephfs$file "
		LOCAL_PATH="$LOCAL_PATH data/`basename $file`"
	fi
done
mkdir data
for url in $INPUT_URLS; do
	echo Acquiring $url...
	globus-url-copy -rst -rst-retries 5 $url file://$_CONDOR_SCRATCH_DIR/data/
done

# Download software
if python -c 'from spt3g import core' 2>/dev/null 1>/dev/null; then
	true
else
	mkdir software
	software_url=http://stash.ci-connect.net/spt/public/nwhitehorn/spt3g_software_autoproc_$OS_ARCH.tgz
	echo 'Downloading software distribution' $software_url
	wget -O software/spt3g_software.tgz $software_url
	cd software
	tar xzf spt3g_software.tgz
	cd ..

	export PYTHONPATH=`pwd`/software/:$PYTHONPATH
	export LD_LIBRARY_PATH=`pwd`/software/spt3g:$LD_LIBRARY_PATH
	export SPT3G_SOFTWARE_BUILD_PATH=`pwd`/software
fi

# Run processing
echo python $script $SPT_SCRIPT_ARGS -o ./`basename $output` $LOCAL_PATH
python $script $SPT_SCRIPT_ARGS -o ./`basename $output` $LOCAL_PATH

# Send back output. If no X509_USER_PROXY variable defined,
# copy directly to a local filesystem since GridFTP won't work
# anyway. Otherwise, the clear intention is to use GridFTP.
if [ -n "$X509_USER_PROXY" ]; then
	globus-url-copy -rst -rst-retries 5 -cd ./`basename $output` gsiftp://gridftp.grid.uchicago.edu:2811/cephfs$output
else
	mkdir -p `dirname $output`
	cp ./`basename $output` $output
fi

trap - EXIT

echo Finished
