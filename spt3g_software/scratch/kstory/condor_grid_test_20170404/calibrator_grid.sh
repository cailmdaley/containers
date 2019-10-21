#!/usr/bin/env bash

echo $HOSTNAME

echo "Hello World"

echo "I will download file $1"

eval `/cvmfs/spt.opensciencegrid.org/py2-v1/setup.sh`

export X509_USER_PROXY=$PWD/kstory_proxy


# Get and unpack spt3g_software
wget http://stash.ci-connect.net/spt/public/kstory/spt3g.tgz
tar xzf spt3g.tgz
ls -l $PWD

# Transfer file
# Note: $1 is the name of the file, specified in "arguments" in calibrator_grid.submit
globus-url-copy -vb $1 file://$PWD/

# Run code
# Note: 
#   1) the script was copied to the execute node in "transfer_input_files" in calibrator_grid.submit
#   2) `basename $1` is the stripped down name of the file that has been copied to the execute node.
echo 
python calibrator_grid.py `basename $1`

# Copy outputfile back
# Note:
#   1) permissions on amundsen are messed up, so I had to make a new 
# directory with globally writeable permissions:
# $ mkdir /spt/user/kstory/test3/
# $ chmod 777 /spt/user/kstory/test3/
#   2) filename must match output from .py script, in this case 'calibrator_grid_0404.g3'
globus-url-copy file://$PWD/calibrator_grid_0404.g3 gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/spt/user/kstory/test3/
