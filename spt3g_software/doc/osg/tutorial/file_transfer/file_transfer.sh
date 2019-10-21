#!/usr/bin/env bash

# Again seeing where we ended up
echo $HOSTNAME

# Saying Hello again
echo "Hello World"

# Telling us what file we will download
echo "I will download file $1"

# Getting the SPT software dependencies
eval `/cvmfs/spt.opensciencegrid.org/py2-v1/setup.sh`

# Copying the file to the local directory
globus-url-copy -vb $1 file://$PWD/

# Checking if our file made it
ls -l $PWD