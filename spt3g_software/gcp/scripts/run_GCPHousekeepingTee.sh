#!/bin/bash

# Must set the following env vars to use this script:
# CRYTHON_DIR = where crython is located
# HWM_DIR = where the hardware map repo lives
# SPT3G_SOFTWARE_DIR = where the spt3g_software repo lives

if [ -z $MANPATH ]; then
  # clustertools, including python
  eval `/software/clustertools/setup.sh`
fi

if [ -z $SPT3G_SOFTWARE_PATH ]; then
  # guess the spt3g_software path location if the env var isn't set
  export SPT3G_SOFTWARE_PATH=$HOME/spt3g_software
fi

if [ -z $CRYTHON_DIR ]; then
  # guess the crython dir location if not set
  export CRYTHON_DIR=$HOME/crython
fi

if [ -z $HWM_DIR ]; then
  # guess the hardware map repo location...
  export HWM_DIR=$HOME/hardware_maps
fi

# Add the right stuff to the path
export PYTHONPATH=$HOME/:$SPT3G_SOFTWARE_PATH/build/

# Run the tee!
echo "Running GCPHousekeepingTee using hardcoded hwm. Better check it!"
python $SPT3G_SOFTWARE_PATH/gcp/python/run_GCPHousekeepingTee.py -v --hardware_map $HWM_DIR/2018/hwm_pole_run2/hwm_squids_onecrate.yml
