#!/bin/sh
export HOME=/home/ndhuang
export SPT3G_SOFTWARE_PATH=$HOME/spt_code/spt3g_software
export SPT3G_SOFTWARE_BUILD_PATH=$HOME/spt_code/spt3g_software/build
export LD_LIBRARY_PATH=$HOME/spt_code/spt3g_software/build/bin:$LD_LIBRARY_PATH
export PYTHONPATH=$HOME/spt_code/spt3g_software/build:$PYTHONPATH
eval `/software/clustertools/py3-v1/setup.sh`

env
python /home/ndhuang/spt_code/spt3g_software/scratch/ndhuang/coadd_LR_sourcemap.py $@
