#!/usr/bin/env bash
eval `/software/clustertools/py3-v1/setup.sh` 2>/dev/null
/home/ndhuang/spt_code/spt3g_software/build/env-shell.sh  'python -W ignore /home/ndhuang/spt_code/spt3g_software/scratch/ndhuang/autoplot.py'
