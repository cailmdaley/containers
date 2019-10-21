#!/bin/sh
if [[ -e ~/mving ]]; then
    echo "locked"
    ps aux|grep mvstuff
    exit
fi
touch ~/mving
files=`find /spt/user/ndhuang/time_check/ -mmin +30 -type f`
if [[ -n $files ]]; then
    mv $files /scratch/ndhuang/time_check
fi
find /scratch/ndhuang/time_check/ -size 0 -delete
rm ~/mving