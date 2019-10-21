#!/bin/sh

cd /spt/user/nwhitehorn/sptpol/converted/ra0hdec-57.5

outdir=$1
shift

for i in *; do
	echo JOB $i makemaps.sub
        echo VARS $i InputFiles=\"/spt/user/nwhitehorn/sptpol/autoproc/calibration/calframe/ra0hdec-57.5/$i.g3 $(for j in $i/*.g3; do echo `pwd`/$j; done)\"
        echo VARS $i OutputFiles=\"$outdir/$i.g3\"
        echo VARS $i ExtraArgs=\"$@\"
	echo VARS $i JobID=\"$i\"
done

