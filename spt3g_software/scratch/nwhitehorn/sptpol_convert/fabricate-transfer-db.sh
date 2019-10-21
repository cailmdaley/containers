#!/bin/sh

echo 'source	observation	status_fullrate	status_downsampled	transfer_fullrate	transfer_downsampled'

cd /spt/user/nwhitehorn/sptpol/converted
for src in *; do (
	cd -- $src
	for obs in *; do
		echo "$src	$obs	verified		1	0"
	done
); done

