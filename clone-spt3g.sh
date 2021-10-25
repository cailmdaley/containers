#!/bin/bash 
  
# can specify clone location
dir=${1:-$SRC} 
 
git -C $dir clone git@github.com:SouthPoleTelescope/spt3g_software.git 