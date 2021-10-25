#!/bin/bash 
  
procs=${1:-1} 
 
[ ! -d build ] && mkdir build 
cd build 
cmake .. 
make -j${procs} 