#!/bin/ksh -x
. /etc/profile
module use modulefiles
. ./build.ver.wcoss2
module load wcoss2.intel
if [ ! -s build ]; then
 mkdir build
fi 
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=.. -DCMAKE_INSTALL_BINDIR=exec
make -j2
make install
