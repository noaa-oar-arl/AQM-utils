help([[
This module loads libraries for building the AQM-utils on
the NOAA RDHPC machine Hera using Intel-2022.1.2
]])

whatis([===[Loads libraries needed for building the AQM-utils on Hera ]===])

prepend_path("MODULEPATH","/contrib/sutils/modulefiles")
load("sutils")

load(pathJoin("cmake", os.getenv("cmake_ver") or "3.20.1"))

prepend_path("MODULEPATH","/scratch2/NCEPDEV/nwprod/hpc-stack/libs/hpc-stack/modulefiles/stack")
load(pathJoin("hpc", os.getenv("hpc_ver") or "1.2.0"))
load(pathJoin("hpc-intel", os.getenv("hpc_intel_ver") or "2022.1.2"))
load(pathJoin("hpc-impi", os.getenv("hpc_impi_ver") or "2022.1.2"))

load("netcdf/4.7.4")
load("bacio/2.4.1")
load("w3nco/2.4.1")
load("nemsio/2.5.4")
load("w3emc/2.9.2")
load("wgrib2/2.0.8")

load("bufr/11.7.1")
load("g2/3.4.5")
load("jasper/2.0.25")
load("jpeg/9.1.0")
load("libpng/1.6.37")
load("nemsiogfs/2.5.3")
load("sigio/2.3.2")
load("zlib/1.2.11")

setenv("CMAKE_C_COMPILER","mpiicc")
setenv("CMAKE_CXX_COMPILER","mpiicpc")
setenv("CMAKE_Fortran_COMPILER","mpiifort")
setenv("CMAKE_Platform","hera.intel")
