help([[
This module loads libraries for building the AQM-utils on
the MSU machine Orion using Intel-2022.1.2
]])

whatis([===[Loads libraries needed for building the AQM-utils on Orion ]===])

load("contrib")
load("noaatools")

load(pathJoin("cmake", os.getenv("cmake_ver") or "3.22.1"))
load(pathJoin("python", os.getenv("python_ver") or "3.9.2"))

prepend_path("MODULEPATH","/apps/contrib/NCEP/libs/hpc-stack/modulefiles/stack")
load(pathJoin("hpc", os.getenv("hpc_ver") or "1.2.0"))
load(pathJoin("hpc-intel", os.getenv("hpc_intel_ver") or "2022.1.2"))
load(pathJoin("hpc-impi", os.getenv("hpc_impi_ver") or "2022.1.2"))

load("netcdf/4.7.4")
load("bacio/2.4.1")
load("w3nco/2.4.1")
load("nemsio/2.5.4")
load("w3emc/2.9.2")
load("wgrib2/2.0.8")

setenv("CMAKE_C_COMPILER","mpiicc")
setenv("CMAKE_CXX_COMPILER","mpiicpc")
setenv("CMAKE_Fortran_COMPILER","mpiifort")
setenv("CMAKE_Platform","orion.intel")
