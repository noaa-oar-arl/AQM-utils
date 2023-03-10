help([[
Load environment to compile UFS_UTILS on Cheyenne using Intel
]])

cmake_ver=os.getenv("cmake_ver") or "3.22.0"
load(pathJoin("cmake", cmake_ver))

ncarenv_ver=os.getenv("ncarenv_ver") or "1.3"
load(pathJoin("ncarenv", ncarenv_ver))

intel_ver=os.getenv("intel_ver") or "2022.1"
load(pathJoin("intel", intel_ver))

mpt_ver=os.getenv("mpt_ver") or "2.25"
load(pathJoin("mpt", mpt_ver))

ncarcompilers_ver=os.getenv("ncarcompilers_ver") or "0.5.0"
load(pathJoin("ncarcompilers", ncarcompilers_ver))



prepend_path("MODULEPATH", "/glade/work/epicufsrt/GMTB/tools/intel/2022.1/hpc-stack-v1.2.0_6eb6/modulefiles/stack")

hpc_ver=os.getenv("hpc_ver") or "1.2.0"
load(pathJoin("hpc", hpc_ver))

hpc_intel_ver=os.getenv("hpc_intel_ver") or "2022.1"
load(pathJoin("hpc-intel", hpc_intel_ver))

hpc_mpt_ver=os.getenv("hpc_mpt_ver") or "2.25"
load(pathJoin("hpc-mpt", hpc_mpt_ver))


load("netcdf/4.7.4")
load("bacio/2.4.1")
load("w3nco/2.4.1")
load("nemsio/2.5.2")
load("w3emc/2.9.2")
load("wgrib2/2.0.8")

load("bufr/11.7.0")
load("g2/3.4.3")
load("jasper/2.0.25")
load("jpeg/9.1.0")
load("libpng/1.6.37")
load("nemsiogfs/2.5.3")
load("sigio/2.3.2")
load("zlib/1.2.11")

setenv("CMAKE_C_COMPILER","icc")
setenv("CMAKE_CXX_COMPILER","icpc")
setenv("CMAKE_Fortran_COMPILER","ifort")
setenv("CMAKE_Platform","cheyenne.intel")

whatis("Description: AQM-utils build environment")

