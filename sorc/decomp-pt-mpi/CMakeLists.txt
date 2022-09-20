set(exe_name decomp-ptemis-mpi)
set(lib_name decomp-ptemis-mpi_lib)
set(exe_src decomp-ptemis-mpi.f90)

add_library(${lib_name} STATIC ${exe_src})
add_executable(${exe_name} ${exe_src})

set(mod_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(${lib_name} PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})
target_include_directories(${lib_name} INTERFACE ${mod_dir})

# find external libraries
find_package(MPI REQUIRED)
find_package(NetCDF REQUIRED COMPONENTS Fortran)

target_link_libraries(
  ${lib_name}
  PUBLIC 
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran
)

target_link_libraries(${exe_name} PRIVATE ${lib_name})

install(TARGETS ${exe_name})
