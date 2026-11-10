if ( NOT DEFINED ENABLE_HDF5 OR ENABLE_HDF5 STREQUAL "" OR ENABLE_HDF5 MATCHES " +" )
  # ENABLE_HDF5 not defined, defined as an empty string or as a string with
  # whitespaces only: variable is not provided by user
  set(HDF5_PREFER_PARALLEL TRUE)

  find_package(HDF5 COMPONENTS C QUIET) # We don't use CXX API of hdf5
  if ( HDF5_FOUND AND NOT HDF5_IS_PARALLEL )
    unset ( HDF5_FOUND )
  endif( )

else ( )
  if ( ENABLE_HDF5 )
    # User wants to use HDF5 I/O
    set(HDF5_PREFER_PARALLEL TRUE) # Only used if CMake Find modules are used

    find_package(HDF5 COMPONENTS C) # We don't use CXX API of hdf5

    if ( NOT HDF5_FOUND )
      message ( FATAL_ERROR "HDF5 library not found."
        " Please set the CMake variable ENABLE_HDF5 to OFF to disable HDF5 I/Os.")
    endif ( )

    if ( NOT HDF5_IS_PARALLEL )
      message ( FATAL_ERROR "HDF5 found but without parallel support: "
        "Please enable parallel support for HDF5 or set the CMake variable "
        "ENABLE_HDF5 to OFF to disable HDF5 I/Os.")
    endif( )

  endif ( )
endif ( )

if ( HDF5_FOUND AND NOT ENABLE_HDF5 MATCHES OFF )
  message ( STATUS "Found HDF5: add HDF5 I/O." )
  set ( ECOGEN_optional_deps "hdf5::hdf5" ) # available if HDF5 and C component are found
  set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_HDF5" )
  set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DECOGEN_SAVE_CARTESIAN_VTKHDF_AS_UNS" ) # For now disable cartesian outputs as it is inefficient

else ( )
  list ( REMOVE_ITEM ECOGEN_source_files
    ${ECOGEN_source_dir}/InputOutput/OutputVTKHDF.cpp
    ${ECOGEN_source_dir}/InputOutput/OutputVTKHDF.h )
endif ( )
