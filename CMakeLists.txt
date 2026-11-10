cmake_minimum_required(VERSION 3.19) # Required for target hdf5::hdf5

set (CMAKE_RELEASE_VERSION_MAJOR "6" )
set (CMAKE_RELEASE_VERSION_MINOR "0" )
set (CMAKE_RELEASE_VERSION_PATCH "0" )
set (CMAKE_RELEASE_DATE "2026-11-10" )

# Set the project name and version
project (ECOGEN
        LANGUAGES C CXX
        VERSION "${CMAKE_RELEASE_VERSION_MAJOR}.${CMAKE_RELEASE_VERSION_MINOR}.${CMAKE_RELEASE_VERSION_PATCH}")


# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Assign default compiler flags to always use
set(CMAKE_C_FLAGS "-Wall -Wextra -Wpedantic -Wundef -Wshadow -Wuninitialized ${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic -Wundef -Wshadow -Wuninitialized ${CMAKE_CXX_FLAGS}")

# Add custom compilation type
set(CMAKE_CXX_FLAGS_RELWITHASSERT "-O3" CACHE STRING
    "Flags used by the CXX compiler during RelWithAssert builds."
    FORCE)
set(CMAKE_C_FLAGS_RELWITHASSERT "-O3" CACHE STRING
    "Flags used by the C compiler during RelWithAssert builds."
    FORCE)
set(CMAKE_CXX_FLAGS_COVERAGE "-O3 --coverage" CACHE STRING
    "Flags used by the CXX compiler during Coverage builds."
    FORCE)
set(CMAKE_C_FLAGS_COVERAGE "-O3 --coverage" CACHE STRING
    "Flags used by the C compiler during Coverage builds."
    FORCE)
set(CMAKE_CXX_FLAGS_PROFILE "-O3 -pg" CACHE STRING
    "Flags used by the CXX compiler during Profile builds."
    FORCE)
set(CMAKE_C_FLAGS_PROFILE "-O3 -pg" CACHE STRING
    "Flags used by the C compiler during Profile builds."
    FORCE)

# List library source files
set(ECOGEN_source_dir ${PROJECT_SOURCE_DIR}/src)

file(
  GLOB_RECURSE
  ECOGEN_source_files
  ${ECOGEN_source_dir}/*.cpp
  )

# Remove main file from library sources
list(REMOVE_ITEM ECOGEN_source_files ${ECOGEN_source_dir}/main.cpp )

# Search for MPI distrib
set(MPI_CXX_SKIP_MPICXX ON CACHE BOOL "" FORCE)
find_package( MPI REQUIRED C) # we don't use CXX API of MPI

# Add a user option allowing to enable/disable hdf5 I/O
# (by default, hdf5 is added if hdf5-mpi is automatically detected)
set ( ENABLE_HDF5 "" CACHE STRING "Enable HDF5 I/O (ON, OFF or <empty>)" )
set_property(CACHE ENABLE_HDF5 PROPERTY STRINGS "ON" "OFF" "")

include(cmake/hdf5.cmake)

# Add ecogen library
set(CMAKE_SHARED_LIBRARY_SUFFIX_CXX ".so") # force .so extension for shared libraries

option(BUILD_SHARED_LIBS "Build using shared libraries" ON) # build ECOGEN using shared library by default

add_library(libECOGEN ${ECOGEN_source_files})
add_library ( ECOGEN::libECOGEN ALIAS libECOGEN )

set_target_properties ( libECOGEN PROPERTIES
  OUTPUT_NAME ECOGEN # Library name (we don't want cmake guess lib namre from target name because it will creates liblibECOGEN)
  POSITION_INDEPENDENT_CODE ON
  VERSION ${CMAKE_RELEASE_VERSION_MAJOR}.${CMAKE_RELEASE_VERSION_MINOR}.${CMAKE_RELEASE_VERSION_PATCH}
  SOVERSION ${CMAKE_RELEASE_VERSION_MAJOR} )

target_link_libraries(libECOGEN MPI::MPI_C ${ECOGEN_optional_deps})
target_compile_definitions(libECOGEN PUBLIC OMPI_SKIP_MPICXX)

# Add the executable
option(ECOGEN_EXE_USE_LIB "Build ECOGEN executable by linking libECOGEN library." OFF)

if (ECOGEN_EXE_USE_LIB)
  # if available on system, set rpath to the relative position of the library compared to the
  # position of the app (allow relocation of app as soon as library is moved preserving its relative position)
  set(CMAKE_BUILD_RPATH_USE_ORIGIN ON)
  add_executable(ECOGEN ${ECOGEN_source_dir}/main.cpp)
  target_link_libraries(ECOGEN libECOGEN MPI::MPI_C ${ECOGEN_optional_deps})

else()
  # Build executable from sources
  add_executable(ECOGEN ${ECOGEN_source_dir}/main.cpp ${ECOGEN_source_files})
  target_link_libraries(ECOGEN MPI::MPI_C ${ECOGEN_optional_deps})
endif()
target_compile_definitions(ECOGEN PUBLIC OMPI_SKIP_MPICXX)
