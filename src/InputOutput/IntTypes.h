//
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-.
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| |
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | |
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  |
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)|
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_)
//      (__)              (_)      (__)     (__)     (__)
//      Official webSite: https://code-mphi.github.io/ECOGEN/
//
//  This file is part of ECOGEN.
//
//  ECOGEN is the legal property of its developers, whose names
//  are listed in the copyright file included with this source
//  distribution.
//
//  ECOGEN is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published
//  by the Free Software Foundation, either version 3 of the License,
//  or (at your option) any later version.
//
//  ECOGEN is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with ECOGEN (file LICENSE).
//  If not, see <http://www.gnu.org/licenses/>.

#ifndef INTTYPES_H
#define INTTYPES_H

#include <cstdint>

static constexpr size_t MAX_SAFE_INTEGER = 9007199254740991; //!<Maximal integer exactly convertible into double for IEEE754 norm (2^53-1)

// Defines preprocessor variables ensuring type compatibility between 'size_t', HDF5 native types, and MPI types based
// on the machine architecture.
#if SIZE_MAX == UINT32_MAX // 32-bit architecture detection
  //! Maps size_t to the appropriate 32-bit unsigned integer HDF5 memory type
  #define ECOGEN_H5_SIZE_T H5T_NATIVE_UINT32
  //! Maps size_t to the appropriate 32-bit unsigned integer MPI type
  #define ECOGEN_MPI_SIZE_T MPI_UINT32_T
#elif SIZE_MAX == UINT64_MAX // 64-bit architecture detection
  //! Maps size_t to the appropriate 64-bit unsigned integer HDF5 memory type
  #define ECOGEN_H5_SIZE_T H5T_NATIVE_UINT64
  //! Maps size_t to the appropriate 64-bit unsigned integer MPI type
  #define ECOGEN_MPI_SIZE_T MPI_UINT64_T
#else
  #error Unsupported size_t size
#endif

#endif
