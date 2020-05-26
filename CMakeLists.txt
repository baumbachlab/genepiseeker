###############################################################################
#                                                                             #
#   Copyright (C) 2020 by David B. Blumenthal                                 #
#                                                                             #
#   This file is part of GenEpiSeeker.                                        #
#                                                                             #
#   GenEpiSeeker is free software: you can redistribute it and/or modify it   #
#   under the terms of the GNU General Public License as published by         #
#   the Free Software Foundation, either version 3 of the License, or         #
#   (at your option) any later version.                                       #
#                                                                             #
#   GenEpiSeeker is distributed in the hope that it will be useful,           #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              #
#   GNU General Public License for more details.                              #
#                                                                             #
#   You should have received a copy of the GNU General Public License         #
#   along with GenEpiSeeker. If not, see <http://www.gnu.org/licenses/>.      #
#                                                                             #
###############################################################################

cmake_minimum_required(VERSION 2.6)
project(GenEpiSeeker)

# Determine build type.
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Building GenEpiSeeker with build type 'Release', as none was specified.")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Release or Debug build type")
else()
  message(STATUS "Building GenEpiSeeker with build type '${CMAKE_BUILD_TYPE}'.")
endif()

# Find doxygen.
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doxyfile @ONLY)
  add_custom_target(docs ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxyfile WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMENT "Generating documentation with Doxygen." VERBATIM)
endif()

# Set up the compiler.
set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra -pedantic -Wredundant-decls -Wdisabled-optimization -Wctor-dtor-privacy -Wnon-virtual-dtor -Werror=old-style-cast -Wsign-promo -Werror=return-type -Wno-unused-parameter")
if(APPLE)
  set(CMAKE_CXX_COMPILER clang++ CACHE STRING "C++ compiler" FORCE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_MACOSX_RPATH ON)

# Add include and link directories.
set(BOOST_ROOT ${CMAKE_SOURCE_DIR}/ext/boost_1_71_0)
set(EIGEN_ROOT ${CMAKE_SOURCE_DIR}/ext/eigen_3_3_7)
set(CLI11_ROOT ${CMAKE_SOURCE_DIR}/ext/cli11_1_9_0)
set(CATCH_ROOT ${CMAKE_SOURCE_DIR}/ext/catch_2_11_0)
include_directories(SYSTEM ${BOOST_ROOT} ${EIGEN_ROOT} ${CLI11_ROOT} ${CATCH_ROOT})
link_directories(${BOOST_ROOT}/stage/lib)
if(APPLE)
  include_directories(SYSTEM ${OMP_ROOT}/include)
  link_directories(${OMP_ROOT}/lib)
endif()

# Add subdirectories.
add_subdirectory(test)