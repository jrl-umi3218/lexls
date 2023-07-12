#
# Copyright 2013-2021 CNRS-UM LIRMM, CNRS-AIST JRL, INRIA
#

cmake_minimum_required(VERSION 3.1.3)

set(CMAKE_CXX_STANDARD 11)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

project(
  lexls
  LANGUAGES CXX
  VERSION 0.5.0
)

include(CTest)
enable_testing()

option(INSTALL_HTML_DOCUMENTATION "Install project HTML documentation" ON)
option(INSTALL_PDF_DOCUMENTATION "Install project PDF documentation" ON)
option(BUILD_OCTAVE_BINDINGS "Build GNU Octave bindings" OFF)
option(BUILD_MATLAB_BINDINGS "Build MATLAB bindings" OFF)

find_package(Eigen3 QUIET NO_CMAKE_PACKAGE_REGISTRY)
if(Eigen3_FOUND)
  find_package(Eigen3 REQUIRED NO_CMAKE_PACKAGE_REGISTRY)
  set(EIGEN3_DEPENDENCY "find_package(Eigen3 REQUIRED NO_CMAKE_PACKAGE_REGISTRY)")
else()
  find_package(Eigen3 MODULE REQUIRED)
  set(EIGEN3_DEPENDENCY "find_package(Eigen3 MODULE REQUIRED)")
endif()

add_subdirectory(include)
add_subdirectory(cmake)
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
add_subdirectory(doc)
add_subdirectory(interfaces/matlab-octave)
