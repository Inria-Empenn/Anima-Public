cmake_minimum_required(VERSION 3.10.2)

cmake_policy(SET CMP0048 NEW)
project(ANIMA_SUPERBUILD VERSION 4.2)

set(${PROJECT_NAME}_CMAKE_DIRS
  ${CMAKE_SOURCE_DIR}/Anima/cmake
  ${CMAKE_SOURCE_DIR}/superbuild
  )

set(CMAKE_MODULE_PATH
  ${${PROJECT_NAME}_CMAKE_DIRS}
  ${CMAKE_MODULE_PATH}
  )

set (CMAKE_CXX_STANDARD 17)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(USE_VTK "Build VTK dependencies" ON)
option(USE_RPI "Build RPI dependencies" ON)
option(USE_NLOPT "Build NLOPT dependencies" ON)
option(BUILD_ANIMA_TOOLS "Build ANIMA tools" ON)
option(BUILD_ANIMA_TESTING "Build ANIMA testing executables" OFF)
option(BUILD_ANIMA_DOCUMENTATION "Build ANIMA doxygen" OFF)

# Define default build type if needed
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Use github SSH connection or https ?

option(USE_GITHUB_SSH
  "Use by default Git SSH addresses, requires public key set on github" OFF
  )
mark_as_advanced(USE_GITHUB_SSH)

if (USE_ANIMA_PRIVATE)
  set(USE_GITHUB_SSH ON)
endif()

include(CheckEnvironment)
# Deal with external dependencies

include(ExternalProject)
include(ConfigureExternalProjectVariables)
include(EP_GeneratePatchCommand)

#Variables for dependencies and update custom, updated in external projects
set(Anima_DEPS "")
set(Update_Repositories "")

# Boost
option(USE_SYSTEM_BOOST "Use system installed Boost" OFF)

if (USE_SYSTEM_BOOST)
  find_package(Boost 1.40.0 REQUIRED)
  set(Boost_SRC_DIR ${Boost_INCLUDE_DIR})
else()
  include(Boost)
endif()

# TCLAP
option(USE_SYSTEM_TCLAP "Use system installed TCLAP" OFF)

if (USE_SYSTEM_TCLAP)
  find_package(TCLAP REQUIRED)
  set(TCLAP_SRC_DIR ${TCLAP_INCLUDE_DIR})
else()
  include(TCLAP)
endif()

# TinyXML2
option(USE_SYSTEM_TinyXML2 "Use system installed TinyXML2" OFF)

if (USE_SYSTEM_TinyXML2)
  find_package(TinyXML2 REQUIRED)
  set(TinyXML2_SRC_DIR ${TinyXML2_INCLUDE_DIR})
  set(TinyXML2_BUILD_DIR ${TinyXML2_LIBRARY_DIR})
else()
  include(TinyXML2)
endif()

# NLOPT
if (USE_NLOPT)
  option(USE_SYSTEM_NLOPT "Use system installed NLOPT" OFF)

  if (USE_SYSTEM_NLOPT)
    find_package(NLOPT REQUIRED)
    set(NLOPT_SRC_DIR ${NLOPT_INCLUDE_DIR})
    set(NLOPT_BUILD_DIR ${NLOPT_LIBRARY_DIR})
  else()
    include(NLOPT)
  endif()
endif()

# VTK
if (USE_VTK)
  option(USE_SYSTEM_VTK "Use system installed VTK" OFF)

  if (USE_SYSTEM_VTK)
    find_package(VTK REQUIRED)
    set(VTK_BUILD_DIR ${VTK_DIR})
  else()
    include(VTK)
  endif()
endif()

# ITK
option(USE_SYSTEM_ITK "Use system installed ITK" OFF)

if (USE_SYSTEM_ITK)
  find_package(ITK REQUIRED)
  set(ITK_BUILD_DIR ${ITK_DIR})
else()
  include(ITK)
endif()

# RPI
if (USE_RPI)
  option(USE_SYSTEM_RPI "Use system installed RPI" OFF)

  if (USE_SYSTEM_RPI)
    find_package(RPI REQUIRED)
    set(RPI_BUILD_DIR ${RPI_DIR})
  else()
    include(RPI)
  endif()
endif()

include(Anima)

configure_file(${CMAKE_SOURCE_DIR}/superbuild/BinariesPackaging.cmake.in
  ${CMAKE_BINARY_DIR}/BinariesPackaging.cmake
  @ONLY IMMEDIATE)

add_custom_target(pack-binaries
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/BinariesPackaging.cmake
  COMMENT "-- Generating binary packages"
  )
