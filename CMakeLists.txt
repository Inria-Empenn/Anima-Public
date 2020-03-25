cmake_minimum_required(VERSION 3.10.2)

cmake_policy(SET CMP0048 NEW)
project(ANIMA_SUPERBUILD VERSION 3.4)

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
option(USE_ARB "Build ARB dependencies" ON)
option(BUILD_ANIMA_TOOLS "Build ANIMA tools" ON)
option(BUILD_ANIMA_TESTING "Build ANIMA testing executables" OFF)
option(BUILD_ANIMA_DOCUMENTATION "Build ANIMA doxygen" OFF)
option(USE_ANIMA_PRIVATE "Use ANIMA private part, requires authorized access" OFF)

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

# GMP
if (USE_GMP)
  find_package(GMP REQUIRED)
  set(GMP_SRC_DIR ${GMP_INCLUDE_DIR})
  set(GMP_BUILD_DIR ${GMP_LIBRARY_DIR})
endif()

# MPFR
if (USE_MPFR)
  find_package(MPFR REQUIRED)
  set(MPFR_SRC_DIR ${MPFR_INCLUDE_DIR})
  set(MPFR_BUILD_DIR ${MPFR_LIBRARY_DIR})
endif()

# FLINT
if (USE_FLINT)
  find_package(FLINT REQUIRED)
  set(FLINT_SRC_DIR ${FLINT_INCLUDE_DIR})
  set(FLINT_BUILD_DIR ${FLINT_LIBRARY_DIR})
endif()

# ARB
if (USE_ARB)
  find_package(ARB REQUIRED)
  set(ARB_SRC_DIR ${ARB_INCLUDE_DIR})
  set(ARB_BUILD_DIR ${ARB_LIBRARY_DIR})
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

if (USE_ANIMA_PRIVATE)
  set(USE_GITHUB_SSH ON CACHE BOOL "Use by default Git SSH addresses, requires public key set on github" FORCE)
  include(Anima_Private)
endif()

configure_file(${CMAKE_SOURCE_DIR}/superbuild/BinariesPackaging.cmake.in
  ${CMAKE_BINARY_DIR}/BinariesPackaging.cmake
  @ONLY IMMEDIATE)

add_custom_target(pack-binaries
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/BinariesPackaging.cmake
  DEPENDS Anima_Private
  COMMENT "-- Generating binary packages"
  )
