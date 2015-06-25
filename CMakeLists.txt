cmake_minimum_required(VERSION 2.8.1)
project(ANIMA-SUPERBUILD)

set(${PROJECT_NAME}_CMAKE_DIRS
  ${CMAKE_SOURCE_DIR}/Projects/Anima-Public/cmake
  ${CMAKE_SOURCE_DIR}/superbuild
  )

set(CMAKE_MODULE_PATH
  ${${PROJECT_NAME}_CMAKE_DIRS}
  ${CMAKE_MODULE_PATH}
  )

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(USE_VTK "Build VTK dependencies" ON)
option(USE_RPI "Build RPI dependencies" ON)
option(USE_NLOPT "Build NLOPT dependencies" ON)

set(cmake_flags "")
if (UNIX AND NOT BUILD_SHARED_LIBS AND "${CMAKE_SYSTEM_PROCESSOR}" MATCHES 64)
  set (cmake_flags "${cmake_flags} -fPIC")
endif()

# Use github SSH connection or https ?

option(USE_GITHUB_SSH 
  "Use by default Git SSH addresses, requires public key set on github" ON
  )
mark_as_advanced(USE_GITHUB_SSH)

if(${USE_GITHUB_SSH})
  set(GITHUB_PREFIX git@github.com:)
else()
  set(GITHUB_PREFIX https://github.com/)
endif()

# Deal with external dependencies

include(ExternalProject)
set(Anima-Public_DEPS "")

# Boost
option(USE_SYSTEM_BOOST "Use system installed Boost" OFF)

if (USE_SYSTEM_BOOST)
  find_package(Boost 1.40.0 REQUIRED)
  set(Boost_SRC_DIR ${Boost_INCLUDE_DIR})
else()
  include(Boost)
  set (Anima-Public_DEPS "${Anima-Public_DEPS};Boost")
endif()

# TCLAP
option(USE_SYSTEM_TCLAP "Use system installed TCLAP" OFF)

if (USE_SYSTEM_TCLAP)
  find_package(TCLAP REQUIRED)
else()
  include(TCLAP)
  set (Anima-Public_DEPS "${Anima-Public_DEPS};TCLAP")
endif()

# TinyXML2
option(USE_SYSTEM_TinyXML2 "Use system installed TinyXML2" OFF)

if (USE_SYSTEM_TinyXML2)
  find_package(TinyXML2 REQUIRED)
  set(TinyXML2_SRC_DIR ${TinyXML2_INCLUDE_DIR})
  set(TinyXML2_BUILD_DIR ${TinyXML2_LIBRARY_DIR})
else()
  include(TinyXML2)
  set (Anima-Public_DEPS "${Anima-Public_DEPS};TinyXML2")
endif()

# NLOPT
if (USE_NLOPT)
  option(USE_SYSTEM_NLOPT "Use system installed NLOPT" OFF)

  if (USE_SYSTEM_NLOPT)
    find_package(NLOPT REQUIRED)
    set(NLOPT_BUILD_DIR ${NLOPT_LIBRARY_DIR})
  else()
    include(NLOPT)
    set (Anima-Public_DEPS "${Anima-Public_DEPS};NLOPT")
  endif()
endif()

# ITK
option(USE_SYSTEM_ITK "Use system installed ITK" OFF)

if (USE_SYSTEM_ITK)
  find_package(ITK REQUIRED)
  set(ITK_BUILD_DIR ${ITK_DIR})
else()
  include(ITK)
  set (Anima-Public_DEPS "${Anima-Public_DEPS};ITK")
endif()

# VTK
if (USE_VTK)
  option(USE_SYSTEM_VTK "Use system installed VTK" OFF)

  if (USE_SYSTEM_VTK)
    find_package(VTK REQUIRED)
    set(VTK_BUILD_DIR ${VTK_DIR})
  else()
    include(VTK)
    set (Anima-Public_DEPS "${Anima-Public_DEPS};VTK")
  endif()
endif()

# RPI
if (USE_RPI)
  option(USE_SYSTEM_RPI "Use system installed RPI" OFF)

  if (USE_SYSTEM_RPI)
    find_package(RPI REQUIRED)
    set(RPI_BUILD_DIR ${RPI_DIR})
  else()
    include(RPI)
    set (Anima-Public_DEPS "${Anima-Public_DEPS};RPI")
  endif()
endif()

include(Anima-Public)