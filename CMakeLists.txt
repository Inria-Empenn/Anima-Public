cmake_minimum_required(VERSION 2.8.1)
project(ANIMA-SUPERBUILD)

set(${PROJECT_NAME}_CMAKE_DIRS
  ${CMAKE_SOURCE_DIR}/superbuild
  )

set(CMAKE_MODULE_PATH
  ${${PROJECT_NAME}_CMAKE_DIRS}
  ${CMAKE_MODULE_PATH}
  )

include(ExternalProject)

include(Boost)
include(TCLAP)
include(TinyXML2)
include(NLOPT)
include(ITK)
include(VTK)
include(RPI)
include(Anima-Public)