cmake_minimum_required(VERSION 2.8.1)
project(ANIMA)


## #############################################################################
## Add additional cmake modules
## #############################################################################

# On Macs, disable MACOSX_RPATH (cmake 3.x)
set(CMAKE_MACOSX_RPATH OFF)

#On windows, add definitions for M_PI
if (WIN32)
  add_definitions(-D_USE_MATH_DEFINES)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /fp:strict")
endif()

set(ANIMA_CMAKE_DOCUMENTATION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cmake/documentation)
set(${PROJECT_NAME}_CMAKE_DIRS
  ${CMAKE_SOURCE_DIR}/ # for CTestConfig file
  ${CMAKE_SOURCE_DIR}/cmake
  ${ANIMA_CMAKE_DOCUMENTATION_DIR}
  ${CMAKE_SOURCE_DIR}/cmake/module
  ${CMAKE_CURRENT_BINARY_DIR}/cmake # for auto generated ${PROJECT_NAME}Exports
  )

set(CMAKE_MODULE_PATH
  ${${PROJECT_NAME}_CMAKE_DIRS}
  ${CMAKE_MODULE_PATH}
  )

include(GenerateExportHeader)
# if GenerateExportHeader is not include explicitly here,
# it may conflict with redefinition of the module in tierce lib (i.e. vtk6)
include(list_source_files)
include(list_header_directories_to_include)
include(set_exe_install_rules)
include(set_lib_install_rules)


################################################################################
# Output directories
################################################################################

set(EXECUTABLE_OUTPUT_PATH
  ${CMAKE_BINARY_DIR}/bin
  CACHE PATH "Single output directory for building all executables."
  )

set(LIBRARY_OUTPUT_PATH
  ${CMAKE_BINARY_DIR}/lib
  CACHE PATH "Single output directory for building all libraries."
  )
  
  
################################################################################
# Setup project
################################################################################

option(BUILD_SHARED_LIBS
  "Build with shared libraries."
  ON
  )

set(${PROJECT_NAME}_LIBRARY_DIRS
  ${LIBRARY_OUTPUT_PATH}
  )

set(${PROJECT_NAME}_USE_FILE
  ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Use.cmake
  )

set(CMAKE_INCLUDE_CURRENT_DIR ON)
configure_file(${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Exports.cmake.in
    ${${PROJECT_NAME}_BINARY_DIR}/cmake/${PROJECT_NAME}Exports.cmake
    @ONLY IMMEDIATE
)

include(${PROJECT_NAME}Exports)
include(AnimaModulesSetup)


################################################################################
# Dashboard
################################################################################

enable_testing()
include(CTest)
include(CTestConfig)


################################################################################
# Find dependencies
################################################################################

include(AnimaDependencies)
include(ITKLibs)


################################################################################
# Documentation
################################################################################

# Create doxygen documentation
option(BUILD_DOCUMENTATION OFF)
if(BUILD_DOCUMENTATION)
    include(UseDoxygen)
    set(DOXYFILE_LATEX "NO")
endif()


################################################################################
# Setup use file
################################################################################

if(EXISTS ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Use.cmake.in)
  configure_file( ## Common use file
    ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Use.cmake.in
    ${${PROJECT_NAME}_BINARY_DIR}/${PROJECT_NAME}Use.cmake
    @ONLY IMMEDIATE
    )
endif(EXISTS ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Use.cmake.in)


################################################################################
# Setup Configuration files
################################################################################

if(EXISTS ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in)
  configure_file( ## Build tree configure file
    ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
    ${${PROJECT_NAME}_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    @ONLY IMMEDIATE
    )
endif(EXISTS ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in)


################################################################################
# add_subdirectories
################################################################################

if (BUILD_MODULE_MATHS)
  add_subdirectory(math-tools)
endif()

if (BUILD_MODULE_FILTERING)
  add_subdirectory(filtering)
endif()

if (BUILD_MODULE_REGISTRATION)
  add_subdirectory(registration)
endif()

if (BUILD_MODULE_QUANTITATIVE_MRI)
  add_subdirectory(quantitative-mri)
endif()

if (BUILD_MODULE_DIFFUSION)
  add_subdirectory(diffusion)
endif()

if (BUILD_MODULE_SEGMENTATION)
  add_subdirectory(segmentation)
endif()
