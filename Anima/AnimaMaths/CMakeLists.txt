project(libAnimaMaths)

## #############################################################################
## List Sources
## #############################################################################

list_source_files_recurse(${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR})
list_include_dirs(${PROJECT_NAME})

## #############################################################################
## add lib
## #############################################################################

add_library(${PROJECT_NAME} ${${PROJECT_NAME}_CFILES})

## ##############################################################################
## set target properties
## ##############################################################################

target_include_directories(${PROJECT_NAME} PUBLIC ${${PROJECT_NAME}_INCLUDE_PATHS})

## #############################################################################
## Link
## #############################################################################

target_link_libraries(${PROJECT_NAME} PUBLIC
  ITKCommon
  ITKOptimizers
  ${NLOPT_LIBRARY}
)

## ##############################################################################
## Auto generate the export file for the libs
## ##############################################################################

generate_export_header(
  ${PROJECT_NAME} STATIC_DEFINE ${PROJECT_NAME}_BUILT_AS_STATIC
  EXPORT_FILE_NAME "${CMAKE_BINARY_DIR}/exports/${PROJECT_NAME}Export.h"
)

## ##############################################################################
## install
## ##############################################################################

set_lib_install_rules(${PROJECT_NAME})
