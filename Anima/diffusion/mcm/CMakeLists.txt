project(AnimaMCM)

## #############################################################################
## List Sources
## #############################################################################

list_source_files(${PROJECT_NAME}
  ${CMAKE_CURRENT_SOURCE_DIR}
  )


## #############################################################################
## add lib
## #############################################################################

add_library(${PROJECT_NAME}
  ${${PROJECT_NAME}_CFILES}
  )

## #############################################################################
## Link
## #############################################################################

target_link_libraries(${PROJECT_NAME}
  AnimaMCMBase
  AnimaSpecialFunctions
  AnimaOptimizers
  ITKCommon
  ${TinyXML2_LIBRARY}
  )


################################################################################
# Auto generate the export file for the libs
################################################################################

generate_export_header(${PROJECT_NAME}
  STATIC_DEFINE ${PROJECT_NAME}_BUILT_AS_STATIC
  EXPORT_FILE_NAME "${PROJECT_NAME}Export.h"
  )


## #############################################################################
## install
## #############################################################################

set_lib_install_rules(${PROJECT_NAME})
