project(AnimaMCMSimilarity)

## #############################################################################
## List Sources
## #############################################################################

list_source_files(${PROJECT_NAME}
                  ${CMAKE_CURRENT_SOURCE_DIR} )


## #############################################################################
## add executable
## #############################################################################

add_library(${PROJECT_NAME}
            ${${PROJECT_NAME}_CFILES})


## #############################################################################
## Link
## #############################################################################

target_link_libraries(${PROJECT_NAME} PUBLIC
  ${ITK_TRANSFORM_LIBRARIES}
  AnimaDiffMCM
  ITKOptimizers
  )


## ##############################################################################
## set target properties
## ##############################################################################

target_include_directories(${PROJECT_NAME}  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})


################################################################################
# Auto generate the export file for the libs
################################################################################

generate_export_header(${PROJECT_NAME}
  STATIC_DEFINE ${PROJECT_NAME}_BUILT_AS_STATIC
  EXPORT_FILE_NAME "${CMAKE_BINARY_DIR}/exports/${PROJECT_NAME}Export.h"
  )

## #############################################################################
## install
## #############################################################################

set_lib_install_rules(${PROJECT_NAME})

