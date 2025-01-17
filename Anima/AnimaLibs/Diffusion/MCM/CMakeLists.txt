project(AnimaDiffMCM)

## #############################################################################
## List Sources
## #############################################################################

list_source_files(${PROJECT_NAME}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/Base
  ${CMAKE_CURRENT_SOURCE_DIR}/ClusteringFilters
  ${CMAKE_CURRENT_SOURCE_DIR}/MCMMain
)

## #############################################################################
## add lib
## #############################################################################

add_library(${PROJECT_NAME}
  ${${PROJECT_NAME}_CFILES}
)

## ##############################################################################
## set target properties
## ##############################################################################

target_include_directories(${PROJECT_NAME} PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/Base
  ${CMAKE_CURRENT_SOURCE_DIR}/ClusteringFilters
  ${CMAKE_CURRENT_SOURCE_DIR}/MCMMain
)

## #############################################################################
## Link
## #############################################################################

target_link_libraries(${PROJECT_NAME} PUBLIC
  ITKCommon
  ${ITKIO_LIBRARIES}
  ITKOptimizers
  ${TinyXML2_LIBRARY}
  AnimaCore
  AnimaInterpolator
  AnimaMathsMatrixOperations
  AnimaMathsOptimizers
  AnimaMathsSpecialFunctions
  AnimaMathsSphereOperations
  AnimaMathsStatisticalDistributions
  AnimaRegistrationTransformations
)

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
