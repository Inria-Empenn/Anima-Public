project(AnimaOptimizers)

## #############################################################################
## List Sources
## #############################################################################

list_source_files(${PROJECT_NAME}
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

if(USE_NLOPT)
  list_source_files(${PROJECT_NAME} nlopt_optimizers)
endif()

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
  ITKOptimizers
  ${NLOPT_LIBRARY}
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

## #############################################################################
## Subdirs exe directories
## #############################################################################

if (BUILD_TESTING)
  add_subdirectory(nnls_test)
  add_subdirectory(nnomp_test)
  add_subdirectory(lasso_test)
endif()
