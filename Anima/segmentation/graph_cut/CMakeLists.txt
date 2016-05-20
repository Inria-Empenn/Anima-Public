project(AnimaGraphCutSegmentation)

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
  ITKCommon
  ITKIOCSV
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
## Subdirs executables
## #############################################################################

add_subdirectory(graph_cut)
add_subdirectory(gc_strem_ms_lesions_segmentation)




if(BUILD_TOOLS)

project(animaRemoveTouchingBorder)

## #############################################################################
## List Sources
## #############################################################################

list_source_files(${PROJECT_NAME}
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

## #############################################################################
## add executable
## #############################################################################

add_executable(${PROJECT_NAME}
  ${${PROJECT_NAME}_CFILES}
  )


## #############################################################################
## Link
## #############################################################################

target_link_libraries(${PROJECT_NAME}
  ${ITKIO_LIBRARIES}
  ${ITK_TRANSFORM_LIBRARIES}
  AnimaGraphCutSegmentation
  AnimaOptimizers
  ITKOptimizers
  )

## #############################################################################
## install
## #############################################################################

set_exe_install_rules(${PROJECT_NAME})

endif()
