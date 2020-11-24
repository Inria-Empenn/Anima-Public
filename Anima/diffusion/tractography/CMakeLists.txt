# Should not be here if VTK not found but you're never too careful
if(USE_VTK AND VTK_FOUND)

project(AnimaTractography)

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
  ${TinyXML2_LIBRARY}
  AnimaDataIO
  AnimaOptimizers
  AnimaSHTools
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

add_subdirectory(dti_probabilistic_tractography)
add_subdirectory(dti_tractography)
add_subdirectory(fibers_counter)
add_subdirectory(fibers_filterer)
add_subdirectory(odf_probabilistic_tractography)

endif() #USE_VTK AND VTK_FOUND
