if(BUILD_TOOLS)
if(USE_RPI AND RPI_FOUND)

project(animaEddyCurrentCorrection)

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
  AnimaOptimizers
  ${ITKIO_LIBRARIES}
  ${ITK_TRANSFORM_LIBRARIES}
  ITKOptimizers
  ITKSmoothing
  ITKStatistics
  )

## #############################################################################
## install
## #############################################################################

set_exe_install_rules(${PROJECT_NAME})

endif()
endif()
