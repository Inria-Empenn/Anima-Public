if(BUILD_TOOLS)

if (USE_VTK AND VTK_FOUND)

project(animaFibersFDRCorrectPValues)

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
  AnimaStatisticalTests
  AnimaDataIO
  ${VTK_PREFIX}CommonCore
  ${VTK_PREFIX}CommonDataModel
  )

## #############################################################################
## install
## #############################################################################

set_exe_install_rules(${PROJECT_NAME})

endif()

endif()
