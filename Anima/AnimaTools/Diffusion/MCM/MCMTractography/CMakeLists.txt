project(${PROJECT_PREFIX}MCMTractography)

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
  ${VTK_PREFIX}CommonCore
  ${VTKSYS_LIBRARY}
  AnimaCore
  AnimaMathsOptimizers
  AnimaTractography
)

## #############################################################################
## install
## #############################################################################

set_exe_install_rules(${PROJECT_NAME})
