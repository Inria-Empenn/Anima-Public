project(${PROJECT_PREFIX}MCMProbabilisticTractography)

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
  ${TinyXML2_LIBRARY}
  ${VTK_PREFIX}CommonCore
  ${VTK_PREFIX}IOLegacy
  ${VTK_PREFIX}IOXML
  AnimaCore
  AnimaMathsOptimizers
  AnimaTractography
)

## #############################################################################
## install
## #############################################################################

set_exe_install_rules(${PROJECT_NAME})
